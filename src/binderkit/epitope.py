"""Epitope definition on the target extracellular region.

The epitope must satisfy three independent constraints simultaneously:

    conserved (human vs mouse)  AND  acidic (Asp/Glu)  AND  solvent exposed

Conservation gates cross-species binding. Acidic side chains are the partners
for the protonated histidines that carry the pH switch. Exposure gates whether
a designed binder can reach the residue at all.

The output is a residue list plus a spatially contiguous patch, both keyed on
author residue numbering of the target chain in the input structure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

ACIDIC = {"ASP", "GLU"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Theoretical maximum solvent accessible area per residue, Tien et al. 2013.
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


@dataclass
class ResidueRecord:
    """One target residue with every quantity the epitope decision needs."""

    resnum: int
    icode: str
    resname: str
    one_letter: str
    rel_sasa: float
    conserved: bool
    mouse_resname: str | None
    ca_coord: np.ndarray = field(repr=False)

    @property
    def is_acidic(self) -> bool:
        return self.resname in ACIDIC

    def passes(self, sasa_cut: float) -> bool:
        return self.is_acidic and self.conserved and self.rel_sasa >= sasa_cut


def _parser_for(path: Path):
    return MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)


def load_structure(path: str | Path, structure_id: str = "target"):
    path = Path(path)
    return _parser_for(path).get_structure(structure_id, str(path))


def chain_residues(structure, chain_id: str):
    """Ordered standard amino acid residues of one chain, hetero records dropped."""
    chain = structure[0][chain_id]
    return [r for r in chain if is_aa(r, standard=True)]


def chain_sequence(residues) -> str:
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in residues)


def read_fasta(path: str | Path) -> str:
    """Single-record FASTA, returned as a bare sequence string."""
    record = next(SeqIO.parse(str(path), "fasta"))
    return str(record.seq).upper()


def _aligner(mode: str = "global") -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def align_index_map(query: str, subject: str) -> dict[int, int]:
    """Map query index -> subject index for aligned, non-gap columns."""
    alignments = _aligner().align(query, subject)
    if len(alignments) == 0:
        return {}
    alignment = alignments[0]
    mapping: dict[int, int] = {}
    for (q_start, q_end), (s_start, s_end) in zip(*alignment.aligned):
        for offset in range(q_end - q_start):
            mapping[q_start + offset] = s_start + offset
    return mapping


def relative_sasa(structure, chain_id: str) -> dict[tuple[int, str], float]:
    """Shrake-Rupley SASA on the full assembly, normalised per residue type.

    SASA is computed on the whole structure so that buried interfaces between
    chains are correctly reported as buried.
    """
    ShrakeRupley().compute(structure[0], level="R")
    out: dict[tuple[int, str], float] = {}
    for residue in chain_residues(structure, chain_id):
        one = THREE_TO_ONE.get(residue.get_resname())
        if one is None:
            continue
        _, resnum, icode = residue.get_id()
        out[(resnum, icode.strip())] = residue.sasa / MAX_ASA[one]
    return out


def build_records(
    structure_path: str | Path,
    chain_id: str,
    human_fasta: str | Path,
    mouse_fasta: str | Path,
    domain_range: tuple[int, int] | None = None,
) -> list[ResidueRecord]:
    """Annotate every target residue with exposure and human/mouse identity.

    domain_range restricts the output to a residue number window on the target
    chain, for example the domain III window of the extracellular region.
    """
    structure = load_structure(structure_path)
    residues = chain_residues(structure, chain_id)
    struct_seq = chain_sequence(residues)

    human_seq = read_fasta(human_fasta)
    mouse_seq = read_fasta(mouse_fasta)

    struct_to_human = align_index_map(struct_seq, human_seq)
    human_to_mouse = align_index_map(human_seq, mouse_seq)

    sasa = relative_sasa(structure, chain_id)

    records: list[ResidueRecord] = []
    for position, residue in enumerate(residues):
        _, resnum, icode = residue.get_id()
        icode = icode.strip()
        if domain_range is not None and not (domain_range[0] <= resnum <= domain_range[1]):
            continue

        human_index = struct_to_human.get(position)
        mouse_index = human_to_mouse.get(human_index) if human_index is not None else None

        if mouse_index is None:
            conserved, mouse_resname = False, None
        else:
            mouse_aa = mouse_seq[mouse_index]
            human_aa = human_seq[human_index]
            conserved = mouse_aa == human_aa
            mouse_resname = mouse_aa

        if "CA" not in residue:
            continue

        records.append(
            ResidueRecord(
                resnum=resnum,
                icode=icode,
                resname=residue.get_resname(),
                one_letter=THREE_TO_ONE.get(residue.get_resname(), "X"),
                rel_sasa=sasa.get((resnum, icode), 0.0),
                conserved=conserved,
                mouse_resname=mouse_resname,
                ca_coord=residue["CA"].get_coord().astype(float),
            )
        )
    return records


def anchor_residues(records: list[ResidueRecord], sasa_cut: float = 0.25) -> list[ResidueRecord]:
    """Residues satisfying conserved AND acidic AND exposed.

    These are the carboxylate partners the binder histidines must reach. If this
    list is empty for the chosen domain window, the pH mechanism has no
    structural basis there and the epitope must move before anything downstream
    is worth running.
    """
    return [r for r in records if r.passes(sasa_cut)]


def patch_around(
    records: list[ResidueRecord],
    centre: ResidueRecord,
    radius: float = 12.0,
    sasa_cut: float = 0.25,
) -> list[ResidueRecord]:
    """All exposed residues whose CA lies within radius of a centre residue.

    This is the hotspot set handed to the backbone generator, not just the
    acidic anchors: the generator needs the surrounding surface to build against.
    """
    return [
        r
        for r in records
        if r.rel_sasa >= sasa_cut
        and float(np.linalg.norm(r.ca_coord - centre.ca_coord)) <= radius
    ]


def rank_patches(
    records: list[ResidueRecord],
    radius: float = 12.0,
    sasa_cut: float = 0.25,
    min_anchors: int = 2,
) -> list[dict]:
    """Score every anchor-centred patch by how well it supports the pH switch.

    Ranking key, in order: number of conserved acidic anchors in the patch,
    then the fraction of the patch that is conserved. A patch with fewer than
    min_anchors carboxylates cannot host two independent histidine salt bridges
    and is dropped.
    """
    anchors = anchor_residues(records, sasa_cut)
    patches = []
    for centre in anchors:
        members = patch_around(records, centre, radius, sasa_cut)
        patch_anchors = [r for r in members if r.passes(sasa_cut)]
        if len(patch_anchors) < min_anchors:
            continue
        conserved_fraction = sum(r.conserved for r in members) / max(len(members), 1)
        patches.append(
            {
                "centre": f"{centre.resname}{centre.resnum}{centre.icode}",
                "centre_resnum": centre.resnum,
                "n_anchors": len(patch_anchors),
                "anchors": [f"{r.resname}{r.resnum}{r.icode}" for r in patch_anchors],
                "n_residues": len(members),
                "conserved_fraction": round(conserved_fraction, 3),
                "mean_rel_sasa": round(float(np.mean([r.rel_sasa for r in members])), 3),
                "hotspot_resnums": sorted(r.resnum for r in members),
            }
        )
    patches.sort(key=lambda p: (p["n_anchors"], p["conserved_fraction"]), reverse=True)
    return patches


def hotspot_string(patch: dict, chain_id: str) -> str:
    """Hotspot specification in the chain-prefixed form the generators take."""
    return ",".join(f"{chain_id}{n}" for n in patch["hotspot_resnums"])
