"""
TfR1 apical-domain epitope selection for dual-species (human + mouse) binder design.

Pipeline:
  1. Global align human (P02786) and mouse (Q62351) TfR1 ectodomains.
  2. Compute per-residue SASA on the human ectodomain homodimer (11QC, chain A
     in the dimer context, so the dimer interface is correctly buried).
  3. Restrict to the apical domain and to residues outside the transferrin
     footprint.
  4. Score each surface residue by (conservation, exposure) and greedily grow a
     spatially contiguous patch.
  5. Emit hotspot residue lists and cropped single-chain target PDBs for
     human and a mouse homology-threaded model.

Residue numbering throughout is UniProt / PDB numbering of mature human TfR1.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from Bio import Align, SeqIO
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

# Lawrence et al. (1999) Science domain boundaries for human TfR1 ectodomain.
APICAL_DOMAIN = (189, 383)
PROTEASE_LIKE = [(121, 188), (384, 606)]
HELICAL_DOMAIN = (607, 760)

# Transferrin footprint on TfR1: helical domain face plus the protease-like
# residues reported as Tf contacts. Anything here is excluded so the designed
# binder does not compete with iron delivery and the receptor keeps recycling.
TF_FOOTPRINT = set(range(607, 761)) | {
    569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580,
    641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652,
    653, 654, 655, 656, 657, 658, 659, 660,
}

SASA_EXPOSED_MIN = 25.0  # A^2, absolute side-chain-inclusive residue SASA
PATCH_RADIUS = 13.0      # A, CB-CB radius for a contiguous epitope patch


@dataclass
class ResidueRecord:
    resnum: int
    human_aa: str
    mouse_aa: str
    conserved: bool
    sasa: float
    domain: str
    in_tf_footprint: bool
    cb: tuple


def domain_of(resnum: int) -> str:
    if APICAL_DOMAIN[0] <= resnum <= APICAL_DOMAIN[1]:
        return "apical"
    if HELICAL_DOMAIN[0] <= resnum <= HELICAL_DOMAIN[1]:
        return "helical"
    for lo, hi in PROTEASE_LIKE:
        if lo <= resnum <= hi:
            return "protease-like"
    return "other"


def read_fasta_seq(path: Path) -> str:
    rec = next(SeqIO.parse(str(path), "fasta"))
    return str(rec.seq)


def align_species(human_full: str, mouse_full: str):
    """Global align full-length sequences; return a human-index -> mouse-aa map.

    Human and mouse TfR1 are ~78% identical over the ectodomain, so a global
    aligner with affine gaps recovers the correct register without profiles.
    Index returned is 1-based UniProt position in the human sequence.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    aln = aligner.align(human_full, mouse_full)[0]

    mapping = {}
    h_idx = 0
    m_idx = 0
    h_row, m_row = aln[0], aln[1]
    for h_char, m_char in zip(h_row, m_row):
        if h_char != "-":
            h_idx += 1
        if m_char != "-":
            m_idx += 1
        if h_char != "-":
            mapping[h_idx] = (h_char, m_char if m_char != "-" else None, m_idx if m_char != "-" else None)
    return mapping, aln


def residue_sasa(pdb_path: Path, chain_id: str):
    """SASA computed on the full biological dimer, reported for one chain.

    Computing on the dimer matters: residues at the A/B interface are not
    available to a binder, and a monomer calculation would wrongly call them
    exposed.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("tfr1", str(pdb_path))
    model = structure[0]

    sr = ShrakeRupley()
    sr.compute(model, level="R")

    out = {}
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        try:
            aa = protein_letters_3to1[residue.get_resname()]
        except KeyError:
            continue
        cb = residue["CB"] if "CB" in residue else residue["CA"]
        out[residue.id[1]] = (aa, float(residue.sasa), tuple(float(x) for x in cb.coord))
    return out, structure


def build_records(sasa_map, species_map) -> list[ResidueRecord]:
    records = []
    for resnum, (aa, sasa, cb) in sorted(sasa_map.items()):
        entry = species_map.get(resnum)
        mouse_aa = entry[1] if entry else None
        records.append(
            ResidueRecord(
                resnum=resnum,
                human_aa=aa,
                mouse_aa=mouse_aa or "-",
                conserved=(mouse_aa == aa),
                sasa=sasa,
                domain=domain_of(resnum),
                in_tf_footprint=resnum in TF_FOOTPRINT,
                cb=cb,
            )
        )
    return records


def select_patch(records, n_hotspots=6):
    """Greedy contiguous patch of conserved, exposed apical residues.

    Seed on the best-scoring candidate, then grow by adding the candidate that
    is within PATCH_RADIUS of the current patch centroid and scores highest.
    This keeps the hotspot set on one face rather than scattering it across the
    domain, which is what makes a hotspot list usable as a design constraint.
    """
    cands = [
        r for r in records
        if r.domain == "apical"
        and not r.in_tf_footprint
        and r.sasa >= SASA_EXPOSED_MIN
        and r.conserved
    ]
    if not cands:
        raise RuntimeError("no conserved exposed apical residues found")

    def score(r):
        return r.sasa

    cands.sort(key=score, reverse=True)
    patch = [cands[0]]
    pool = cands[1:]

    while len(patch) < n_hotspots and pool:
        centroid = np.mean([r.cb for r in patch], axis=0)
        reachable = [r for r in pool if np.linalg.norm(np.array(r.cb) - centroid) <= PATCH_RADIUS]
        if not reachable:
            break
        best = max(reachable, key=score)
        patch.append(best)
        pool.remove(best)

    return sorted(patch, key=lambda r: r.resnum), cands


class _CropSelect(Select):
    def __init__(self, chain_id, lo, hi):
        self.chain_id = chain_id
        self.lo = lo
        self.hi = hi

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " " and self.lo <= residue.id[1] <= self.hi

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def write_crop(structure, chain_id, lo, hi, out_path: Path, new_chain="A"):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), _CropSelect(chain_id, lo, hi))
    if new_chain != chain_id:
        text = out_path.read_text().splitlines()
        fixed = []
        for line in text:
            if line.startswith(("ATOM", "HETATM", "TER")) and len(line) > 21:
                line = line[:21] + new_chain + line[22:]
            fixed.append(line)
        out_path.write_text("\n".join(fixed) + "\n")


def thread_mouse(structure, chain_id, lo, hi, species_map, out_path: Path):
    """Write a mouse apical-domain model by threading onto the human backbone.

    At ~78% ectodomain identity the apical backbone is conserved, so renaming
    residues and truncating to the common atom set gives a target good enough
    for hotspot-directed design. Side chains at substituted positions are
    trimmed to CB and left for the design model to rebuild, rather than being
    left as the human rotamer, which would encode the wrong chemistry at
    exactly the positions that differ between species.

    Written through Biopython's PDBIO rather than hand-formatted ATOM records.
    A hand-rolled format string is how the residue name ends up one column
    left of where the spec puts it (cols 18-20), which parsers then read as an
    altLoc plus a truncated two-letter residue name. The file still looks
    plausible in a text editor and still loads far enough to crash a design
    run several minutes in.
    """
    import copy

    aa3 = {v: k for k, v in protein_letters_3to1.items()}
    backbone = {"N", "CA", "C", "O", "CB"}

    work = copy.deepcopy(structure)
    chain = work[0][chain_id]

    substituted = []
    drop_residues = []
    for residue in list(chain):
        if residue.id[0] != " " or not (lo <= residue.id[1] <= hi):
            drop_residues.append(residue.id)
            continue
        resnum = residue.id[1]
        entry = species_map.get(resnum)
        if entry is None or entry[1] is None:
            # No mouse counterpart (e.g. the V210 deletion): drop the residue.
            drop_residues.append(residue.id)
            continue
        human_aa, mouse_aa = entry[0], entry[1]
        if human_aa == mouse_aa:
            continue
        resname = aa3.get(mouse_aa)
        if resname is None:
            drop_residues.append(residue.id)
            continue
        substituted.append((resnum, human_aa, mouse_aa))
        residue.resname = resname
        for atom in list(residue):
            name = atom.get_name()
            if name not in backbone or (mouse_aa == "G" and name == "CB"):
                residue.detach_child(atom.get_id())

    for rid in drop_residues:
        if rid in [r.id for r in chain]:
            chain.detach_child(rid)

    io = PDBIO()
    io.set_structure(work)
    io.save(str(out_path), _CropSelect(chain_id, lo, hi))

    verify_pdb(out_path)
    return substituted


def verify_pdb(path: Path) -> None:
    """Re-parse a written PDB and assert it round-trips to standard residues.

    Cheap, and it is the check that catches a column-misaligned writer before
    the file is handed to a design platform.
    """
    standard = set(protein_letters_3to1)
    bad_names, bad_altloc = set(), set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            resname = line[17:20].strip()
            altloc = line[16]
            if resname not in protein_letters_3to1:
                bad_names.add(resname)
            if altloc not in (" ", "A"):
                bad_altloc.add(altloc)
    if bad_names or bad_altloc:
        raise ValueError(
            f"{path.name} is malformed: non-standard residue names {sorted(bad_names)}, "
            f"unexpected altLoc codes {sorted(bad_altloc)} — check ATOM column alignment"
        )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("check", str(path))
    n = sum(1 for r in structure.get_residues() if r.id[0] == " ")
    if n == 0:
        raise ValueError(f"{path.name} parsed to zero standard residues")


def main(raw_dir: Path, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    human_full = read_fasta_seq(raw_dir / "P02786.fasta")
    mouse_full = read_fasta_seq(raw_dir / "Q62351_mouse.fasta")

    species_map, aln = align_species(human_full, mouse_full)
    sasa_map, structure = residue_sasa(raw_dir / "11QC.pdb", "A")
    records = build_records(sasa_map, species_map)

    patch, candidates = select_patch(records)

    lo, hi = APICAL_DOMAIN
    write_crop(structure, "A", lo, hi, out_dir / "human_TfR1_apical.pdb")
    verify_pdb(out_dir / "human_TfR1_apical.pdb")
    substituted = thread_mouse(structure, "A", lo, hi, species_map, out_dir / "mouse_TfR1_apical.pdb")
    write_crop(structure, "A", 121, 760, out_dir / "human_TfR1_ectodomain_chainA.pdb")
    verify_pdb(out_dir / "human_TfR1_ectodomain_chainA.pdb")

    apical_records = [r for r in records if r.domain == "apical"]
    apical_surface = [r for r in apical_records if r.sasa >= SASA_EXPOSED_MIN]
    identity = np.mean([r.conserved for r in records if r.mouse_aa != "-"])

    summary = {
        "ectodomain_identity_human_mouse": round(float(identity), 4),
        "apical_domain": list(APICAL_DOMAIN),
        "apical_residues_modelled": len(apical_records),
        "apical_surface_residues": len(apical_surface),
        "apical_surface_conserved": sum(r.conserved for r in apical_surface),
        "apical_surface_divergent": sum(not r.conserved for r in apical_surface),
        "hotspots": [f"A{r.resnum}" for r in patch],
        "hotspot_detail": [
            {"res": f"{r.human_aa}{r.resnum}", "sasa": round(r.sasa, 1), "mouse": r.mouse_aa}
            for r in patch
        ],
        "mouse_model_substitutions_in_apical": len(substituted),
    }
    (out_dir / "epitope_summary.json").write_text(json.dumps(summary, indent=2))

    import pandas as pd
    pd.DataFrame([asdict(r) for r in records]).to_csv(out_dir / "tfr1_residue_table.csv", index=False)

    return summary, patch, records


if __name__ == "__main__":
    import sys
    raw = Path(sys.argv[1])
    out = Path(sys.argv[2])
    s, patch, records = main(raw, out)
    print(json.dumps(s, indent=2))
