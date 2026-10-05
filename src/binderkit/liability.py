"""Developability liabilities and novelty for designed binder sequences.

Everything here is a filter, not a ranking signal. A design either carries a
liability or it does not. Thresholds are defaults and are meant to be
recalibrated against whatever reference set is available rather than trusted
as published constants.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydropathy.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

HYDROPHOBIC = set("AVILMFWYC")

MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Motifs that cause heterogeneity or degradation in expressed protein.
MOTIFS = {
    "n_glycosylation": r"N[^P][ST]",
    "deamidation_ng": r"NG",
    "isomerisation_dg": r"DG",
    "fragmentation_dp": r"DP",
}


@dataclass
class LiabilityReport:
    name: str
    length: int
    n_cys: int
    net_charge_ph7: float
    fraction_hydrophobic: float
    gravy: float
    max_hydrophobic_run: int
    n_his: int
    motif_hits: dict
    surface_hydrophobic_patch: float | None
    max_identity_to_reference: float | None
    passes: bool
    failures: list

    def as_row(self) -> dict:
        row = asdict(self)
        row["motif_hits"] = ";".join(f"{k}:{v}" for k, v in self.motif_hits.items() if v)
        row["failures"] = ";".join(self.failures)
        return row


def net_charge(sequence: str, ph: float = 7.0) -> float:
    """Approximate net charge from side chain counts at the given pH.

    Histidine is treated by its Henderson-Hasselbalch fraction so that the
    charge reported at 6.5 and 7.4 actually differs, which is the point.
    """
    pka = {"D": 3.9, "E": 4.3, "C": 8.3, "Y": 10.1, "H": 6.2, "K": 10.5, "R": 12.5}
    negative = sum(
        sequence.count(aa) * (1.0 / (1.0 + 10.0 ** (pka[aa] - ph))) for aa in "DECY"
    )
    positive = sum(
        sequence.count(aa) * (1.0 / (1.0 + 10.0 ** (ph - pka[aa]))) for aa in "HKR"
    )
    return positive - negative


def gravy(sequence: str) -> float:
    return float(np.mean([HYDROPATHY[a] for a in sequence if a in HYDROPATHY]))


def max_hydrophobic_run(sequence: str) -> int:
    best = current = 0
    for residue in sequence:
        current = current + 1 if residue in HYDROPHOBIC else 0
        best = max(best, current)
    return best


def motif_counts(sequence: str) -> dict[str, int]:
    return {name: len(re.findall(pattern, sequence)) for name, pattern in MOTIFS.items()}


def surface_hydrophobic_patch(
    structure_path: str | Path,
    chain_id: str,
    sasa_cut: float = 0.25,
    radius: float = 8.0,
) -> float:
    """Largest exposed hydrophobic cluster area, square angstrom.

    Exposed hydrophobic residues are clustered by CA proximity and the summed
    SASA of the largest cluster is returned. Large patches drive aggregation and
    non-specific binding.
    """
    model = PDBParser(QUIET=True).get_structure("b", str(structure_path))[0]
    ShrakeRupley().compute(model, level="R")

    exposed = []
    for residue in model[chain_id]:
        one = THREE_TO_ONE.get(residue.get_resname())
        if one is None or one not in HYDROPHOBIC or "CA" not in residue:
            continue
        if residue.sasa / MAX_ASA[one] >= sasa_cut:
            exposed.append((residue["CA"].get_coord(), residue.sasa))
    if not exposed:
        return 0.0

    coords = np.array([c for c, _ in exposed])
    areas = np.array([a for _, a in exposed])
    adjacency = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1) <= radius

    seen: set[int] = set()
    best = 0.0
    for start in range(len(exposed)):
        if start in seen:
            continue
        stack, cluster = [start], []
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            cluster.append(node)
            stack.extend(np.where(adjacency[node])[0].tolist())
        best = max(best, float(areas[cluster].sum()))
    return round(best, 1)


def _identity_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    return aligner


def sequence_identity(a: str, b: str) -> float:
    """Percent identity over the aligned region of the shorter sequence."""
    alignments = _identity_aligner().align(a, b)
    if len(alignments) == 0:
        return 0.0
    alignment = alignments[0]
    matches = aligned = 0
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        for offset in range(a_end - a_start):
            aligned += 1
            matches += a[a_start + offset] == b[b_start + offset]
    return 100.0 * matches / aligned if aligned else 0.0


def max_identity_to_reference(sequence: str, reference_fasta: str | Path | None) -> float | None:
    """Highest identity against a reference set, for the novelty requirement.

    The reference set is whatever the de novo rule makes relevant: known binders
    of the target, prior competition entries, natural homologues. No set, no
    novelty claim.
    """
    if reference_fasta is None:
        return None
    identities = [
        sequence_identity(sequence, str(record.seq).upper())
        for record in SeqIO.parse(str(reference_fasta), "fasta")
    ]
    return round(max(identities), 1) if identities else None


def assess(
    name: str,
    sequence: str,
    structure_path: str | Path | None = None,
    chain_id: str = "A",
    reference_fasta: str | Path | None = None,
    length_range: tuple[int, int] = (40, 100),
    max_cys: int = 0,
    max_hydrophobic_run_allowed: int = 5,
    max_patch_area: float = 400.0,
    max_identity: float = 40.0,
    min_his: int = 2,
) -> LiabilityReport:
    """Run every liability check on one design and collect the failures.

    max_cys defaults to zero: free cysteines in a small de novo binder give
    disulfide scrambling and dimerisation with no compensating benefit, and the
    assay format does not need them.
    """
    sequence = sequence.upper().strip()
    failures: list[str] = []

    invalid = sorted(set(sequence) - VALID_AA)
    if invalid:
        failures.append(f"invalid_residues:{''.join(invalid)}")

    length = len(sequence)
    if not (length_range[0] <= length <= length_range[1]):
        failures.append(f"length:{length}")

    n_cys = sequence.count("C")
    if n_cys > max_cys:
        failures.append(f"cys:{n_cys}")

    n_his = sequence.count("H")
    if n_his < min_his:
        failures.append(f"his:{n_his}")

    run = max_hydrophobic_run(sequence)
    if run > max_hydrophobic_run_allowed:
        failures.append(f"hydrophobic_run:{run}")

    motifs = motif_counts(sequence)
    for motif, count in motifs.items():
        if count:
            failures.append(f"{motif}:{count}")

    patch = None
    if structure_path is not None:
        patch = surface_hydrophobic_patch(structure_path, chain_id)
        if patch > max_patch_area:
            failures.append(f"hydrophobic_patch:{patch}")

    identity = max_identity_to_reference(sequence, reference_fasta)
    if identity is not None and identity > max_identity:
        failures.append(f"identity_to_reference:{identity}")

    hydrophobic_fraction = sum(sequence.count(a) for a in HYDROPHOBIC) / max(length, 1)

    return LiabilityReport(
        name=name,
        length=length,
        n_cys=n_cys,
        net_charge_ph7=round(net_charge(sequence, 7.0), 2),
        fraction_hydrophobic=round(hydrophobic_fraction, 3),
        gravy=round(gravy(sequence), 3) if sequence else 0.0,
        max_hydrophobic_run=run,
        n_his=n_his,
        motif_hits=motifs,
        surface_hydrophobic_patch=patch,
        max_identity_to_reference=identity,
        passes=not failures,
        failures=failures,
    )
