"""
Scoring and filtering of candidate designs.

Two layers, kept separate on purpose.

Hard filters are the constraints that disqualify a design regardless of how
good it looks: length bounds from the challenge rules, hydrophobicity (the
receptor is in detergent micelles, so a design with large exposed hydrophobic
surface partitions into detergent and reads as aggregation rather than
binding), and free cysteines that can scramble against the native disulfides.

The composite score is a weighted combination of co-folding confidence and
pocket insertion. The weights are a judgement call, not a calibrated model -
the point of benchmark_controls.py is to find out whether the combination
separates known agonists from known antagonists before any of it is trusted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import targets

# Kyte-Doolittle hydropathy.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

VALID_RESIDUES = set(HYDROPATHY)

# Composite weights. Pocket recovery carries the most weight because it is the
# only term that bears on the activation objective; co-folding confidence is
# necessary but, on this target class, not sufficient.
WEIGHTS = {
    "pocket_contact_recovery": 0.40,
    "insertion_depth": 0.25,
    "iptm": 0.25,
    "groove_contact_recovery": 0.10,
}


@dataclass
class FilterResult:
    passed: bool
    reasons: list[str]


def gravy(sequence: str) -> float:
    """Mean hydropathy. Positive means net hydrophobic."""
    values = [HYDROPATHY[residue] for residue in sequence if residue in HYDROPATHY]
    return float(np.mean(values)) if values else 0.0


def max_hydrophobic_run(sequence: str) -> int:
    """Longest consecutive stretch of strongly hydrophobic residues.

    A long run is a stronger detergent-partitioning signal than mean
    hydropathy, because it is what forms a membrane-insertable segment.
    """
    longest = current = 0
    for residue in sequence:
        if HYDROPATHY.get(residue, 0.0) > 1.5:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def hard_filters(sequence: str) -> FilterResult:
    reasons = []

    if not sequence or set(sequence) - VALID_RESIDUES:
        reasons.append("non-standard residues")

    if not targets.MIN_LENGTH <= len(sequence) <= targets.MAX_LENGTH:
        reasons.append(f"length {len(sequence)} outside challenge bounds")

    if gravy(sequence) > targets.MAX_GRAVY:
        reasons.append(f"GRAVY {gravy(sequence):.2f} too hydrophobic for detergent assay")

    if max_hydrophobic_run(sequence) >= 8:
        reasons.append("long hydrophobic run, membrane-partitioning risk")

    if sequence.count("C") > targets.MAX_FREE_CYSTEINES:
        reasons.append("free cysteine")

    return FilterResult(passed=not reasons, reasons=reasons)


def composite_score(row: pd.Series) -> float:
    """Weighted score over available metric columns.

    Missing metrics are treated as zero rather than dropped, so a design that
    has not been through the geometry stage never outranks one that has.
    """
    total = 0.0
    for column, weight in WEIGHTS.items():
        value = row.get(column, 0.0)
        if pd.isna(value):
            value = 0.0
        total += weight * float(value)
    return total


def apply_filters(frame: pd.DataFrame, sequence_column: str = "sequence") -> pd.DataFrame:
    """Annotate a design table with filter outcome and sequence properties."""
    out = frame.copy()
    results = out[sequence_column].map(hard_filters)
    out["passed_filters"] = [result.passed for result in results]
    out["filter_reasons"] = ["; ".join(result.reasons) for result in results]
    out["gravy"] = out[sequence_column].map(gravy)
    out["max_hydrophobic_run"] = out[sequence_column].map(max_hydrophobic_run)
    return out


def rank(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["composite_score"] = out.apply(composite_score, axis=1)
    return out.sort_values("composite_score", ascending=False).reset_index(drop=True)
