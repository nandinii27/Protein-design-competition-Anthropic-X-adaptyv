"""
Novelty guard: reject designs too close to the known agonists.

The challenge accepts de novo designs only. Nothing in this pipeline copies a
natural sequence, but two routes can drift toward one without being told to.
A language model conditioned on the receptor has seen its natural partners
during pretraining, and a generative model given hotspots drawn from the
native complex is being pointed at the pose the native ligand occupies. Either
can land on something close enough to GLP-1 or exendin-4 that "de novo"
becomes arguable.

Cheaper to check than to argue about after the fact. The comparison is against
the same four control peptides used for the functional gate.

Threshold is a judgement call, not a rule from the organisers. 40 percent over
a short peptide is already a strong resemblance; raise it only with a reason.
"""

from __future__ import annotations

import pandas as pd

import targets
from make_submission import identity

MAX_IDENTITY_TO_NATURAL = 0.40


def closest_natural(sequence: str) -> tuple[str, float]:
    """Nearest control peptide and the identity to it."""
    best_name, best_score = "", 0.0
    for name, record in targets.CONTROL_PEPTIDES.items():
        value = identity(sequence, record["sequence"])
        if value > best_score:
            best_name, best_score = name, value
    return best_name, best_score


def annotate(frame: pd.DataFrame,
             sequence_column: str = "sequence",
             threshold: float = MAX_IDENTITY_TO_NATURAL) -> pd.DataFrame:
    """Add nearest-natural columns and a pass flag to a design table."""
    out = frame.copy()
    pairs = out[sequence_column].map(closest_natural)
    out["nearest_natural"] = [name for name, _ in pairs]
    out["identity_to_natural"] = [value for _, value in pairs]
    out["passes_novelty"] = out["identity_to_natural"] <= threshold
    return out


def report(frame: pd.DataFrame) -> None:
    flagged = frame[~frame["passes_novelty"]]
    if flagged.empty:
        print("all designs below the identity threshold")
        return
    print(f"{len(flagged)} designs exceed {MAX_IDENTITY_TO_NATURAL:.0%} identity "
          "to a known agonist or antagonist:")
    print(
        flagged[["design_id", "nearest_natural", "identity_to_natural"]]
        .sort_values("identity_to_natural", ascending=False)
        .to_string(index=False)
    )
    print("\nThese are excluded from the submission. If many are flagged, the "
          "generation was anchored too hard on the native pose and that belongs "
          "in the writeup.")
