"""
Assemble the final submission from scored designs.

Selection is not simply the top N by score. Generative runs produce clusters of
near-identical sequences, and spending twenty submission slots on twenty
variants of one design wastes the experiment. This script enforces a sequence
identity ceiling within the selected set and allocates slots across generation
tracks so that each track returns interpretable information about itself.

    python make_submission.py --scored outputs/scored.csv --slots 20
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import score

# Fraction of the submission budget per track. Pure groove binders are included
# as insurance for the affinity objective only; they are not expected to
# activate the receptor, so they get the smallest share.
TRACK_ALLOCATION = {
    "boltzgen_agonist": 0.40,
    "rfpeptides": 0.30,
    "pepmlm_ecd": 0.20,
    "groove_binder": 0.10,
}

MAX_IDENTITY = 0.70


def identity(a: str, b: str) -> float:
    """Ungapped identity over the shorter sequence.

    Crude by design. Generative duplicates are near-exact, so an alignment-free
    measure catches them, and borderline cases are not worth the dependency.
    """
    if not a or not b:
        return 0.0
    shorter, longer = (a, b) if len(a) <= len(b) else (b, a)
    best = 0
    for offset in range(len(longer) - len(shorter) + 1):
        window = longer[offset:offset + len(shorter)]
        matches = sum(1 for x, y in zip(shorter, window) if x == y)
        best = max(best, matches)
    return best / len(shorter)


def select_diverse(frame: pd.DataFrame, slots: int,
                   max_identity: float = MAX_IDENTITY) -> pd.DataFrame:
    """Greedy selection down the ranking, rejecting near-duplicates."""
    chosen: list[pd.Series] = []
    for _, candidate in frame.iterrows():
        if len(chosen) >= slots:
            break
        if any(identity(candidate["sequence"], picked["sequence"]) > max_identity
               for picked in chosen):
            continue
        chosen.append(candidate)
    return pd.DataFrame(chosen)


def allocate(frame: pd.DataFrame, slots: int) -> pd.DataFrame:
    selected = []
    for track, share in TRACK_ALLOCATION.items():
        budget = max(1, round(slots * share))
        subset = frame[frame["source"].astype(str).str.startswith(track)]
        if subset.empty:
            continue
        selected.append(select_diverse(subset, budget))

    if not selected:
        return select_diverse(frame, slots)

    combined = pd.concat(selected, ignore_index=True)

    # Backfill from the global ranking if some tracks returned nothing.
    if len(combined) < slots:
        remaining = frame[~frame["design_id"].isin(combined["design_id"])]
        combined = pd.concat(
            [combined, select_diverse(remaining, slots - len(combined))],
            ignore_index=True,
        )

    return combined.head(slots)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scored", type=Path, required=True)
    parser.add_argument("--slots", type=int, required=True,
                        help="submission budget for the track entered")
    parser.add_argument("--out", type=Path, default=Path("submission.csv"))
    args = parser.parse_args()

    frame = pd.read_csv(args.scored)
    frame = score.apply_filters(frame)

    rejected = frame[~frame["passed_filters"]]
    if not rejected.empty:
        print(f"{len(rejected)} designs rejected by hard filters")
        print(rejected["filter_reasons"].value_counts().to_string())

    kept = score.rank(frame[frame["passed_filters"]])
    if kept.empty:
        raise SystemExit("no designs passed the hard filters")

    final = allocate(kept, args.slots)

    submission = final[["design_id", "sequence"]].copy()
    submission.to_csv(args.out, index=False)

    print(f"\n{len(submission)} designs written")
    print(final.groupby("source").size().to_string())
    print(f"\nmean composite score: {final['composite_score'].mean():.3f}")
    print(f"mean insertion depth: {final.get('insertion_depth', pd.Series()).mean():.3f}")


if __name__ == "__main__":
    main()
