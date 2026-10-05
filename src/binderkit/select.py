"""Pareto selection over the three objectives, then diversity picking.

The objectives are not commensurable and the stated ranking order is pH
selectivity, then cross-species binding, then affinity. Collapsing them into a
weighted sum invents exchange rates that do not exist, so selection is a Pareto
front plus a diversity pass over the front.

Submission slots are limited, and a slot spent on a near-duplicate of another
submitted design buys no additional information from the assay. Diversity is
therefore part of the objective, not a tiebreak.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .liability import sequence_identity


def pareto_mask(values: np.ndarray) -> np.ndarray:
    """Boolean mask of non-dominated rows. All columns are maximised.

    A row is dominated when another row is at least equal on every objective
    and strictly better on one.
    """
    n = values.shape[0]
    mask = np.ones(n, dtype=bool)
    for i in range(n):
        if not mask[i]:
            continue
        dominates = np.all(values >= values[i], axis=1) & np.any(values > values[i], axis=1)
        if dominates.any():
            mask[i] = False
    return mask


def pareto_rank(values: np.ndarray) -> np.ndarray:
    """Non-dominated sorting rank, 0 for the first front.

    Designs beyond the first front are still usable when the front is smaller
    than the submission budget.
    """
    ranks = np.full(values.shape[0], -1, dtype=int)
    remaining = np.arange(values.shape[0])
    front = 0
    while remaining.size:
        mask = pareto_mask(values[remaining])
        ranks[remaining[mask]] = front
        remaining = remaining[~mask]
        front += 1
    return ranks


def add_pareto_rank(
    df: pd.DataFrame,
    objectives: dict[str, str],
) -> pd.DataFrame:
    """Attach a pareto_rank column.

    objectives maps column name to direction, either "max" or "min". Minimised
    columns are negated before ranking.
    """
    matrix = []
    for column, direction in objectives.items():
        series = df[column].astype(float).to_numpy()
        matrix.append(-series if direction == "min" else series)
    values = np.column_stack(matrix)
    out = df.copy()
    out["pareto_rank"] = pareto_rank(values)
    return out


def identity_distance_matrix(sequences: list[str]) -> np.ndarray:
    """Pairwise 1 - identity/100, computed once and reused by the picker."""
    n = len(sequences)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = 1.0 - sequence_identity(sequences[i], sequences[j]) / 100.0
            distances[i, j] = distances[j, i] = d
    return distances


def diverse_pick(
    df: pd.DataFrame,
    n: int,
    sequence_column: str = "sequence",
    rank_column: str = "pareto_rank",
    tiebreak_column: str | None = None,
    min_distance: float = 0.4,
) -> pd.DataFrame:
    """Greedy max-min diversity pick, front by front.

    Seeds with the best-ranked design, then repeatedly adds the candidate whose
    minimum distance to the already-selected set is largest, preferring lower
    Pareto rank. Candidates closer than min_distance to a selected design are
    skipped unless the budget cannot otherwise be filled.
    """
    if df.empty:
        return df.copy()

    sort_columns = [rank_column] + ([tiebreak_column] if tiebreak_column else [])
    ordered = df.sort_values(sort_columns, ascending=True).reset_index(drop=True)
    sequences = ordered[sequence_column].tolist()
    distances = identity_distance_matrix(sequences)

    selected = [0]
    blocked: set[int] = set()

    while len(selected) < min(n, len(ordered)):
        best_index, best_key = None, None
        for candidate in range(len(ordered)):
            if candidate in selected or candidate in blocked:
                continue
            min_dist = distances[candidate, selected].min()
            if min_dist < min_distance:
                continue
            key = (-ordered.loc[candidate, rank_column], min_dist)
            if best_key is None or key > best_key:
                best_index, best_key = candidate, key
        if best_index is None:
            # Relax the diversity floor rather than return an underfilled set.
            relaxed = [
                (distances[c, selected].min(), c)
                for c in range(len(ordered))
                if c not in selected
            ]
            if not relaxed:
                break
            best_index = max(relaxed)[1]
        selected.append(best_index)

    picked = ordered.loc[selected].copy()
    picked["min_distance_to_selected"] = [
        round(float(distances[i, [s for s in selected if s != i]].min()), 3)
        if len(selected) > 1
        else 1.0
        for i in selected
    ]
    return picked.reset_index(drop=True)


def select_submission(
    df: pd.DataFrame,
    objectives: dict[str, str],
    n: int = 20,
    sequence_column: str = "sequence",
    hard_filters: list[str] | None = None,
    primary_objective: str | None = None,
    min_distance: float = 0.4,
) -> pd.DataFrame:
    """Apply hard filters, rank by Pareto front, then pick a diverse subset.

    hard_filters names boolean columns that must all be true. These are the
    pass or fail criteria, the geometric pH filter and the liability checks,
    which are not tradeable against the objectives.
    """
    working = df.copy()
    for column in hard_filters or []:
        working = working[working[column].astype(bool)]
    if working.empty:
        return working

    working = add_pareto_rank(working, objectives)

    tiebreak = None
    if primary_objective:
        direction = objectives[primary_objective]
        tiebreak = "_tiebreak"
        series = working[primary_objective].astype(float)
        working[tiebreak] = series if direction == "min" else -series

    picked = diverse_pick(
        working,
        n=n,
        sequence_column=sequence_column,
        tiebreak_column=tiebreak,
        min_distance=min_distance,
    )
    if tiebreak and tiebreak in picked.columns:
        picked = picked.drop(columns=[tiebreak])
    return picked
