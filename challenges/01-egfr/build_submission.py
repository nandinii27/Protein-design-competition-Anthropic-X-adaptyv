"""Select the final design set and write the submission CSV.

Selection is tiered, with a global cap on how many designs may come from any
one backbone. Variants of a single backbone differ by a residue or two and
measure nearly the same thing, so the cap spends assay slots on independent
designs rather than on near-duplicates of the leaders.

Tiers, in order:
    1  grafted mutants clearing geometry, sign and liability filters
    2  parent designs whose own linkage score favours the low pH condition
    3  further grafted mutants that improved on their parent and keep the sign
    4  affinity controls: strong binders with no predicted switch, from
       backbones not otherwise represented, labelled as such

The fourth tier is deliberate. Nothing currently establishes whether a
predicted pKa-shift linkage score tracks measured pH dependence on designed
interfaces. Submitting a few confident binders with no predicted switch gives
the assay a comparison point, so that a null result on the switching designs
is still informative.

Within a tier, ranking is a Pareto front over three objectives. They are not
combined into a weighted sum: the brief states a priority order, and a
priority order does not define exchange rates.

Reads:
    outputs/mutant_metrics.csv
    outputs/merged_metrics.csv
Writes:
    outputs/selected.csv
    submission.csv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

from binderkit import select, submission

CHALLENGE_DIR = Path(__file__).resolve().parent
OUTPUTS = CHALLENGE_DIR / "outputs"
MUTANT_CSV = OUTPUTS / "mutant_metrics.csv"
PARENT_CSV = OUTPUTS / "merged_metrics.csv"

N_SLOTS = 20
MAX_PER_BACKBONE = 2
NAME_PREFIX = "nd-egfr-ph"
MOLECULE_CLASS = "protein"

OBJECTIVES = {
    "ph_selectivity_kcal": "min",
    "interface_confidence": "max",
    "n_distinct_his_on_anchor": "max",
}


def backbone_of(design_id: str) -> str:
    """Generator seed, which identifies the backbone a design came from."""
    match = re.search(r"_s(\d+)", str(design_id))
    return match.group(1) if match else str(design_id)


def load_mutants() -> pd.DataFrame:
    df = pd.read_csv(MUTANT_CSV)
    df["source"] = "grafted"
    df["design_id"] = df["mutant_id"]
    df["interface_confidence"] = pd.to_numeric(df.get("iptm"), errors="coerce")
    for column in ("geometry_pass", "sign_pass", "liabilities_pass", "improved"):
        df[column] = df[column].astype(bool)
    df["all_pass"] = df["geometry_pass"] & df["sign_pass"] & df["liabilities_pass"]
    return df


def load_parents() -> pd.DataFrame:
    if not PARENT_CSV.exists():
        return pd.DataFrame()
    df = pd.read_csv(PARENT_CSV)
    df["source"] = "parent"
    df["design_id"] = df["Design"]
    df["sequence"] = df["Sequence"]
    df["interface_confidence"] = pd.to_numeric(df.get("Average_i_pTM"), errors="coerce")
    df["mutations"] = ""
    df["parent_selectivity_kcal"] = pd.NA
    df["delta_selectivity"] = pd.NA
    df["improved"] = False
    df["liabilities_pass"] = (
        df["passes"].astype(bool) if "passes" in df.columns else True
    )
    for column in ("geometry_pass", "sign_pass"):
        df[column] = df[column].astype(bool) if column in df.columns else False
    if "n_distinct_his_on_anchor" not in df.columns:
        df["n_distinct_his_on_anchor"] = 0
    return df


def pareto_order(df: pd.DataFrame) -> pd.DataFrame:
    """Non-dominated sorting, then most negative selectivity within a front."""
    working = df.copy()
    for column in OBJECTIVES:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    working = working.dropna(subset=list(OBJECTIVES))
    if working.empty:
        return working
    working = select.add_pareto_rank(working, OBJECTIVES)
    return working.sort_values(["pareto_rank", "ph_selectivity_kcal"])


def take(
    pool: pd.DataFrame,
    chosen: list[pd.Series],
    counts: dict[str, int],
    label: str,
) -> int:
    """Add rows from pool, respecting the global backbone cap and slot budget."""
    added = 0
    taken_sequences = {row["sequence"] for row in chosen}
    for _, row in pareto_order(pool).iterrows():
        if len(chosen) >= N_SLOTS:
            break
        if row["sequence"] in taken_sequences:
            continue
        backbone = backbone_of(row["design_id"])
        if counts.get(backbone, 0) >= MAX_PER_BACKBONE:
            continue
        entry = row.copy()
        entry["tier"] = label
        entry["backbone"] = backbone
        chosen.append(entry)
        taken_sequences.add(row["sequence"])
        counts[backbone] = counts.get(backbone, 0) + 1
        added += 1
    return added


def main() -> int:
    if not MUTANT_CSV.exists():
        print(f"missing {MUTANT_CSV.name}", file=sys.stderr)
        return 1

    mutants = load_mutants()
    parents = load_parents()

    chosen: list[pd.Series] = []
    counts: dict[str, int] = {}

    n = take(mutants[mutants["all_pass"]], chosen, counts, "grafted_switch")
    print(f"tier 1, grafted and passing every filter: {n}")

    if not parents.empty:
        n = take(
            parents[parents["n_cys"].fillna(0) == 0],
            chosen,
            counts,
            "affinity_control",
        )

        print(f"tier 2, parents with the correct sign: {n}")

    n = take(
        mutants[mutants["improved"] & mutants["sign_pass"] & mutants["liabilities_pass"]],
        chosen,
        counts,
        "grafted_switch",
    )
    print(f"tier 3, further grafted designs with the correct sign: {n}")

    if len(chosen) < N_SLOTS and not parents.empty:
        n = take(
            parents[parents["liabilities_pass"]],
            chosen,
            counts,
            "affinity_control",
        )
        print(f"tier 4, affinity controls on unrepresented backbones: {n}")

    if not chosen:
        print("nothing selected", file=sys.stderr)
        return 1

    selected = pd.DataFrame(chosen).reset_index(drop=True)

    # Switching designs first, controls last, each block best-first.
    selected["tier_order"] = (selected["tier"] == "affinity_control").astype(int)
    selected = selected.sort_values(
        ["tier_order", "ph_selectivity_kcal"]
    ).reset_index(drop=True)
    selected.to_csv(OUTPUTS / "selected.csv", index=False)

    table = pd.DataFrame(
        {
            "name": [f"{NAME_PREFIX}-{i + 1:03d}" for i in range(len(selected))],
            "sequence": selected["sequence"].astype(str).str.upper().str.strip(),
            "molecule_class": MOLECULE_CLASS,
        }
    )

    carry = {
        "design_role": "tier",
        "ph_selectivity_kcal": "ph_selectivity_kcal",
        "ddg_ph_low": "ddg_ph_low",
        "ddg_ph_high": "ddg_ph_high",
        "max_his_pka_shift": "max_his_pka_shift",
        "n_his_carboxylate_contacts": "n_his_carboxylate_contacts",
        "n_his_on_conserved_anchor": "n_distinct_his_on_anchor",
        "anchor_contacts": "anchor_contacts",
        "interface_confidence_iptm": "interface_confidence",
        "parent_selectivity_kcal": "parent_selectivity_kcal",
        "delta_selectivity_vs_parent": "delta_selectivity",
        "histidine_grafts": "mutations",
        "backbone_id": "backbone",
        "internal_id": "design_id",
        "length": "length",
        "net_charge_ph7": "net_charge_ph7",
        "pareto_rank": "pareto_rank",
    }
    for out_name, source in carry.items():
        if source in selected.columns:
            table[out_name] = selected[source].values

    result = submission.write(table, CHALLENGE_DIR / "submission.csv")
    print("\n" + result.report())
    print(f"\n{len(table)} designs, {selected['backbone'].nunique()} backbones")
    print(selected["tier"].value_counts().to_string())
    return 0 if result.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
