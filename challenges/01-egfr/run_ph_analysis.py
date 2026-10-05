"""Score generated designs for pH selectivity, liabilities and anchor contacts.

Expects, relative to the challenge directory:
    final_design_stats.csv              generator output, one row per design
    outputs/complexes/<Design>.pdb      predicted complex, target chain A, binder chain B

Writes:
    outputs/ph_scores.csv               per-design linkage score and contacts
    outputs/liabilities.csv             per-design developability report
    outputs/merged_metrics.csv          everything joined, ready for selection

Run from the challenge directory with the project environment active.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from binderkit import liability, phscore

CHALLENGE_DIR = Path(__file__).resolve().parent
OUTPUTS = CHALLENGE_DIR / "outputs"
COMPLEX_DIR = OUTPUTS / "complexes"
STATS_CSV = CHALLENGE_DIR / "final_design_stats.csv"

TARGET_CHAIN = "A"
BINDER_CHAIN = "B"
PH_LOW, PH_HIGH = 6.5, 7.4

# Conserved, exposed carboxylates defining the epitope.
ANCHORS = (21, 24, 65, 68)

# Liability thresholds. Relaxed relative to the module defaults because the
# generator was run with relaxed filters and the pool is small.
LENGTH_RANGE = (60, 110)
MAX_PATCH_AREA = 600.0
MIN_HIS = 1


def load_designs() -> pd.DataFrame:
    df = pd.read_csv(STATS_CSV)
    df = df[df["Design"].notna()].copy()
    df["Design"] = df["Design"].astype(str).str.strip()
    df["Sequence"] = df["Sequence"].astype(str).str.strip().str.upper()
    return df


def complex_path(design: str) -> Path | None:
    direct = COMPLEX_DIR / f"{design}.pdb"
    if direct.exists():
        return direct
    matches = sorted(COMPLEX_DIR.glob(f"{design}*.pdb"))
    return matches[0] if matches else None


def anchor_contacts(contacts: list[dict]) -> list[dict]:
    """Restrict histidine contacts to the four epitope carboxylates."""
    return [
        c
        for c in contacts
        if any(str(a) in c["target_residue"] for a in ANCHORS)
    ]


def main() -> int:
    if not STATS_CSV.exists():
        print(f"missing {STATS_CSV.name}", file=sys.stderr)
        return 1
    if not COMPLEX_DIR.exists():
        print(f"missing complex directory {COMPLEX_DIR}", file=sys.stderr)
        return 1

    designs = load_designs()
    print(f"{len(designs)} designs in stats file")

    ph_rows, liability_rows, missing = [], [], []

    for _, row in designs.iterrows():
        design = row["Design"]
        path = complex_path(design)
        if path is None:
            missing.append(design)
            continue

        try:
            result = phscore.score_design(
                complex_path=path,
                binder_chain=BINDER_CHAIN,
                target_chain=TARGET_CHAIN,
                ph_low=PH_LOW,
                ph_high=PH_HIGH,
                min_his_contacts=1,
                workdir=OUTPUTS / "propka" / design,
            )
        except Exception as exc:
            print(f"{design}: pH scoring failed, {exc}", file=sys.stderr)
            continue

        on_anchor = anchor_contacts(result["contacts"])
        flat = {k: v for k, v in result.items() if k not in ("per_site", "contacts")}
        flat["Design"] = design
        flat["n_anchor_his_contacts"] = len(on_anchor)
        flat["n_distinct_his_on_anchor"] = len({c["his_resnum"] for c in on_anchor})
        flat["anchor_contacts"] = ";".join(
            f"H{c['his_resnum']}-{c['target_residue']}@{c['distance']}" for c in on_anchor
        )
        ph_rows.append(flat)

        report = liability.assess(
            name=design,
            sequence=row["Sequence"],
            structure_path=path,
            chain_id=BINDER_CHAIN,
            reference_fasta=None,
            length_range=LENGTH_RANGE,
            max_patch_area=MAX_PATCH_AREA,
            min_his=MIN_HIS,
        )
        liability_rows.append({**report.as_row(), "Design": design})

        print(
            f"{design}: selectivity {flat['ph_selectivity_kcal']:+.3f} kcal/mol, "
            f"{flat['n_distinct_his_on_anchor']} His on anchors, "
            f"liabilities {'pass' if report.passes else report.failures}"
        )

    if missing:
        print(f"\nno complex found for: {missing}", file=sys.stderr)

    if not ph_rows:
        print("no designs scored", file=sys.stderr)
        return 1

    ph_df = pd.DataFrame(ph_rows)
    liability_df = pd.DataFrame(liability_rows).drop(columns=["name"], errors="ignore")

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    ph_df.to_csv(OUTPUTS / "ph_scores.csv", index=False)
    liability_df.to_csv(OUTPUTS / "liabilities.csv", index=False)

    keep = [
        "Design",
        "Sequence",
        "Length",
        "Average_i_pTM",
        "Average_i_pAE",
        "Average_pLDDT",
        "Average_ShapeComplementarity",
        "Average_dG",
        "Average_dSASA",
        "Average_n_InterfaceResidues",
        "Average_Hotspot_RMSD",
        "Average_Binder_RMSD",
    ]
    keep = [c for c in keep if c in designs.columns]

    merged = (
        designs[keep]
        .merge(ph_df, on="Design", how="inner")
        .merge(liability_df, on="Design", how="inner", suffixes=("", "_liab"))
    )
    merged.to_csv(OUTPUTS / "merged_metrics.csv", index=False)

    print(f"\nscored {len(merged)} designs")
    print(f"geometry pass (His on any carboxylate): {int(merged['geometry_pass'].sum())}")
    print(f"sign pass (favours low pH): {int(merged['sign_pass'].sum())}")
    print(f"His contacting a chosen anchor: {int((merged['n_distinct_his_on_anchor'] > 0).sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
