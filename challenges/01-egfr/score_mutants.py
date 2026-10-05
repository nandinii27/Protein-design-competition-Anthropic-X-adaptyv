"""Score re-predicted histidine mutants and compare them to their parents.

Reads:
    outputs/mutants.csv              mutant sequences and their mutations
    outputs/colabfold/**/*rank_001*.pdb   re-predicted complexes
    outputs/colabfold/**/*scores_rank_001*.json
    outputs/ph_scores.csv            parent linkage scores, for the delta

Writes:
    outputs/mutant_metrics.csv       every mutant, scored and compared

The delta against the parent is the quantity of interest: an absolute
selectivity score mixes the grafted histidine with whatever incidental
protonation effects the parent already had. The change isolates the graft.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pandas as pd

from binderkit import liability, phscore

CHALLENGE_DIR = Path(__file__).resolve().parent
OUTPUTS = CHALLENGE_DIR / "outputs"
FOLD_DIR = OUTPUTS / "colabfold"
MUTANTS_CSV = OUTPUTS / "mutants.csv"
PARENT_PH_CSV = OUTPUTS / "ph_scores.csv"

TARGET_CHAIN = "A"
BINDER_CHAIN = "B"
PH_LOW, PH_HIGH = 6.5, 7.4

ANCHORS = (21, 24, 65, 68)

LENGTH_RANGE = (60, 110)
MAX_PATCH_AREA = 600.0
MIN_HIS = 1


def find_structure(mutant_id: str) -> Path | None:
    matches = sorted(FOLD_DIR.rglob(f"{mutant_id}_unrelaxed_rank_001*.pdb"))
    if not matches:
        matches = sorted(FOLD_DIR.rglob(f"{mutant_id}*rank_001*.pdb"))
    return matches[0] if matches else None


def find_scores(mutant_id: str) -> dict:
    matches = sorted(FOLD_DIR.rglob(f"{mutant_id}*scores_rank_001*.json"))
    if not matches:
        return {}
    data = json.loads(matches[0].read_text())
    plddt = data.get("plddt") or []
    return {
        "iptm": data.get("iptm"),
        "ptm": data.get("ptm"),
        "mean_plddt": round(sum(plddt) / len(plddt), 2) if plddt else None,
        "pdockq": data.get("pdockq"),
    }


def anchor_contacts(contacts: list[dict]) -> list[dict]:
    return [c for c in contacts if any(str(a) in c["target_residue"] for a in ANCHORS)]


def main() -> int:
    if not MUTANTS_CSV.exists():
        print(f"missing {MUTANTS_CSV.name}", file=sys.stderr)
        return 1
    if not FOLD_DIR.exists():
        print(f"missing {FOLD_DIR}", file=sys.stderr)
        return 1

    mutants = pd.read_csv(MUTANTS_CSV)

    parents = {}
    if PARENT_PH_CSV.exists():
        parent_df = pd.read_csv(PARENT_PH_CSV)
        parents = dict(
            zip(parent_df["Design"], parent_df["ph_selectivity_kcal"].astype(float))
        )

    rows, missing = [], []

    # Structures are staged under short names before scoring. The predictor
    # filenames plus a per-mutant directory exceed the Windows path limit,
    # which fails as a missing-file error rather than anything informative.
    scratch = OUTPUTS / "pk"
    if scratch.exists():
        shutil.rmtree(scratch, ignore_errors=True)
    scratch.mkdir(parents=True, exist_ok=True)

    for index, (_, row) in enumerate(mutants.iterrows()):
        mutant_id = str(row["mutant_id"])
        path = find_structure(mutant_id)
        if path is None:
            missing.append(mutant_id)
            continue

        workdir = scratch / f"m{index:03d}"
        workdir.mkdir(parents=True, exist_ok=True)
        staged = workdir / "c.pdb"
        shutil.copy(path, staged)

        try:
            result = phscore.score_design(
                complex_path=staged,
                binder_chain=BINDER_CHAIN,
                target_chain=TARGET_CHAIN,
                ph_low=PH_LOW,
                ph_high=PH_HIGH,
                min_his_contacts=1,
                workdir=workdir,
            )
        except Exception as exc:
            print(f"{mutant_id}: pH scoring failed, {exc}", file=sys.stderr)
            continue

        on_anchor = anchor_contacts(result["contacts"])
        report = liability.assess(
            name=mutant_id,
            sequence=str(row["sequence"]),
            structure_path=staged,
            chain_id=BINDER_CHAIN,
            reference_fasta=None,
            length_range=LENGTH_RANGE,
            max_patch_area=MAX_PATCH_AREA,
            min_his=MIN_HIS,
        )

        parent_score = parents.get(row["Design"])
        selectivity = result["ph_selectivity_kcal"]
        delta = round(selectivity - parent_score, 3) if parent_score is not None else None

        entry = {
            "mutant_id": mutant_id,
            "Design": row["Design"],
            "mutations": row["mutations"],
            "anchors_targeted": row["anchors"],
            "n_mutations": row["n_mutations"],
            "sequence": row["sequence"],
            "ph_selectivity_kcal": selectivity,
            "parent_selectivity_kcal": parent_score,
            "delta_selectivity": delta,
            "improved": (delta is not None and delta < 0),
            "ddg_ph_low": result["ddg_ph_low"],
            "ddg_ph_high": result["ddg_ph_high"],
            "max_his_pka_shift": result["max_his_pka_shift"],
            "n_his_carboxylate_contacts": result["n_his_carboxylate_contacts"],
            "n_distinct_his_on_anchor": len({c["his_resnum"] for c in on_anchor}),
            "anchor_contacts": ";".join(
                f"H{c['his_resnum']}-{c['target_residue']}@{c['distance']}" for c in on_anchor
            ),
            "geometry_pass": result["geometry_pass"],
            "sign_pass": result["sign_pass"],
            "liabilities_pass": report.passes,
            "liability_failures": ";".join(report.failures),
            "length": report.length,
            "n_cys": report.n_cys,
            "net_charge_ph7": report.net_charge_ph7,
            "surface_hydrophobic_patch": report.surface_hydrophobic_patch,
        }
        entry.update(find_scores(mutant_id))
        rows.append(entry)

        print(
            f"{mutant_id}: {selectivity:+.3f} kcal/mol "
            f"(parent {parent_score:+.3f}, delta {delta:+.3f}), "
            f"{entry['n_distinct_his_on_anchor']} His on anchors, "
            f"ipTM {entry.get('iptm')}"
            if parent_score is not None
            else f"{mutant_id}: {selectivity:+.3f} kcal/mol"
        )

    if missing:
        print(f"\nno structure found for: {missing}", file=sys.stderr)
    if not rows:
        print("nothing scored", file=sys.stderr)
        return 1

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUTS / "mutant_metrics.csv", index=False)

    print(f"\nscored {len(df)} mutants")
    print(f"improved on parent: {int(df['improved'].sum())}")
    print(f"sign pass (favours low pH): {int(df['sign_pass'].sum())}")
    print(f"His contacting an anchor: {int((df['n_distinct_his_on_anchor'] > 0).sum())}")
    print(f"geometry + sign + liabilities: "
          f"{int((df['geometry_pass'] & df['sign_pass'] & df['liabilities_pass']).sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
