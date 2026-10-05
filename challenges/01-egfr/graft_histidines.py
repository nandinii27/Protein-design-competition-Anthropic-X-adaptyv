"""Identify binder positions able to reach epitope carboxylates and write His mutants.

The generator optimises binding with no pH objective, so any histidine at the
interface is incidental. This step places histidine deliberately: it finds
binder positions whose side chain could form an imidazole-carboxylate contact
with a chosen anchor, mutates them, and emits sequences for re-prediction.

Geometric criteria, set by histidine side chain dimensions:
    CB to anchor carboxylate oxygen within [MIN_CB_DIST, MAX_CB_DIST]
    CA->CB direction pointing towards that oxygen

Closer than the lower bound and the imidazole clashes; further than the upper
bound and it cannot reach. The direction test rejects positions whose side
chain points into the binder core.

Reads:
    final_design_stats.csv
    outputs/complexes/<Design>*.pdb      target chain A, binder chain B
    outputs/merged_metrics.csv           optional, to restrict to clean designs

Writes:
    outputs/graft_candidates.csv         every position considered, with geometry
    outputs/mutants.csv                  mutant sequences and their mutations
    outputs/colabfold_input.csv          id,sequence for batch re-prediction
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

CHALLENGE_DIR = Path(__file__).resolve().parent
OUTPUTS = CHALLENGE_DIR / "outputs"
COMPLEX_DIR = OUTPUTS / "complexes"
STATS_CSV = CHALLENGE_DIR / "final_design_stats.csv"
MERGED_CSV = OUTPUTS / "merged_metrics.csv"

TARGET_CHAIN = "A"
BINDER_CHAIN = "B"

# Anchor carboxylates in the renumbered complex. These correspond to the
# conserved exposed acidic residues selected during epitope definition.
ANCHORS = (21, 24, 65, 68)

MIN_CB_DIST = 4.5
MAX_CB_DIST = 10.0
MIN_COSINE = 0.0

# Residues not worth mutating: no CB, backbone constrained, or already His.
SKIP_RESIDUES = {"GLY", "PRO", "CYS", "HIS"}

MAX_SINGLES_PER_DESIGN = 4
MAX_PAIRS_PER_DESIGN = 3

# Only graft designs that passed the liability checks, when that table exists.
RESTRICT_TO_CLEAN = True

CARBOXYLATE_OXYGENS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def complex_path(design: str) -> Path | None:
    direct = COMPLEX_DIR / f"{design}.pdb"
    if direct.exists():
        return direct
    matches = sorted(COMPLEX_DIR.glob(f"{design}*.pdb"))
    return matches[0] if matches else None


def chain_sequence(model, chain_id: str) -> str:
    return "".join(
        THREE_TO_ONE.get(r.get_resname(), "X")
        for r in model[chain_id]
        if r.get_id()[0] == " "
    )


def anchor_oxygens(model) -> dict[int, np.ndarray]:
    """Carboxylate oxygen coordinates for each anchor residue."""
    out: dict[int, np.ndarray] = {}
    for residue in model[TARGET_CHAIN]:
        resnum = residue.get_id()[1]
        if resnum not in ANCHORS:
            continue
        names = CARBOXYLATE_OXYGENS.get(residue.get_resname())
        if names is None:
            continue
        coords = [residue[n].get_coord() for n in names if n in residue]
        if coords:
            out[resnum] = np.array(coords, dtype=float)
    return out


def find_candidates(model) -> list[dict]:
    """Binder positions geometrically able to reach an anchor carboxylate."""
    anchors = anchor_oxygens(model)
    if not anchors:
        return []

    candidates: list[dict] = []
    for residue in model[BINDER_CHAIN]:
        if residue.get_id()[0] != " ":
            continue
        resname = residue.get_resname()
        if resname in SKIP_RESIDUES or "CA" not in residue or "CB" not in residue:
            continue

        ca = residue["CA"].get_coord().astype(float)
        cb = residue["CB"].get_coord().astype(float)
        direction = cb - ca
        norm = np.linalg.norm(direction)
        if norm < 1e-6:
            continue
        direction = direction / norm

        for anchor, oxygens in anchors.items():
            offsets = oxygens - cb
            distances = np.linalg.norm(offsets, axis=1)
            index = int(np.argmin(distances))
            distance = float(distances[index])
            if not (MIN_CB_DIST <= distance <= MAX_CB_DIST):
                continue
            cosine = float(np.dot(direction, offsets[index] / distance))
            if cosine < MIN_COSINE:
                continue
            candidates.append(
                {
                    "binder_resnum": residue.get_id()[1],
                    "binder_resname": resname,
                    "anchor": anchor,
                    "cb_distance": round(distance, 2),
                    "cosine": round(cosine, 3),
                    # Prefer short, well oriented positions.
                    "score": round(cosine / distance, 4),
                }
            )
    candidates.sort(key=lambda c: c["score"], reverse=True)
    return candidates


def mutate(sequence: str, positions: list[int]) -> str:
    """Substitute histidine at 1-indexed binder positions."""
    residues = list(sequence)
    for position in positions:
        residues[position - 1] = "H"
    return "".join(residues)


def best_per_position(candidates: list[dict]) -> list[dict]:
    """Keep the best anchor per binder position."""
    seen: dict[int, dict] = {}
    for candidate in candidates:
        current = seen.get(candidate["binder_resnum"])
        if current is None or candidate["score"] > current["score"]:
            seen[candidate["binder_resnum"]] = candidate
    return sorted(seen.values(), key=lambda c: c["score"], reverse=True)


def main() -> int:
    if not STATS_CSV.exists():
        print(f"missing {STATS_CSV.name}", file=sys.stderr)
        return 1

    designs = pd.read_csv(STATS_CSV)
    designs = designs[designs["Design"].notna()].copy()
    designs["Design"] = designs["Design"].astype(str).str.strip()
    designs["Sequence"] = designs["Sequence"].astype(str).str.strip().str.upper()

    if RESTRICT_TO_CLEAN and MERGED_CSV.exists():
        merged = pd.read_csv(MERGED_CSV)
        if "passes" in merged.columns:
            clean = set(merged.loc[merged["passes"].astype(bool), "Design"])
            kept = designs[designs["Design"].isin(clean)]
            print(f"{len(kept)} of {len(designs)} designs passed liability checks")
            if len(kept) >= 4:
                designs = kept
            else:
                print("too few clean designs, grafting all of them instead")

    parser = PDBParser(QUIET=True)
    candidate_rows: list[dict] = []
    mutant_rows: list[dict] = []
    target_sequence: str | None = None

    for _, row in designs.iterrows():
        design = row["Design"]
        path = complex_path(design)
        if path is None:
            print(f"{design}: no complex found", file=sys.stderr)
            continue

        model = parser.get_structure("c", str(path))[0]
        if target_sequence is None:
            target_sequence = chain_sequence(model, TARGET_CHAIN)

        binder_sequence = chain_sequence(model, BINDER_CHAIN)
        if binder_sequence != row["Sequence"]:
            print(
                f"{design}: structure sequence differs from stats file, skipping",
                file=sys.stderr,
            )
            continue

        candidates = find_candidates(model)
        for candidate in candidates:
            candidate_rows.append({"Design": design, **candidate})

        unique = best_per_position(candidates)
        if not unique:
            print(f"{design}: no graftable position")
            continue

        singles = unique[:MAX_SINGLES_PER_DESIGN]
        for candidate in singles:
            position = candidate["binder_resnum"]
            mutant_rows.append(
                {
                    "Design": design,
                    "mutant_id": f"{design}_H{position}",
                    "mutations": f"{candidate['binder_resname']}{position}H",
                    "anchors": str(candidate["anchor"]),
                    "n_mutations": 1,
                    "parent_sequence": binder_sequence,
                    "sequence": mutate(binder_sequence, [position]),
                }
            )

        # Pairs engaging two different anchors: one salt bridge carries too
        # little protonation-linked free energy to switch binding on its own.
            pairs = [
            (a, b)
            for a, b in combinations(unique, 2)
            if a["binder_resnum"] != b["binder_resnum"]
        ][:MAX_PAIRS_PER_DESIGN]
        for first, second in pairs:
            positions = sorted((first["binder_resnum"], second["binder_resnum"]))
            mutant_rows.append(
                {
                    "Design": design,
                    "mutant_id": f"{design}_H{positions[0]}_H{positions[1]}",
                    "mutations": (
                        f"{first['binder_resname']}{first['binder_resnum']}H+"
                        f"{second['binder_resname']}{second['binder_resnum']}H"
                    ),
                    "anchors": f"{first['anchor']}+{second['anchor']}",
                    "n_mutations": 2,
                    "parent_sequence": binder_sequence,
                    "sequence": mutate(binder_sequence, positions),
                }
            )

        print(f"{design}: {len(unique)} positions, {len(singles)} singles, {len(pairs)} pairs")

    if not mutant_rows:
        print("no mutants generated", file=sys.stderr)
        return 1

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(candidate_rows).to_csv(OUTPUTS / "graft_candidates.csv", index=False)

    mutants = pd.DataFrame(mutant_rows).drop_duplicates(subset=["sequence"])
    mutants.to_csv(OUTPUTS / "mutants.csv", index=False)

    colabfold = pd.DataFrame(
        {
            "id": mutants["mutant_id"],
            "sequence": target_sequence + ":" + mutants["sequence"],
        }
    )
    colabfold.to_csv(OUTPUTS / "colabfold_input.csv", index=False)

    print(f"\n{len(mutants)} unique mutants from {mutants['Design'].nunique()} designs")
    print(f"singles {int((mutants['n_mutations'] == 1).sum())}, "
          f"pairs {int((mutants['n_mutations'] == 2).sum())}")
    print(f"target chain length {len(target_sequence)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
