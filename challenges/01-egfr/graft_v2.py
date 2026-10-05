"""Burial-aware multi-histidine grafting with explicit rotamer placement.

The first grafting pass selected positions on reach and orientation alone and
allowed at most two substitutions. Both limits cap the achievable switch.

A pKa shift comes from the environment changing around the residue. A histidine
that stays solvent-exposed in both the free and bound states sees almost no
change, so its pKa barely moves. Large shifts need the side chain to become
buried on binding, where desolvation and a nearby carboxylate stabilise the
protonated form. This pass therefore requires a drop in side-chain accessibility
on complex formation.

Magnitude also forces multiple sites. For a single titratable site, even a
two-unit pKa shift yields roughly 1 kcal/mol of pH selectivity between 6.5 and
7.4, about five-fold. Reaching ten- to hundred-fold needs two or three
histidines acting together, so combinations up to three are enumerated.

Search is staged by cost. Candidate positions and their rotamers are scored
geometrically, which is instant; combinations are ranked on that; only the
shortlist is passed to the pKa predictor; only what survives that is written
out for re-prediction.

Reads:
    final_design_stats.csv
    outputs/complexes/<Design>*.pdb      target chain A, binder chain B
Writes:
    outputs/graft_v2_candidates.csv      positions, burial, rotamer quality
    outputs/graft_v2_combos.csv          every combination scored
    outputs/mutants_v2.csv               shortlist with sequences
    outputs/colabfold_input_v2.csv       id,sequence for re-prediction
"""

from __future__ import annotations

import math
import shutil
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.SASA import ShrakeRupley

from binderkit import phscore

CHALLENGE_DIR = Path(__file__).resolve().parent
OUTPUTS = CHALLENGE_DIR / "outputs"
COMPLEX_DIR = OUTPUTS / "complexes"
STATS_CSV = CHALLENGE_DIR / "final_design_stats.csv"
SCRATCH = OUTPUTS / "graft_v2_work"

TARGET_CHAIN = "A"
BINDER_CHAIN = "B"
ANCHORS = (21, 24, 65, 68)
PH_LOW, PH_HIGH = 6.5, 7.4

# Candidate position criteria.
MIN_CB_DIST = 4.0
MAX_CB_DIST = 11.0
MIN_COSINE = -0.1
MIN_BURIAL = 15.0        # square angstrom of side-chain SASA lost on binding
MIN_FREE_SASA = 20.0     # must be surface exposed in the free binder

MAX_CANDIDATES_PER_DESIGN = 8
MAX_HIS_PER_DESIGN = 3
SHORTLIST_PER_DESIGN = 6
MAX_TOTAL_MUTANTS = 48

# Imidazole nitrogen to carboxylate oxygen, the contact being designed for.
IDEAL_CONTACT = 3.0
CONTACT_TOLERANCE = 1.2
CLASH_DISTANCE = 2.9

SKIP_RESIDUES = {"GLY", "PRO", "CYS", "HIS"}
CARBOXYLATE_OXYGENS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
BACKBONE = {"N", "CA", "C", "O", "OXT"}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Histidine side chain internal coordinates: (atom, parents, bond, angle, dihedral)
# Dihedral of None means the value is taken from the sampled chi angles.
HIS_GEOMETRY = [
    ("CG", ("N", "CA", "CB"), 1.497, 113.8, "chi1"),
    ("ND1", ("CA", "CB", "CG"), 1.378, 122.7, "chi2"),
    ("CD2", ("CA", "CB", "CG"), 1.354, 131.0, "chi2+180"),
    ("CE1", ("CB", "CG", "ND1"), 1.321, 109.0, 180.0),
    ("NE2", ("CB", "CG", "CD2"), 1.374, 107.2, 180.0),
]

CHI1_ROTAMERS = (-177.0, -67.0, 62.0)
CHI2_ROTAMERS = (-165.0, -80.0, -75.0, 60.0, 80.0, 100.0)

RING_ATOMS = ("CG", "ND1", "CD2", "CE1", "NE2")
NITROGENS = ("ND1", "NE2")


class _ChainSelect(Select):
    def __init__(self, keep: set[str]):
        self.keep = keep

    def accept_chain(self, chain):
        return chain.get_id() in self.keep

    def accept_residue(self, residue):
        return residue.get_id()[0] == " "


def place_atom(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, bond: float, angle: float, dihedral: float
) -> np.ndarray:
    """Natural extension reference frame: place D given A-B-C and internal coords.

    Standard construction used to build side chains from backbone geometry.
    """
    angle = math.radians(angle)
    dihedral = math.radians(dihedral)

    bc = c - b
    bc /= np.linalg.norm(bc)
    ab = b - a
    normal = np.cross(ab, bc)
    norm = np.linalg.norm(normal)
    if norm < 1e-8:
        normal = np.array([0.0, 0.0, 1.0])
    else:
        normal /= norm
    cross = np.cross(normal, bc)

    local = np.array(
        [
            -bond * math.cos(angle),
            bond * math.sin(angle) * math.cos(dihedral),
            bond * math.sin(angle) * math.sin(dihedral),
        ]
    )
    basis = np.array([bc, cross, normal]).T
    return c + basis.dot(local)


def build_histidine(residue, chi1: float, chi2: float) -> dict[str, np.ndarray] | None:
    """Side chain coordinates for one histidine rotamer on this backbone."""
    try:
        coords = {
            "N": residue["N"].get_coord().astype(float),
            "CA": residue["CA"].get_coord().astype(float),
            "CB": residue["CB"].get_coord().astype(float),
        }
    except KeyError:
        return None

    for name, parents, bond, angle, dihedral in HIS_GEOMETRY:
        if dihedral == "chi1":
            value = chi1
        elif dihedral == "chi2":
            value = chi2
        elif dihedral == "chi2+180":
            value = chi2 + 180.0
        else:
            value = dihedral
        try:
            a, b, c = (coords[p] for p in parents)
        except KeyError:
            return None
        coords[name] = place_atom(a, b, c, bond, angle, value)

    return {k: v for k, v in coords.items() if k not in ("N", "CA")}


def side_chain_sasa(structure_model, chain_id: str) -> dict[int, float]:
    """Per-residue SASA for one chain of whatever model is passed in."""
    ShrakeRupley().compute(structure_model, level="R")
    return {
        residue.get_id()[1]: float(residue.sasa)
        for residue in structure_model[chain_id]
        if residue.get_id()[0] == " "
    }


def anchor_oxygen_coords(model) -> dict[int, np.ndarray]:
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


def environment_atoms(model, exclude_resnum: int) -> np.ndarray:
    """Heavy atom coordinates near which a placed ring must not clash."""
    coords = []
    for chain in model:
        for residue in chain:
            if residue.get_id()[0] != " ":
                continue
            if chain.get_id() == BINDER_CHAIN and residue.get_id()[1] == exclude_resnum:
                continue
            for atom in residue:
                if atom.element != "H":
                    coords.append(atom.get_coord())
    return np.array(coords, dtype=float)


def score_rotamer(
    ring: dict[str, np.ndarray],
    anchors: dict[int, np.ndarray],
    environment: np.ndarray,
) -> dict | None:
    """Clash and contact quality for one placed rotamer.

    Returns None when the rotamer clashes badly enough to be unusable.
    """
    ring_coords = np.array([ring[a] for a in RING_ATOMS if a in ring])
    if len(ring_coords) < len(RING_ATOMS):
        return None

    distances = np.linalg.norm(
        ring_coords[:, None, :] - environment[None, :, :], axis=-1
    )
    clashes = int((distances < CLASH_DISTANCE).sum())
    if clashes > 2:
        return None

    best = None
    for anchor, oxygens in anchors.items():
        for nitrogen in NITROGENS:
            if nitrogen not in ring:
                continue
            separations = np.linalg.norm(oxygens - ring[nitrogen], axis=1)
            distance = float(separations.min())
            deviation = abs(distance - IDEAL_CONTACT)
            if deviation > CONTACT_TOLERANCE:
                continue
            quality = 1.0 - deviation / CONTACT_TOLERANCE
            if best is None or quality > best["contact_quality"]:
                best = {
                    "anchor": anchor,
                    "nitrogen": nitrogen,
                    "contact_distance": round(distance, 2),
                    "contact_quality": round(quality, 3),
                    "clashes": clashes,
                }
    return best


def candidate_positions(model, free_sasa, bound_sasa, anchors) -> list[dict]:
    """Positions that are exposed when free, buried on binding, and can reach."""
    found: list[dict] = []
    for residue in model[BINDER_CHAIN]:
        if residue.get_id()[0] != " ":
            continue
        resname = residue.get_resname()
        resnum = residue.get_id()[1]
        if resname in SKIP_RESIDUES or "CA" not in residue or "CB" not in residue:
            continue

        free = free_sasa.get(resnum, 0.0)
        bound = bound_sasa.get(resnum, 0.0)
        burial = free - bound
        if free < MIN_FREE_SASA or burial < MIN_BURIAL:
            continue

        ca = residue["CA"].get_coord().astype(float)
        cb = residue["CB"].get_coord().astype(float)
        direction = cb - ca
        direction /= max(np.linalg.norm(direction), 1e-8)

        reachable = False
        for oxygens in anchors.values():
            offsets = oxygens - cb
            distances = np.linalg.norm(offsets, axis=1)
            index = int(np.argmin(distances))
            distance = float(distances[index])
            if not (MIN_CB_DIST <= distance <= MAX_CB_DIST):
                continue
            if float(np.dot(direction, offsets[index] / distance)) >= MIN_COSINE:
                reachable = True
                break
        if not reachable:
            continue

        environment = environment_atoms(model, resnum)
        best = None
        for chi1 in CHI1_ROTAMERS:
            for chi2 in CHI2_ROTAMERS:
                ring = build_histidine(residue, chi1, chi2)
                if ring is None:
                    continue
                scored = score_rotamer(ring, anchors, environment)
                if scored is None:
                    continue
                if best is None or scored["contact_quality"] > best["contact_quality"]:
                    best = {**scored, "chi1": chi1, "chi2": chi2, "ring": ring}
        if best is None:
            continue

        found.append(
            {
                "binder_resnum": resnum,
                "binder_resname": resname,
                "free_sasa": round(free, 1),
                "bound_sasa": round(bound, 1),
                "burial": round(burial, 1),
                "anchor": best["anchor"],
                "nitrogen": best["nitrogen"],
                "contact_distance": best["contact_distance"],
                "contact_quality": best["contact_quality"],
                "clashes": best["clashes"],
                "chi1": best["chi1"],
                "chi2": best["chi2"],
                "ring": best["ring"],
                # Burial is capped so that one deeply buried position cannot
                # outweigh good contact geometry at several others.
                "score": round(
                    best["contact_quality"] * (1.0 + min(burial, 60.0) / 60.0), 4
                ),
            }
        )

    found.sort(key=lambda c: c["score"], reverse=True)
    return found[:MAX_CANDIDATES_PER_DESIGN]


def combo_score(members: list[dict]) -> float:
    """Geometric score for a set of simultaneous substitutions.

    Distinct anchors are rewarded: independent salt bridges contribute more
    protonation-linked free energy than several contacts to one carboxylate.
    """
    base = sum(m["score"] for m in members)
    distinct = len({m["anchor"] for m in members})
    return round(base * (1.0 + 0.25 * (distinct - 1)), 4)


def mutate_structure(model, members: list[dict]):
    """Replace side chains with the chosen histidine rotamers, in place."""
    for member in members:
        residue = model[BINDER_CHAIN][(" ", member["binder_resnum"], " ")]
        for atom in [a.get_name() for a in residue if a.get_name() not in BACKBONE]:
            residue.detach_child(atom)
        for index, (name, coord) in enumerate(member["ring"].items()):
            element = "N" if name.startswith("N") else "C"
            residue.add(
                Atom(
                    name,
                    np.array(coord, dtype=float),
                    20.0,
                    1.0,
                    " ",
                    name,
                    1000 + index,
                    element,
                )
            )
        residue.resname = "HIS"
    return model


def write_structure(structure, path: Path, chains: set[str] | None = None) -> Path:
    io = PDBIO()
    io.set_structure(structure)
    path.parent.mkdir(parents=True, exist_ok=True)
    if chains is None:
        io.save(str(path))
    else:
        io.save(str(path), _ChainSelect(chains))
    return path


def chain_sequence(model, chain_id: str) -> str:
    return "".join(
        THREE_TO_ONE.get(r.get_resname(), "X")
        for r in model[chain_id]
        if r.get_id()[0] == " "
    )


def complex_path(design: str) -> Path | None:
    direct = COMPLEX_DIR / f"{design}.pdb"
    if direct.exists():
        return direct
    matches = sorted(COMPLEX_DIR.glob(f"{design}*.pdb"))
    return matches[0] if matches else None


def score_with_propka(structure, workdir: Path) -> dict | None:
    """Linkage score for a mutated complex, via pKa prediction free and bound."""
    workdir.mkdir(parents=True, exist_ok=True)
    bound = write_structure(structure, workdir / "b.pdb", {TARGET_CHAIN, BINDER_CHAIN})
    free = write_structure(structure, workdir / "f.pdb", {BINDER_CHAIN})

    try:
        pka_bound = phscore.run_propka(bound, workdir / "wb")
        pka_free = phscore.run_propka(free, workdir / "wf")
    except Exception:
        return None

    sites = phscore.interface_sites(bound, BINDER_CHAIN, TARGET_CHAIN)
    low = high = 0.0
    shifts = []
    for site in sites:
        f = pka_free.get(site)
        b = pka_bound.get(site)
        if f is None or b is None:
            continue
        low += phscore.linkage_ddg(f, b, PH_LOW)
        high += phscore.linkage_ddg(f, b, PH_HIGH)
        if site.resname == "HIS":
            shifts.append(b - f)

    contacts = phscore.his_carboxylate_contacts(bound, BINDER_CHAIN, TARGET_CHAIN)
    return {
        "ph_selectivity_kcal": round(low - high, 3),
        "ddg_ph_low": round(low, 3),
        "ddg_ph_high": round(high, 3),
        "max_his_pka_shift": round(max(shifts), 2) if shifts else None,
        "n_his_carboxylate_contacts": len(contacts),
        "n_distinct_his_in_contact": len({c["his_resnum"] for c in contacts}),
    }


def main() -> int:
    if not STATS_CSV.exists() or not COMPLEX_DIR.exists():
        print("missing stats file or complex directory", file=sys.stderr)
        return 1

    designs = pd.read_csv(STATS_CSV)
    designs = designs[designs["Design"].notna()].copy()
    designs["Design"] = designs["Design"].astype(str).str.strip()
    designs["Sequence"] = designs["Sequence"].astype(str).str.strip().str.upper()

    if SCRATCH.exists():
        shutil.rmtree(SCRATCH, ignore_errors=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True)
    candidate_rows, combo_rows, mutant_rows = [], [], []
    target_sequence = None
    work_index = 0

    for _, row in designs.iterrows():
        design = row["Design"]
        path = complex_path(design)
        if path is None:
            continue

        model = parser.get_structure("c", str(path))[0]
        if target_sequence is None:
            target_sequence = chain_sequence(model, TARGET_CHAIN)
        binder_sequence = chain_sequence(model, BINDER_CHAIN)
        if binder_sequence != row["Sequence"]:
            print(f"{design}: sequence mismatch, skipping", file=sys.stderr)
            continue

        bound_sasa = side_chain_sasa(
            parser.get_structure("c", str(path))[0], BINDER_CHAIN
        )
        binder_only = write_structure(
            parser.get_structure("c", str(path)), SCRATCH / design / "binder.pdb",
            {BINDER_CHAIN},
        )
        free_sasa = side_chain_sasa(
            parser.get_structure("b", str(binder_only))[0], BINDER_CHAIN
        )

        anchors = anchor_oxygen_coords(model)
        candidates = candidate_positions(model, free_sasa, bound_sasa, anchors)
        for candidate in candidates:
            candidate_rows.append(
                {"Design": design, **{k: v for k, v in candidate.items() if k != "ring"}}
            )
        if not candidates:
            print(f"{design}: no buried reachable position")
            continue

        combos = []
        for size in range(1, MAX_HIS_PER_DESIGN + 1):
            for members in combinations(candidates, size):
                combos.append((combo_score(list(members)), list(members)))
        combos.sort(key=lambda pair: pair[0], reverse=True)

        for score, members in combos:
            combo_rows.append(
                {
                    "Design": design,
                    "positions": "+".join(str(m["binder_resnum"]) for m in members),
                    "anchors": "+".join(str(m["anchor"]) for m in members),
                    "n_his": len(members),
                    "total_burial": round(sum(m["burial"] for m in members), 1),
                    "geometric_score": score,
                }
            )

        shortlist = combos[:SHORTLIST_PER_DESIGN]
        print(f"{design}: {len(candidates)} positions, {len(combos)} combinations, "
              f"scoring top {len(shortlist)}")

        for score, members in shortlist:
            work_index += 1
            structure = parser.get_structure("m", str(path))
            mutate_structure(structure[0], members)
            result = score_with_propka(structure, SCRATCH / f"w{work_index:04d}")
            if result is None:
                continue

            positions = sorted(m["binder_resnum"] for m in members)
            sequence = list(binder_sequence)
            for position in positions:
                sequence[position - 1] = "H"
            sequence = "".join(sequence)

            mutant_rows.append(
                {
                    "Design": design,
                    "mutant_id": f"{design}_v2_" + "_".join(f"H{p}" for p in positions),
                    "mutations": "+".join(
                        f"{m['binder_resname']}{m['binder_resnum']}H" for m in members
                    ),
                    "anchors": "+".join(str(m["anchor"]) for m in members),
                    "n_mutations": len(members),
                    "total_burial": round(sum(m["burial"] for m in members), 1),
                    "mean_contact_quality": round(
                        float(np.mean([m["contact_quality"] for m in members])), 3
                    ),
                    "geometric_score": score,
                    "parent_sequence": binder_sequence,
                    "sequence": sequence,
                    **result,
                }
            )
            print(
                f"  {'+'.join(str(p) for p in positions)}: "
                f"{result['ph_selectivity_kcal']:+.3f} kcal/mol, "
                f"burial {sum(m['burial'] for m in members):.0f} A^2"
            )

    if not mutant_rows:
        print("no mutants produced", file=sys.stderr)
        return 1

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(candidate_rows).to_csv(OUTPUTS / "graft_v2_candidates.csv", index=False)
    pd.DataFrame(combo_rows).to_csv(OUTPUTS / "graft_v2_combos.csv", index=False)

    mutants = pd.DataFrame(mutant_rows).drop_duplicates(subset=["sequence"])
    mutants = mutants.sort_values("ph_selectivity_kcal").head(MAX_TOTAL_MUTANTS)
    mutants.to_csv(OUTPUTS / "mutants_v2.csv", index=False)

    pd.DataFrame(
        {
            "id": mutants["mutant_id"],
            "sequence": target_sequence + ":" + mutants["sequence"],
        }
    ).to_csv(OUTPUTS / "colabfold_input_v2.csv", index=False)

    negative = mutants[mutants["ph_selectivity_kcal"] < 0]
    print(f"\n{len(mutants)} mutants kept from {mutants['Design'].nunique()} designs")
    print(f"negative selectivity: {len(negative)}")
    if not negative.empty:
        print(f"best: {negative['ph_selectivity_kcal'].min():+.3f} kcal/mol")
        print(negative.groupby("n_mutations")["ph_selectivity_kcal"].min().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
