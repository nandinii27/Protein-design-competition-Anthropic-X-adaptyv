"""
Build the pH-switch motif specification for the TNF-alpha cation triad.

For each conserved target cation at the exposed protomer boundary, computes the
outward solvent direction and places an ideal histidine imidazole position
opposite it: close enough that a neutral histidine hydrogen bonds, close enough
that a protonated histidine suffers like-charge repulsion.

Outputs the switch geometry, a hotspot string for binder design, and a PDB of
pseudo-atoms marking the intended histidine positions for visual inspection and
for use as motif-scaffolding targets.

Usage:
    python make_switch_motif.py --pdb 1TNF.cif --out-prefix switch_motif

Dependencies:
    pip install biopython
"""

import argparse
import json
import warnings

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser

warnings.filterwarnings("ignore")

# Conserved positive residues on the exposed inter-protomer boundary.
# Tip atom is the formal charge centre of the side chain.
DEFAULT_CATIONS = [
    ("A", 103, "ARG", "CZ"),
    ("B", 112, "LYS", "NZ"),
    ("B", 103, "ARG", "CZ"),
]

# Separation between the target cation charge centre and the histidine ring
# centroid. At this range a neutral imidazole can hydrogen bond through its
# ring nitrogens, while a protonated imidazole carries a full positive charge
# against the target cation.
HIS_STANDOFF = 5.0

# Minimum allowed distance from any target heavy atom, to keep the placed
# histidine in solvent rather than inside the protein.
MIN_CLEARANCE = 3.6


def load_structure(path):
    if path.lower().endswith((".cif", ".mmcif")):
        return MMCIFParser(QUIET=True).get_structure("target", path)
    return PDBParser(QUIET=True).get_structure("target", path)


def outward_direction(model, point, radius=11.0):
    """
    Unit vector pointing from the local protein mass out into solvent.
    Computed as the direction from the centroid of nearby heavy atoms to the
    point of interest.
    """
    local = []
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            for atom in res:
                if atom.element == "H":
                    continue
                d = np.linalg.norm(atom.coord - point)
                if 0.1 < d <= radius:
                    local.append(atom.coord)
    if not local:
        raise RuntimeError("No neighbouring atoms found; check the input coordinates.")
    centroid = np.mean(local, axis=0)
    vec = point - centroid
    norm = np.linalg.norm(vec)
    if norm < 1e-6:
        raise RuntimeError("Degenerate outward direction.")
    return vec / norm


def clearance(model, point):
    best = 1e9
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            for atom in res:
                if atom.element == "H":
                    continue
                best = min(best, float(np.linalg.norm(atom.coord - point)))
    return best


def place_histidine_site(model, tip, standoff):
    """
    Place the histidine ring centroid along the outward normal from the cation
    tip, pushing outward until it clears the protein surface.
    """
    direction = outward_direction(model, tip)
    for extra in np.arange(0.0, 4.01, 0.25):
        candidate = tip + direction * (standoff + extra)
        if clearance(model, candidate) >= MIN_CLEARANCE:
            return candidate, direction, standoff + extra
    candidate = tip + direction * (standoff + 4.0)
    return candidate, direction, standoff + 4.0


def write_pseudo_pdb(path, sites):
    with open(path, "w") as fh:
        serial = 1
        for idx, site in enumerate(sites, start=1):
            x, y, z = site["target_tip"]
            fh.write(
                f"HETATM{serial:>5}  ZN  CAT X{idx:>4}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          ZN\n"
            )
            serial += 1
            x, y, z = site["his_centroid"]
            fh.write(
                f"HETATM{serial:>5}  O   HSP X{idx + 50:>4}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00           O\n"
            )
            serial += 1
        fh.write("END\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--standoff", type=float, default=HIS_STANDOFF)
    ap.add_argument("--out-prefix", default="switch_motif")
    args = ap.parse_args()

    structure = load_structure(args.pdb)
    model = structure[0]

    sites = []
    for chain_id, resnum, expect_name, tip_atom in DEFAULT_CATIONS:
        res = model[chain_id][(" ", resnum, " ")]
        if res.get_resname() != expect_name:
            raise SystemExit(
                f"Expected {expect_name} at {chain_id}{resnum}, "
                f"found {res.get_resname()}. Check the numbering of the input file."
            )
        if tip_atom not in res:
            raise SystemExit(f"Atom {tip_atom} missing from {chain_id}{resnum}.")
        tip = res[tip_atom].coord.astype(float)
        centroid, direction, used = place_histidine_site(model, tip, args.standoff)
        sites.append({
            "target_chain": chain_id,
            "target_resnum": resnum,
            "target_resname": expect_name,
            "target_tip": [round(float(v), 3) for v in tip],
            "his_centroid": [round(float(v), 3) for v in centroid],
            "outward_unit_vector": [round(float(v), 4) for v in direction],
            "standoff_used": round(float(used), 2),
            "clearance": round(float(clearance(model, centroid)), 2),
        })

    print("pH-switch motif sites")
    print()
    print(f"{'target':<10}{'standoff':>10}{'clearance':>11}   histidine ring centroid")
    for s in sites:
        label = f"{s['target_chain']}{s['target_resnum']} {s['target_resname']}"
        xyz = ", ".join(f"{v:8.3f}" for v in s["his_centroid"])
        print(f"{label:<10}{s['standoff_used']:>10.2f}{s['clearance']:>11.2f}   [{xyz}]")

    print()
    print("Separation between histidine sites (A):")
    n = len(sites)
    for i in range(n):
        for j in range(i + 1, n):
            a = np.array(sites[i]["his_centroid"])
            b = np.array(sites[j]["his_centroid"])
            li = f"{sites[i]['target_chain']}{sites[i]['target_resnum']}"
            lj = f"{sites[j]['target_chain']}{sites[j]['target_resnum']}"
            print(f"  {li} - {lj}: {np.linalg.norm(a - b):.1f}")

    span = max(
        float(np.linalg.norm(np.array(sites[i]["his_centroid"]) - np.array(sites[j]["his_centroid"])))
        for i in range(n) for j in range(n)
    )
    print(f"\nMaximum span across switch sites: {span:.1f} A")

    hotspots = ",".join(f"{s['target_chain']}{s['target_resnum']}" for s in sites)
    print(f"\nHotspot string (switch sites): {hotspots}")

    spec_path = f"{args.out_prefix}.json"
    with open(spec_path, "w") as fh:
        json.dump({"standoff_target": args.standoff, "sites": sites,
                   "hotspots": hotspots, "max_span": round(span, 2)}, fh, indent=2)

    pdb_path = f"{args.out_prefix}_sites.pdb"
    write_pseudo_pdb(pdb_path, sites)
    print(f"\nWrote {spec_path} and {pdb_path}")
    print("Load the pseudo-atom file alongside the target to inspect the geometry.")


if __name__ == "__main__":
    main()
