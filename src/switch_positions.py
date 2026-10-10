"""
Identify binder positions that face the target cation triad, and emit a
ProteinMPNN bias file that drives histidine into those positions.

Takes a designed binder-target complex and the switch geometry produced by
make_switch_motif.py. For each switch site, finds binder residues whose side
chain projects toward the target cation within the distance window where a
neutral imidazole can hydrogen bond and a protonated one would clash.

Reports which switch sites a given backbone actually covers, so backbones that
cannot support the switch are discarded before sequence design rather than
after.

Usage:
    python switch_positions.py --complex design.pdb --motif switch_motif.json \\
        --binder-chain B --out-prefix design_switch

Dependencies:
    pip install biopython
"""

import argparse
import glob
import json
import os
import warnings

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser

warnings.filterwarnings("ignore")

# ProteinMPNN alphabet order.
MPNN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
HIS_INDEX = MPNN_ALPHABET.index("H")

# A histidine engages the target cation when its side-chain centre sits in this
# range of the cation charge centre. Closer than the minimum and the neutral
# form clashes sterically; further and the protonated form carries no
# meaningful repulsion.
MIN_ENGAGE = 3.5
MAX_ENGAGE = 8.5


def load_structure(path):
    if path.lower().endswith((".cif", ".mmcif")):
        return MMCIFParser(QUIET=True).get_structure("s", path)
    return PDBParser(QUIET=True).get_structure("s", path)


def ordered_ca(model, chain_ids):
    """CA coordinates of the named chains, in file order."""
    coords = []
    for cid in chain_ids:
        if cid not in [c.id for c in model]:
            continue
        for res in model[cid]:
            if res.id[0] == " " and "CA" in res:
                coords.append(res["CA"].coord.astype(float))
    return np.array(coords)


def align_transform(ref_pts, out_pts):
    """
    Rigid transform mapping reference coordinates onto the design's frame.
    Design pipelines commonly recentre the system, so motif coordinates
    recorded against the original structure must be moved before use.
    """
    ref_c = ref_pts.mean(axis=0)
    out_c = out_pts.mean(axis=0)
    H = (ref_pts - ref_c).T @ (out_pts - out_c)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    rmsd = float(np.sqrt((((ref_pts - ref_c) @ R.T - (out_pts - out_c)) ** 2).sum(axis=1).mean()))
    return R, ref_c, out_c, rmsd


def apply_transform(point, transform):
    R, ref_c, out_c = transform
    return R @ (np.asarray(point, dtype=float) - ref_c) + out_c


def virtual_cb(res):
    """
    Build a CB position from backbone N, CA and C using ideal tetrahedral
    geometry. Generated backbones are often poly-glycine, which has no CB, so
    the real atom cannot be relied on.
    """
    if not all(a in res for a in ("N", "CA", "C")):
        return None
    n = res["N"].coord.astype(float)
    ca = res["CA"].coord.astype(float)
    c = res["C"].coord.astype(float)
    b = ca - n
    cc = c - ca
    a = np.cross(b, cc)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * cc + ca


def side_chain_point(res):
    """
    Approximate where a side chain placed at this position would sit: the CB
    position, projected one bond length further along the CA to CB direction so
    the estimate reflects a side-chain centre rather than its first atom.
    """
    if "CA" not in res:
        return None
    ca = res["CA"].coord.astype(float)
    cb = res["CB"].coord.astype(float) if "CB" in res else virtual_cb(res)
    if cb is None:
        return None
    direction = cb - ca
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        return cb
    return cb + (direction / norm) * 2.0


def analyse(complex_path, sites, binder_chain, reference=None, target_chains=("A", "B")):
    structure = load_structure(complex_path)
    model = structure[0]
    if binder_chain not in [c.id for c in model]:
        raise SystemExit(
            f"Chain {binder_chain} not found in {complex_path}. "
            f"Chains present: {','.join(c.id for c in model)}")

    transform = None
    align_rmsd = None
    if reference is not None:
        ref_pts = ordered_ca(reference[0], target_chains)
        out_pts = ordered_ca(model, target_chains)
        if len(ref_pts) and len(ref_pts) == len(out_pts):
            R, ref_c, out_c, align_rmsd = align_transform(ref_pts, out_pts)
            transform = (R, ref_c, out_c)

    assignments = []
    for res in model[binder_chain]:
        if res.id[0] != " ":
            continue
        point = side_chain_point(res)
        if point is None:
            continue
        for idx, site in enumerate(sites):
            tip = np.array(site["target_tip"], dtype=float)
            if transform is not None:
                tip = apply_transform(tip, transform)
            dist = float(np.linalg.norm(point - tip))
            if MIN_ENGAGE <= dist <= MAX_ENGAGE:
                assignments.append({
                    "binder_resnum": res.id[1],
                    "binder_resname": res.get_resname(),
                    "site_index": idx,
                    "site_label": f"{site['target_chain']}{site['target_resnum']}",
                    "distance": round(dist, 2),
                })

    # Keep the single best binder position per switch site.
    best = {}
    for a in assignments:
        cur = best.get(a["site_index"])
        if cur is None or a["distance"] < cur["distance"]:
            best[a["site_index"]] = a
    return (structure, sorted(best.values(), key=lambda a: a["site_index"]),
            assignments, align_rmsd)


def write_bias_jsonl(path, entries):
    with open(path, "w") as fh:
        for name, chains in entries.items():
            fh.write(json.dumps({"name": name, **chains}) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complex", nargs="+", required=True,
                    help="designed binder-target complex PDB files (globs accepted)")
    ap.add_argument("--motif", required=True, help="switch_motif.json")
    ap.add_argument("--binder-chain", default="C")
    ap.add_argument("--reference", default=None,
                    help="original target structure the motif was built from; "
                         "used to undo any recentring applied by the design pipeline")
    ap.add_argument("--target-chains", default="A,B")
    ap.add_argument("--bias-strength", type=float, default=3.0,
                    help="logit bias added to histidine at covered positions")
    ap.add_argument("--min-sites", type=int, default=2,
                    help="discard backbones covering fewer switch sites than this")
    ap.add_argument("--out-prefix", default="switch_positions")
    args = ap.parse_args()

    with open(args.motif) as fh:
        motif = json.load(fh)
    sites = motif["sites"]

    paths = []
    for pattern in args.complex:
        paths.extend(sorted(glob.glob(pattern)))
    if not paths:
        raise SystemExit("No input structures matched.")

    reference = load_structure(args.reference) if args.reference else None
    target_chains = tuple(c for c in args.target_chains.split(",") if c)

    bias_entries = {}
    report = []
    kept = 0

    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        try:
            structure, best, _, align_rmsd = analyse(
                path, sites, args.binder_chain, reference, target_chains)
        except SystemExit:
            raise
        except Exception as exc:
            report.append({"design": name, "sites_covered": 0,
                           "positions": "", "note": f"failed: {exc}"})
            continue

        covered = len(best)
        if align_rmsd is not None and align_rmsd > 1.0:
            print(f"warning: {name} target does not superpose on the reference "
                  f"(rmsd {align_rmsd:.2f} A)")
        positions = [a["binder_resnum"] for a in best]
        labels = ";".join(f"{a['site_label']}->{a['binder_resnum']}@{a['distance']}A"
                          for a in best)

        if covered >= args.min_sites:
            kept += 1
            model = structure[0]
            chain_len = len([r for r in model[args.binder_chain] if r.id[0] == " "])
            resnum_to_index = {}
            for i, r in enumerate(r for r in model[args.binder_chain] if r.id[0] == " "):
                resnum_to_index[r.id[1]] = i

            matrix = [[0.0] * len(MPNN_ALPHABET) for _ in range(chain_len)]
            for rn in positions:
                i = resnum_to_index.get(rn)
                if i is not None:
                    matrix[i][HIS_INDEX] = args.bias_strength
            bias_entries[name] = {args.binder_chain: matrix}

        report.append({"design": name, "sites_covered": covered,
                       "positions": " ".join(str(p) for p in positions),
                       "note": labels})

    bias_path = f"{args.out_prefix}_bias.jsonl"
    write_bias_jsonl(bias_path, bias_entries)

    csv_path = f"{args.out_prefix}_report.csv"
    header = ["design", "sites_covered", "positions", "note"]
    with open(csv_path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in sorted(report, key=lambda r: -r["sites_covered"]):
            fh.write(",".join('"' + str(r[c]).replace('"', "'") + '"' for c in header) + "\n")

    print(f"Analysed {len(report)} structures against {len(sites)} switch sites.")
    print(f"Backbones covering at least {args.min_sites} sites: {kept}")
    print()
    print(f"{'design':<34}{'sites':>7}  positions")
    for r in sorted(report, key=lambda r: -r["sites_covered"])[:25]:
        print(f"{r['design']:<34}{r['sites_covered']:>7}  {r['positions']}")
    print()
    print(f"Wrote {bias_path} and {csv_path}")


if __name__ == "__main__":
    main()
