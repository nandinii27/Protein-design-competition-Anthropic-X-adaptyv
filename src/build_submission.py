"""
Assemble the competition submission from ProteinMPNN output.

Extracts binder sequences, applies histidine at the switch positions identified
against the target cation triad, pairs each switch variant with its unmodified
counterpart as a control, and ranks by the structural metrics available.

Usage:
    python build_submission.py --mpnn mpnn_results.csv \\
        --switch-report beta_switch_report.csv \\
        --out submission.csv

Dependencies: none beyond the standard library.
"""

import argparse
import csv
import json
import os

# Order in which surviving backbones were renumbered for sequence design:
# position in this list is the design index MPNN saw, the value is the
# original RFdiffusion design number the switch report is keyed to.
RENUMBER = [1, 3, 7, 10, 12, 14, 15, 16, 17, 19, 25, 28]

PREFIX = "tnfa_switch_1r9qo_"


def load_switch_positions(path):
    """Map original design number -> list of binder positions facing a cation."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            name = row["design"]
            if not name.startswith(PREFIX):
                continue
            original = int(name[len(PREFIX):])
            pos = [int(p) for p in row["positions"].split() if p]
            if pos:
                out[original] = sorted(pos)
    return out


def binder_sequence(seq_field):
    """MPNN writes target/binder; the designed chain is the final segment."""
    return seq_field.split("/")[-1].strip()


def apply_switch(seq, positions):
    """
    Place histidine at the switch positions. Positions are 1-indexed residue
    numbers in the binder chain.
    """
    chars = list(seq)
    applied = []
    for p in positions:
        i = p - 1
        if 0 <= i < len(chars):
            applied.append((p, chars[i]))
            chars[i] = "H"
    return "".join(chars), applied


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpnn", required=True)
    ap.add_argument("--switch-report", required=True)
    ap.add_argument("--out", default="submission.csv")
    ap.add_argument("--limit", type=int, default=40)
    ap.add_argument("--controls", type=int, default=8,
                    help="unmodified counterparts included as paired controls")
    args = ap.parse_args()

    switch_pos = load_switch_positions(args.switch_report)

    records = []
    with open(args.mpnn) as fh:
        for row in csv.DictReader(fh):
            design = int(row["design"])
            if design >= len(RENUMBER):
                continue
            original = RENUMBER[design]
            seq = binder_sequence(row["seq"])
            positions = switch_pos.get(original, [])
            switched, applied = apply_switch(seq, positions)
            records.append({
                "design": design,
                "original": original,
                "n": int(row["n"]),
                "seq": seq,
                "switched": switched,
                "positions": positions,
                "applied": applied,
                "n_switch": len(applied),
                "rmsd": float(row["rmsd"]),
                "plddt": float(row["plddt"]),
                "i_ptm": float(row["i_ptm"]),
                "mpnn": float(row["mpnn"]),
            })

    # No design validated, so ranking uses the weak signals that remain:
    # self-consistency first, then fold confidence, then interface score.
    records.sort(key=lambda r: (r["rmsd"], -r["plddt"], -r["i_ptm"]))

    n_switch = args.limit - args.controls
    chosen = []
    seen = set()

    for r in records:
        if len(chosen) >= n_switch:
            break
        if r["switched"] in seen or r["n_switch"] == 0:
            continue
        seen.add(r["switched"])
        chosen.append({
            "id": f"tnfa_sw_d{r['original']}_n{r['n']}",
            "sequence": r["switched"],
            "arm": "switch",
            "switch_positions": " ".join(str(p) for p in r["positions"]),
            "replaced": " ".join(f"{aa}{p}H" for p, aa in r["applied"]),
            "rmsd": r["rmsd"],
            "plddt": r["plddt"],
            "i_ptm": r["i_ptm"],
        })

    # Paired controls: the same sequences without the histidines, so the
    # effect of the substitution is separable if anything binds.
    control_sources = [c for c in chosen[:args.controls]]
    by_id = {(r["original"], r["n"]): r for r in records}
    for c in control_sources:
        key = c["id"].replace("tnfa_sw_d", "").split("_n")
        orig, n = int(key[0]), int(key[1])
        src = by_id.get((orig, n))
        if src is None or src["seq"] in seen:
            continue
        seen.add(src["seq"])
        chosen.append({
            "id": f"tnfa_ctl_d{orig}_n{n}",
            "sequence": src["seq"],
            "arm": "control",
            "switch_positions": "",
            "replaced": "",
            "rmsd": src["rmsd"],
            "plddt": src["plddt"],
            "i_ptm": src["i_ptm"],
        })

    chosen = chosen[:args.limit]

    cols = ["id", "sequence", "arm", "switch_positions", "replaced",
            "rmsd", "plddt", "i_ptm"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for c in chosen:
            w.writerow(c)

    # Minimal two-column file in case the portal wants only id and sequence.
    simple = os.path.splitext(args.out)[0] + "_sequences_only.csv"
    with open(simple, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "sequence"])
        for c in chosen:
            w.writerow([c["id"], c["sequence"]])

    lengths = {len(c["sequence"]) for c in chosen}
    n_sw = sum(1 for c in chosen if c["arm"] == "switch")
    print(f"designs written: {len(chosen)}  (switch {n_sw}, control {len(chosen)-n_sw})")
    print(f"sequence lengths: {sorted(lengths)}")
    print(f"unique sequences: {len({c['sequence'] for c in chosen})}")
    print(f"non-standard characters: "
          f"{sorted({ch for c in chosen for ch in c['sequence']} - set('ACDEFGHIKLMNPQRSTVWY'))}")
    print()
    print(f"{'id':<22}{'arm':<9}{'rmsd':>8}{'plddt':>8}  substitutions")
    for c in chosen[:12]:
        print(f"{c['id']:<22}{c['arm']:<9}{c['rmsd']:>8.1f}{c['plddt']:>8.3f}  {c['replaced']}")
    print(f"\nWrote {args.out} and {simple}")


if __name__ == "__main__":
    main()
