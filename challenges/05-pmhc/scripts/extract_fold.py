"""
Recover on-target confidence from BoltzGen's folding step.

BoltzGen's `analysis` step can fail after `folding` has already succeeded
(assert len(sample_ids) > 0, "Computed metrics successfully for 0 out of N").
When that happens the confidence metrics look lost, but they are not: the
folding step writes one .npz per design under

    run/intermediate_designs_inverse_folded/fold_out_npz/

and those files hold the full confidence array set. This script reads them
and merges with the geometry table, which is how the metrics for all 25
designs were recovered after the analysis step crashed.

Scores reported, and why these ones:

  design_to_target_ipsae  ipSAE restricted to the design-target interface.
                          Preferred over ipTM, which is inflated by a
                          well-predicted target regardless of the interface.
  design_to_target_iptm   the same interface under ipTM, for comparison.
  min_interaction_pae     the most confident single contact in the interface;
                          a useful floor when ipSAE is near zero.
  complex_plddt           whole-assembly confidence, mostly a sanity check.

    python extract_fold.py <designs_dir> <geometry.csv> <out.csv>
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

FIELDS = [
    ("d2t_ipsae", ("design_to_target_ipsae",)),
    ("d2t_iptm", ("design_to_target_iptm",)),
    ("iptm", ("iptm",)),
    ("design_ptm", ("design_ptm",)),
    ("complex_plddt", ("complex_plddt",)),
    ("min_ipae", ("min_interaction_pae",)),
]


def scalar(data, keys):
    """Mean over the diffusion samples. Each array holds one value per
    sample (5 by default); a single sample is not the design's score."""
    for k in keys:
        if k in data:
            v = np.asarray(data[k]).flatten()
            if v.size:
                return round(float(v.mean()), 4)
    return None


def main() -> None:
    designs = Path(sys.argv[1])
    geometry = Path(sys.argv[2])
    out = Path(sys.argv[3] if len(sys.argv) > 3 else "ontarget_metrics.csv")

    npz_dir = designs / "fold_out_npz"
    if not npz_dir.is_dir():
        found = list(designs.glob("**/fold_out_npz"))
        if not found:
            sys.exit(f"no fold_out_npz under {designs}")
        npz_dir = found[0]

    folded = {}
    for f in sorted(npz_dir.glob("*.npz")):
        data = np.load(f)
        folded[f.stem] = {name: scalar(data, keys) for name, keys in FIELDS}
    print(f"recovered {len(folded)} designs from {npz_dir}")

    geo = {r["design"]: r for r in csv.DictReader(geometry.open())}
    missing = set(folded) - set(geo)
    if missing:
        print(f"  warning: no geometry row for {sorted(missing)[:5]}")

    rows = []
    for name, metrics in folded.items():
        row = {"design": name, **metrics}
        row.update({k: v for k, v in geo.get(name, {}).items() if k != "design"})
        rows.append(row)

    rows.sort(key=lambda r: -(r["d2t_ipsae"] or -1))
    fields = list(rows[0])
    with out.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"\n{'design':<18}{'ipSAE':>8}{'ipTM':>8}{'iPAE':>8}")
    for r in rows[:8]:
        print(f"{r['design']:<18}{r['d2t_ipsae']:>8.3f}{r['d2t_iptm']:>8.3f}"
              f"{r['min_ipae']:>8.2f}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
