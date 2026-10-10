"""
Rank designed binder-target complexes by predicted pH-selectivity.

Computes the binding-linked proton release from PROPKA pKa values in the bound
and free states, then integrates the Wyman linkage relation over the assay
window to predict how much binding changes between pH 7.4 and pH 6.0.

    d log Ka / d pH = -dNu(pH)
    dNu(pH) = sum_i [ f_i_bound(pH) - f_i_free(pH) ]
    f_i(pH) = 1 / (1 + 10^(pH - pKa_i))

A design that binds at 7.4 and releases at 6.0 needs log Ka(7.4) > log Ka(6.0),
which requires net proton release on binding across the window.

Usage:
    python linkage_score.py --complexes designs/*.pdb --binder-chain B --out ranked.csv

Dependencies:
    pip install propka biopython
"""

import argparse
import glob
import math
import os
import subprocess
import sys
import tempfile
import warnings

from Bio.PDB import PDBParser, PDBIO, Select

warnings.filterwarnings("ignore")

RT_LN10 = 1.364  # kcal/mol at 298 K
TITRATABLE_CODES = {"HIS", "ASP", "GLU", "CYS", "TYR", "LYS", "ARG", "N+", "C-"}


class ChainSelect(Select):
    def __init__(self, chains):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains


def write_subset(structure, chains, path):
    io_writer = PDBIO()
    io_writer.set_structure(structure)
    io_writer.save(path, ChainSelect(chains))


def run_propka(pdb_path, workdir):
    """Run propka3 and return {(chain, resnum, resname): pKa}."""
    subprocess.run(
        [sys.executable, "-m", "propka", os.path.basename(pdb_path)],
        cwd=workdir, check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    pka_path = os.path.splitext(pdb_path)[0] + ".pka"
    if not os.path.exists(pka_path):
        raise RuntimeError(f"PROPKA produced no output for {pdb_path}")

    values = {}
    in_summary = False
    with open(pka_path) as fh:
        for line in fh:
            if line.startswith("SUMMARY OF THIS PREDICTION"):
                in_summary = True
                continue
            if in_summary:
                if line.strip().startswith("-") or not line.strip():
                    if values:
                        break
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                resname, resnum, chain, pka = parts[0], parts[1], parts[2], parts[3]
                if resname not in TITRATABLE_CODES:
                    continue
                try:
                    values[(chain, int(resnum), resname)] = float(pka)
                except ValueError:
                    continue
    return values


def fraction_protonated(pka, ph):
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def delta_nu(bound, free, ph):
    """Net protons bound on complex formation at a given pH."""
    total = 0.0
    for key, pka_bound in bound.items():
        pka_free = free.get(key)
        if pka_free is None:
            continue
        total += fraction_protonated(pka_bound, ph) - fraction_protonated(pka_free, ph)
    return total


def integrate_log_ka(bound, free, ph_low, ph_high, steps=140):
    """
    log Ka(ph_high) - log Ka(ph_low) = -integral of dNu dpH over the window.
    Positive result means tighter binding at the higher pH.
    """
    step = (ph_high - ph_low) / steps
    acc = 0.0
    for i in range(steps):
        a = ph_low + i * step
        b = a + step
        acc += 0.5 * (delta_nu(bound, free, a) + delta_nu(bound, free, b)) * step
    return -acc


def effective_sites(bound, free, ph_low, ph_high):
    """Count titratable sites whose protonation changes materially across the window."""
    count = 0
    for key, pka_bound in bound.items():
        pka_free = free.get(key)
        if pka_free is None:
            continue
        swing = abs(
            (fraction_protonated(pka_bound, ph_low) - fraction_protonated(pka_free, ph_low))
            - (fraction_protonated(pka_bound, ph_high) - fraction_protonated(pka_free, ph_high))
        )
        if swing > 0.25:
            count += 1
    return count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexes", nargs="+", required=True,
                    help="PDB files of binder-target complexes (globs accepted)")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--target-chains", default="A")
    ap.add_argument("--ph-low", type=float, default=6.0)
    ap.add_argument("--ph-high", type=float, default=7.4)
    ap.add_argument("--out", default="ranked.csv")
    args = ap.parse_args()

    paths = []
    for pattern in args.complexes:
        paths.extend(sorted(glob.glob(pattern)))
    if not paths:
        sys.exit("No input structures matched.")

    binder_chains = [c for c in args.binder_chain.split(",") if c]
    target_chains = [c for c in args.target_chains.split(",") if c]
    all_chains = binder_chains + target_chains

    parser = PDBParser(QUIET=True)
    results = []

    for path in paths:
        name = os.path.basename(path)
        try:
            structure = parser.get_structure("cplx", path)
            with tempfile.TemporaryDirectory() as workdir:
                cplx_path = os.path.join(workdir, "complex.pdb")
                binder_path = os.path.join(workdir, "binder.pdb")
                target_path = os.path.join(workdir, "target.pdb")

                write_subset(structure, all_chains, cplx_path)
                write_subset(structure, binder_chains, binder_path)
                write_subset(structure, target_chains, target_path)

                bound = run_propka(cplx_path, workdir)
                free_binder = run_propka(binder_path, workdir)
                free_target = run_propka(target_path, workdir)

            free = dict(free_target)
            free.update(free_binder)

            dlog = integrate_log_ka(bound, free, args.ph_low, args.ph_high)
            fold = 10.0 ** dlog
            ddg = -RT_LN10 * dlog
            n_sites = effective_sites(bound, free, args.ph_low, args.ph_high)
            nu_high = delta_nu(bound, free, args.ph_high)
            nu_low = delta_nu(bound, free, args.ph_low)

            results.append({
                "design": name,
                "delta_log_Ka": round(dlog, 3),
                "fold_selectivity": round(fold, 1),
                "ddG_kcal_per_mol": round(ddg, 3),
                "effective_sites": n_sites,
                "delta_nu_at_high_pH": round(nu_high, 3),
                "delta_nu_at_low_pH": round(nu_low, 3),
            })
        except Exception as exc:
            results.append({
                "design": name,
                "delta_log_Ka": "",
                "fold_selectivity": "",
                "ddG_kcal_per_mol": "",
                "effective_sites": "",
                "delta_nu_at_high_pH": "",
                "delta_nu_at_low_pH": f"failed: {exc}",
            })

    scored = [r for r in results if r["delta_log_Ka"] != ""]
    scored.sort(key=lambda r: r["delta_log_Ka"], reverse=True)
    failed = [r for r in results if r["delta_log_Ka"] == ""]
    ordered = scored + failed

    header = ["design", "delta_log_Ka", "fold_selectivity", "ddG_kcal_per_mol",
              "effective_sites", "delta_nu_at_high_pH", "delta_nu_at_low_pH"]
    with open(args.out, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in ordered:
            fh.write(",".join(str(r[c]) for c in header) + "\n")

    print(f"Scored {len(scored)} of {len(results)} structures.")
    print()
    print(f"{'design':<42}{'fold':>10}{'sites':>8}")
    for r in scored[:20]:
        print(f"{r['design']:<42}{r['fold_selectivity']:>10}{r['effective_sites']:>8}")
    if failed:
        print(f"\n{len(failed)} structures failed; see the output file.")


if __name__ == "__main__":
    main()
