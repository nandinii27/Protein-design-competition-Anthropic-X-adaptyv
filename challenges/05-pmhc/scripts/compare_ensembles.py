"""
Side-by-side flexibility comparison of the G12D and wild-type free-state
AF2 ensembles.

Each set of five ColabFold models is superposed on the reference a1-a2
groove, so the peptide is measured in a common frame. Two readouts:

  spread  root-mean-square deviation of the peptide backbone across the
          five models at each position -- how much the predictions disagree
  pLDDT   mean per-residue confidence -- how certain the model is locally

Both are proxies for local conformational freedom. They are not MD, but
they are free and they are two independent views of the same question.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, Superimposer

REF = Path("targets/g12d_pmhc.pdb")
GROOVE = 182
BB = ("N", "CA", "C", "O")
P = PDBParser(QUIET=True)


def res(st, cid):
    return [r for r in st[0][cid] if r.id[0] == " "]


def analyse(paths):
    ref = P.get_structure("ref", str(REF))
    fixed = [r["CA"] for r in res(ref, "A")[:GROOVE] if "CA" in r]

    coords_per_model, plddt_per_model, seq = [], [], None
    for path in paths:
        st = P.get_structure("m", path)
        moving = [r["CA"] for r in res(st, "A")[:GROOVE] if "CA" in r]
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        sup.apply(list(st.get_atoms()))
        pep = res(st, "C")
        if seq is None:
            seq = [r.resname for r in pep]
        coords_per_model.append(
            [np.concatenate([r[a].coord for a in BB if a in r]) for r in pep]
        )
        plddt_per_model.append([np.mean([a.bfactor for a in r]) for r in pep])

    n = len(seq)
    spread, plddt = [], []
    for i in range(n):
        arr = np.vstack([m[i] for m in coords_per_model])
        c = arr.mean(axis=0)
        spread.append(float(np.sqrt(np.mean(np.sum((arr - c) ** 2, axis=1) / len(BB)))))
        plddt.append(float(np.mean([m[i] for m in plddt_per_model])))
    return seq, spread, plddt


def main(g12d_paths, wt_paths):
    seq_d, spr_d, pl_d = analyse(g12d_paths)
    seq_w, spr_w, pl_w = analyse(wt_paths)

    three_to_one = {"ALA": "A", "GLY": "G", "VAL": "V", "ASP": "D", "LYS": "K"}
    print("free-state peptide: backbone spread across 5 models, and mean pLDDT")
    print(f"{'pos':>4} {'G12D':>6} {'WT':>4} | {'spread_D':>9} {'spread_W':>9} "
          f"{'ratio':>6} | {'pLDDT_D':>8} {'pLDDT_W':>8} {'drop':>6}")
    for i in range(len(seq_d)):
        d = three_to_one.get(seq_d[i], "?")
        w = three_to_one.get(seq_w[i], "?")
        mark = "  <--" if i + 1 == 6 else ""
        print(f"{i+1:>4} {d+str(i+1):>6} {w:>4} | {spr_d[i]:>9.2f} {spr_w[i]:>9.2f} "
              f"{spr_w[i]/spr_d[i]:>6.2f} | {pl_d[i]:>8.1f} {pl_w[i]:>8.1f} "
              f"{pl_d[i]-pl_w[i]:>6.1f}{mark}")

    hs = slice(3, 8)
    anchor = [0, 1, 2, 8, 9]
    hd, hw = np.mean(spr_d[hs]), np.mean(spr_w[hs])
    ad, aw = np.mean([spr_d[i] for i in anchor]), np.mean([spr_w[i] for i in anchor])
    print(f"\n{'':20} {'G12D':>8} {'WT':>8}")
    print(f"{'mean spread p4-p8':20} {hd:>8.2f} {hw:>8.2f}   ({hw/hd:.2f}x)")
    print(f"{'mean spread anchors':20} {ad:>8.2f} {aw:>8.2f}")
    print(f"{'bulge/anchor ratio':20} {hd/ad:>8.1f} {hw/aw:>8.1f}")
    print(f"{'mean pLDDT p4-p8':20} {np.mean(pl_d[hs]):>8.1f} {np.mean(pl_w[hs]):>8.1f}")
    print(f"{'mean pLDDT anchors':20} "
          f"{np.mean([pl_d[i] for i in anchor]):>8.1f} "
          f"{np.mean([pl_w[i] for i in anchor]):>8.1f}")


if __name__ == "__main__":
    k = (len(sys.argv) - 1) // 2
    main(sys.argv[1:1 + k], sys.argv[1 + k:])
