"""
Inter-model spread of the AF2 free-state ensemble, and the trimmed design target.

The five ColabFold models are five independent samples of the free complex.
After superposing each on the reference groove, the spread of the peptide
backbone across them is a zero-cost proxy for local flexibility -- a preview
of what MD measures properly.

Also writes the trimmed free-state design target (a1-a2 + peptide).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.Polypeptide import three_to_index, index_to_one

REF = Path("targets/g12d_pmhc.pdb")
GROOVE = 182
BB = ("N", "CA", "C", "O")
A1A2_LAST = 182


class DesignSelect(Select):
    def accept_chain(self, c):
        return c.id in ("A", "C")

    def accept_residue(self, r):
        if r.id[0] != " ":
            return False
        if r.get_parent().id == "A":
            return r.id[1] <= A1A2_LAST
        return True

    def accept_atom(self, a):
        return a.element != "H" and a.altloc in (" ", "A")


def res(st, cid):
    return [r for r in st[0][cid] if r.id[0] == " "]


def main(paths: list[str], out_target: str) -> None:
    p = PDBParser(QUIET=True)
    ref = p.get_structure("ref", str(REF))
    ref_groove = res(ref, "A")[:GROOVE]
    fixed = [r["CA"] for r in ref_groove if "CA" in r]

    aligned, plddts = [], []
    for path in paths:
        st = p.get_structure("m", path)
        moving = [r["CA"] for r in res(st, "A")[:GROOVE] if "CA" in r]
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        sup.apply(list(st.get_atoms()))
        aligned.append(st)
        pep = res(st, "C")
        plddts.append([np.mean([a.bfactor for a in r]) for r in pep])

    ref_pep = res(ref, "C")
    n = len(ref_pep)
    print("peptide backbone: spread across the 5 free-state models, and mean pLDDT")
    print(f"{'pos':>4} {'res':>5} {'spread':>8} {'pLDDT':>7}")

    spreads = []
    for i in range(n):
        coords = []
        for st in aligned:
            r = res(st, "C")[i]
            coords.append(np.concatenate([r[a].coord for a in BB if a in r]))
        arr = np.vstack(coords)
        centroid = arr.mean(axis=0)
        spread = float(np.sqrt(np.mean(np.sum((arr - centroid) ** 2, axis=1) / len(BB))))
        spreads.append(spread)
        label = index_to_one(three_to_index(ref_pep[i].resname)) + str(i + 1)
        pl = np.mean([x[i] for x in plddts])
        print(f"{i+1:>4} {label:>5} {spread:>8.2f} {pl:>7.1f}")

    hs = np.mean(spreads[3:8])
    anch = np.mean([spreads[0], spreads[1], spreads[2], spreads[8], spreads[9]])
    print(f"\nmean spread, hotspot p4-p8   {hs:.2f} A")
    print(f"mean spread, anchor region   {anch:.2f} A")
    print(f"ratio                        {hs/anch:.1f}x")

    io = PDBIO()
    io.set_structure(aligned[0])
    io.save(out_target, DesignSelect())
    chk = p.get_structure("c", out_target)
    counts = {c.id: sum(1 for r in c if r.id[0] == " ") for c in chk[0]}
    print(f"\nwrote {out_target}  {counts}")


if __name__ == "__main__":
    main(sys.argv[1:-1], sys.argv[-1])
