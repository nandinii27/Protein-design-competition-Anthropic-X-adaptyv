"""
Trim the assay targets to the design construct.

Generation cost scales with total residue count, and b2-microglobulin sits
underneath the groove contributing nothing to a binder that lands on the
peptide. Removing it and the a3 domain cuts the target from ~290 residues
to ~192, which is roughly a 40 percent speedup per trajectory on the same
GPU. This is the standard pMHC design construct.

The full trimer is kept for final scoring: anything that survives filtering
gets re-scored against the complete assay construct before submission.

Inputs  (from prepare_targets.py)
    targets/g12d_pmhc.pdb
    targets/wt_pmhc.pdb

Outputs
    targets/design_g12d.pdb   a1-a2 + G12D peptide     generation + on-target
    targets/design_wt.pdb     a1-a2 + WT peptide       counter-screen
    targets/design_apo.pdb    a1-a2, no peptide        counter-screen

Run from the repository root:
    .\\.venv\\Scripts\\python.exe challenges\\05-pmhc\\make_design_targets.py
"""

from __future__ import annotations

import sys
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select

CHAIN_MHC = "A"
CHAIN_PEPTIDE = "C"
A1A2_LAST_RESIDUE = 182  # a1 and a2 domains; a3 begins at DPPKTHM

HERE = Path(__file__).resolve().parent
TARGETS = HERE / "targets"


class DesignSelect(Select):
    def __init__(self, keep_peptide: bool):
        self.keep_peptide = keep_peptide

    def accept_chain(self, chain):
        if chain.id == CHAIN_MHC:
            return True
        return chain.id == CHAIN_PEPTIDE and self.keep_peptide

    def accept_residue(self, residue):
        if residue.id[0] != " ":
            return False
        if residue.get_parent().id == CHAIN_MHC:
            return residue.id[1] <= A1A2_LAST_RESIDUE
        return True

    def accept_atom(self, atom):
        return atom.element != "H" and atom.altloc in (" ", "A")


def trim(source: Path, out: Path, keep_peptide: bool) -> None:
    if not source.exists():
        sys.exit(f"missing {source.name} - run prepare_targets.py first")
    structure = PDBParser(QUIET=True).get_structure("t", str(source))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out), DesignSelect(keep_peptide))

    check = PDBParser(QUIET=True).get_structure("c", str(out))
    counts = {c.id: sum(1 for r in c if r.id[0] == " ") for c in check[0]}
    total = sum(counts.values())
    print(f"{out.name:<20} {counts}  total {total} residues")


def main() -> None:
    trim(TARGETS / "g12d_pmhc.pdb", TARGETS / "design_g12d.pdb", keep_peptide=True)
    trim(TARGETS / "wt_pmhc.pdb", TARGETS / "design_wt.pdb", keep_peptide=True)
    trim(TARGETS / "g12d_pmhc.pdb", TARGETS / "design_apo.pdb", keep_peptide=False)

    print()
    print("hotspots for generation: C4,C5,C6,C7,C8")
    print("binder length range:     55-90")
    print("re-score survivors against targets/g12d_pmhc.pdb (full trimer)")


if __name__ == "__main__":
    main()
