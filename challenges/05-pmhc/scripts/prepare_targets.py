"""
Challenge 05 target preparation.

Builds the three assay targets from the reference ternary structure and
answers the question that decides the design strategy: how much of the
mutated aspartate side chain is actually reachable from solvent in the
free peptide-MHC complex, which is the state the binder meets in the assay.

Outputs
    targets/ternary.pdb      MHC + b2m + peptide + TCR, cleaned
    targets/g12d_pmhc.pdb    MHC + b2m + G12D peptide        (on-target)
    targets/wt_pmhc.pdb      MHC + b2m + wild-type peptide   (off-target 1)
    targets/apo_mhc.pdb      MHC + b2m, no peptide           (off-target 2)
    analysis/peptide_sasa.csv

Run from the repository root:
    .\\.venv\\Scripts\\python.exe challenges\\05-pmhc\\prepare_targets.py
"""

from __future__ import annotations

import csv
import sys
import urllib.request
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.SASA import ShrakeRupley

PDB_ID = "7OW6"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"

CHAIN_MHC = "A"
CHAIN_B2M = "B"
CHAIN_PEPTIDE = "C"
CHAIN_TCR_A = "D"
CHAIN_TCR_B = "E"

EXPECTED_PEPTIDE = "VVVGADGVGK"
MUTATED_RESNAME = "ASP"
GLYCINE_ATOMS = {"N", "CA", "C", "O"}

# Theoretical maximum solvent-accessible surface area per residue, Tien et al. 2013.
MAX_ASA = {"A": 129.0, "G": 104.0, "V": 174.0, "D": 193.0, "K": 236.0}

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
TARGETS = HERE / "targets"
ANALYSIS = HERE / "analysis"


class ChainSelect(Select):
    """Keep standard amino-acid residues from the named chains only."""

    def __init__(self, chains: set[str]):
        self.chains = chains

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        return residue.id[0] == " "

    def accept_atom(self, atom):
        return atom.element != "H" and (atom.altloc in (" ", "A"))


def fetch_structure() -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    path = DATA / f"{PDB_ID.lower()}.pdb"
    if not path.exists():
        print(f"downloading {PDB_ID}")
        urllib.request.urlretrieve(PDB_URL, path)
    return path


def load(path: Path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(PDB_ID, str(path))


def one_letter(residue) -> str:
    return index_to_one(three_to_index(residue.resname))


def peptide_residues(structure):
    chain = structure[0][CHAIN_PEPTIDE]
    return [r for r in chain if r.id[0] == " "]


def verify(structure) -> None:
    present = {c.id for c in structure[0]}
    required = {CHAIN_MHC, CHAIN_B2M, CHAIN_PEPTIDE, CHAIN_TCR_A, CHAIN_TCR_B}
    missing = required - present
    if missing:
        sys.exit(f"expected chains {sorted(required)}, missing {sorted(missing)}")

    residues = peptide_residues(structure)
    sequence = "".join(one_letter(r) for r in residues)
    print(f"peptide chain {CHAIN_PEPTIDE}: {sequence}")
    print("numbering: " + ", ".join(f"{one_letter(r)}{r.id[1]}" for r in residues))
    if sequence != EXPECTED_PEPTIDE:
        sys.exit(f"peptide sequence is {sequence}, expected {EXPECTED_PEPTIDE}")


def write_subset(structure, chains: set[str], out: Path) -> None:
    io = PDBIO()
    io.set_structure(structure)
    out.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out), ChainSelect(chains))


def write_wildtype(path_in: Path, path_out: Path) -> int:
    """Aspartate to glycine by deletion. Glycine's atoms are a strict subset
    of aspartate's, so no rebuilding or relaxation is required and the
    backbone geometry is preserved exactly between the two targets."""
    structure = load(path_in)
    chain = structure[0][CHAIN_PEPTIDE]
    targets = [r for r in chain if r.resname == MUTATED_RESNAME]
    if len(targets) != 1:
        sys.exit(f"found {len(targets)} aspartates in the peptide, expected exactly 1")

    residue = targets[0]
    position = [r.id[1] for r in chain if r.id[0] == " "].index(residue.id[1]) + 1
    for atom in [a for a in residue if a.get_name() not in GLYCINE_ATOMS]:
        residue.detach_child(atom.get_id())
    residue.resname = "GLY"

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(path_out))
    return position


def residue_sasa(path: Path) -> dict[int, dict[str, float]]:
    """Per-residue and side-chain-only SASA for the peptide chain."""
    structure = load(path)
    ShrakeRupley().compute(structure, level="A")

    out: dict[int, dict[str, float]] = {}
    for residue in structure[0][CHAIN_PEPTIDE]:
        if residue.id[0] != " ":
            continue
        total = sum(a.sasa for a in residue)
        side = sum(a.sasa for a in residue if a.get_name() not in GLYCINE_ATOMS)
        out[residue.id[1]] = {
            "resname": one_letter(residue),
            "total": total,
            "sidechain": side,
            "atoms": {a.get_name(): a.sasa for a in residue},
        }
    return out


def main() -> None:
    structure = load(fetch_structure())
    verify(structure)

    TARGETS.mkdir(parents=True, exist_ok=True)
    ANALYSIS.mkdir(parents=True, exist_ok=True)

    ternary = TARGETS / "ternary.pdb"
    g12d = TARGETS / "g12d_pmhc.pdb"
    apo = TARGETS / "apo_mhc.pdb"
    wt = TARGETS / "wt_pmhc.pdb"

    write_subset(structure, {CHAIN_MHC, CHAIN_B2M, CHAIN_PEPTIDE,
                             CHAIN_TCR_A, CHAIN_TCR_B}, ternary)
    write_subset(structure, {CHAIN_MHC, CHAIN_B2M, CHAIN_PEPTIDE}, g12d)
    write_subset(structure, {CHAIN_MHC, CHAIN_B2M}, apo)
    mutated_position = write_wildtype(g12d, wt)

    bound = residue_sasa(ternary)
    free = residue_sasa(g12d)

    print()
    print("peptide solvent accessibility, square angstroms")
    print(f"{'pos':>4} {'res':>4} {'free total':>11} {'free side':>10} "
          f"{'rel free':>9} {'bound side':>11} {'buried by TCR':>14}")

    rows = []
    for order, resseq in enumerate(sorted(free), start=1):
        f = free[resseq]
        b = bound[resseq]
        relative = f["total"] / MAX_ASA[f["resname"]]
        delta = b["sidechain"] - f["sidechain"]
        print(f"{order:>4} {f['resname'] + str(resseq):>4} {f['total']:>11.1f} "
              f"{f['sidechain']:>10.1f} {relative:>9.2f} {b['sidechain']:>11.1f} "
              f"{-delta:>14.1f}")
        rows.append({
            "peptide_position": order,
            "pdb_resseq": resseq,
            "residue": f["resname"],
            "free_total_sasa": round(f["total"], 2),
            "free_sidechain_sasa": round(f["sidechain"], 2),
            "free_relative_sasa": round(relative, 3),
            "tcr_bound_sidechain_sasa": round(b["sidechain"], 2),
            "sidechain_buried_by_tcr": round(-delta, 2),
        })

    with (ANALYSIS / "peptide_sasa.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    mutated = next(r for r in rows if r["peptide_position"] == mutated_position)
    resseq = mutated["pdb_resseq"]

    print()
    print(f"MUTATED RESIDUE: {mutated['residue']}{resseq}, "
          f"peptide position {mutated_position}")
    print(f"  free side-chain SASA        {mutated['free_sidechain_sasa']:.1f}")
    print(f"  free relative SASA          {mutated['free_relative_sasa']:.2f}")
    print(f"  TCR-bound side-chain SASA   {mutated['tcr_bound_sidechain_sasa']:.1f}")
    print("  per-atom, free state:")
    for name in ("CB", "CG", "OD1", "OD2"):
        value = free[resseq]["atoms"].get(name)
        if value is not None:
            print(f"    {name:<4} {value:>6.1f}")

    print()
    print("up-facing candidates, free state, ranked by side-chain SASA")
    ranked = sorted(rows, key=lambda r: r["free_sidechain_sasa"], reverse=True)
    for r in ranked:
        if r["free_sidechain_sasa"] > 5.0:
            print(f"  p{r['peptide_position']:<2} {r['residue']}{r['pdb_resseq']:<4} "
                  f"{r['free_sidechain_sasa']:>6.1f}")


if __name__ == "__main__":
    main()
