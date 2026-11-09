"""
Decide which conformation to generate against.

The reference structure is TCR-bound, so its peptide carries the induced-fit
conformation. The assay analyte is free trimer. This superposes a predicted
free-state model onto the reference groove and measures how far the peptide
backbone has moved across the hotspot region.

Superposition uses the a1-a2 heavy-chain CA atoms only. The peptide is then
measured in place, without a second fit, so the number reported is genuine
displacement of the peptide relative to the groove rather than a best-fit
overlay of the peptide onto itself.

Chain order in a ColabFold model follows the FASTA: A heavy chain, B b2m,
C peptide. Residues are renumbered from 1 per chain, so everything here is
matched by position within the chain after an explicit sequence check.

Run from the repository root:
    .\\.venv\\Scripts\\python.exe challenges\\05-pmhc\\compare_free_state.py ^
        --predicted path\\to\\free_g12d_model_1.pdb
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import three_to_index, index_to_one

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "targets" / "g12d_pmhc.pdb"

GROOVE_RESIDUES = 182          # a1-a2, used for superposition
HOTSPOT_POSITIONS = (4, 5, 6, 7, 8)
BACKBONE = ("N", "CA", "C", "O")


def residues(structure, chain_id: str) -> list:
    chain = structure[0][chain_id]
    return [r for r in chain if r.id[0] == " "]


def sequence(chain_residues) -> str:
    return "".join(index_to_one(three_to_index(r.resname)) for r in chain_residues)


def load(path: Path, name: str):
    if not path.exists():
        sys.exit(f"missing {path}")
    return PDBParser(QUIET=True).get_structure(name, str(path))


def backbone_atoms(residue) -> list:
    return [residue[a] for a in BACKBONE if a in residue]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predicted", required=True,
                        help="ColabFold model PDB for the free complex")
    parser.add_argument("--peptide-chain", default="C")
    parser.add_argument("--heavy-chain", default="A")
    args = parser.parse_args()

    reference = load(REFERENCE, "ref")
    predicted = load(Path(args.predicted), "pred")

    ref_heavy = residues(reference, "A")[:GROOVE_RESIDUES]
    pred_heavy = residues(predicted, args.heavy_chain)[:GROOVE_RESIDUES]
    if len(pred_heavy) < GROOVE_RESIDUES:
        sys.exit(f"predicted chain {args.heavy_chain} has only {len(pred_heavy)} "
                 f"residues; is the chain order A/B/C as in the FASTA?")

    if sequence(ref_heavy) != sequence(pred_heavy):
        matches = sum(a == b for a, b in zip(sequence(ref_heavy), sequence(pred_heavy)))
        print(f"warning: groove sequences differ "
              f"({matches}/{GROOVE_RESIDUES} identical) - check chain assignment")

    ref_pep = residues(reference, "C")
    pred_pep = residues(predicted, args.peptide_chain)
    print(f"reference peptide  {sequence(ref_pep)}")
    print(f"predicted peptide  {sequence(pred_pep)}")
    if len(ref_pep) != len(pred_pep):
        sys.exit("peptide lengths differ")

    # Superpose on the groove, then move the whole predicted model with it.
    fixed = [r["CA"] for r in ref_heavy if "CA" in r]
    moving = [r["CA"] for r in pred_heavy if "CA" in r]
    sup = Superimposer()
    sup.set_atoms(fixed, moving)
    sup.apply(list(predicted.get_atoms()))
    print(f"\ngroove superposition RMSD  {sup.rms:.2f} A  ({len(fixed)} CA atoms)")

    print("\nper-residue peptide backbone displacement")
    print(f"{'pos':>4} {'res':>5} {'rmsd':>8}")
    squares, hotspot_squares, counts, hotspot_counts = 0.0, 0.0, 0, 0

    for index, (ref_res, pred_res) in enumerate(zip(ref_pep, pred_pep), start=1):
        ref_atoms = backbone_atoms(ref_res)
        pred_atoms = backbone_atoms(pred_res)
        names = {a.get_name() for a in ref_atoms} & {a.get_name() for a in pred_atoms}
        if not names:
            continue
        diffs = [np.linalg.norm(ref_res[n].coord - pred_res[n].coord) for n in names]
        value = float(np.sqrt(np.mean(np.square(diffs))))

        label = index_to_one(three_to_index(ref_res.resname)) + str(index)
        marker = " *" if index in HOTSPOT_POSITIONS else ""
        print(f"{index:>4} {label:>5} {value:>8.2f}{marker}")

        squares += sum(d * d for d in diffs)
        counts += len(diffs)
        if index in HOTSPOT_POSITIONS:
            hotspot_squares += sum(d * d for d in diffs)
            hotspot_counts += len(diffs)

    whole = float(np.sqrt(squares / counts))
    hotspot = float(np.sqrt(hotspot_squares / hotspot_counts))

    print(f"\nwhole peptide backbone RMSD   {whole:.2f} A")
    print(f"hotspot p4-p8 backbone RMSD   {hotspot:.2f} A")

    print("\nDECISION")
    if hotspot < 1.0:
        print("  < 1.0 A - induced-fit remodelling is small.")
        print("  Generate against targets/design_g12d.pdb.")
    elif hotspot > 2.0:
        print("  > 2.0 A - the free state is genuinely different.")
        print("  Trim the predicted model and generate against that instead.")
    else:
        print("  1.0-2.0 A - ambiguous.")
        print("  Split generation across both and let filtering decide.")


if __name__ == "__main__":
    main()
