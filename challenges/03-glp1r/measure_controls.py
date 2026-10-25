"""
Measure the co-folded control complexes against the native pocket frame.

Supersedes benchmark_controls.py for the hosted-prediction data shape. Two
differences matter.

Chain identifiers. The co-folding tool takes the receptor and the peptide as
one colon-separated input and labels them by input order, so the receptor is
chain A and the peptide chain B, not R and P as in the reference structure.

Residue numbering. The co-folder numbers every chain from 1, so the predicted
receptor runs 1-409 while the reference runs 24-432. Both ranges overlap, so a
naive alignment by residue number finds hundreds of shared positions, succeeds,
and is wrong by a constant 23 throughout. The measurements that come out of
that look entirely reasonable and mean nothing. Renumbering before alignment is
not optional.

    python measure_controls.py --complexes outputs/control_complexes
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import geometry
import score
import targets

# Chain labels in the co-folder's output: receptor first, peptide second.
PREDICTED_RECEPTOR_CHAIN = "A"
PREDICTED_BINDER_CHAIN = "B"

# The predicted receptor starts at the first residue of the construct, which is
# UniProt 24. Adding this to a predicted residue number gives the UniProt number
# the reference structure uses.
NUMBERING_OFFSET = targets.UNIPROT_OFFSET - 1


def renumber(structure, chain_id: str, offset: int):
    """Shift a chain's residue numbers so they match reference numbering.

    Done in descending order so an intermediate number never collides with one
    not yet shifted.
    """
    chain = next(structure.get_models())[chain_id]
    residues = sorted(chain, key=lambda r: r.id[1], reverse=True)
    for residue in residues:
        het, number, insertion = residue.id
        residue.id = (het, number + offset, insertion)
    return structure


def match_file(directory: Path, name: str) -> Path | None:
    """Find the predicted complex for one control by fuzzy name match."""
    stem = name.lower().replace("(", "").replace(")", "").replace("-", "")
    for path in sorted(directory.rglob("*.pdb")):
        candidate = path.name.lower().replace("_", "").replace("-", "")
        if stem.replace(".", "") in candidate:
            return path
    return None


def collect(directory: Path, frame: geometry.PocketFrame,
            reference) -> pd.DataFrame:
    rows = []
    for name, record in targets.CONTROL_PEPTIDES.items():
        path = match_file(directory, name)
        if path is None:
            print(f"no predicted complex found for {name}")
            continue

        predicted = geometry.load_structure(path)
        renumber(predicted, PREDICTED_RECEPTOR_CHAIN, NUMBERING_OFFSET)
        geometry.align_to_reference(
            predicted, reference,
            PREDICTED_RECEPTOR_CHAIN, targets.RECEPTOR_CHAIN,
        )

        aligned_path = path.with_name(path.stem + "_aligned.pdb")
        from Bio.PDB import PDBIO
        io = PDBIO()
        io.set_structure(predicted)
        io.save(str(aligned_path))

        metrics = geometry.measure(
            name, aligned_path, frame,
            PREDICTED_RECEPTOR_CHAIN, PREDICTED_BINDER_CHAIN,
        )
        row = metrics.to_dict()
        row["sequence"] = record["sequence"]
        row["label"] = record["label"]
        row["source"] = "control"
        rows.append(row)
        print(f"measured {name} from {path.name}")

    return pd.DataFrame(rows)


def report(frame: pd.DataFrame) -> None:
    if frame.empty:
        print("no controls measured")
        return

    scored = score.rank(score.apply_filters(frame))

    columns = [
        "design_id", "label", "composite_score",
        "insertion_depth", "nterm_depth",
        "pocket_contact_recovery", "groove_contact_recovery",
    ]
    print("\nRanking (highest composite first)")
    print(scored[[c for c in columns if c in scored]].to_string(index=False))

    agonists = scored[scored["label"] == "agonist"]["composite_score"]
    antagonists = scored[scored["label"] == "antagonist"]["composite_score"]

    if agonists.empty or antagonists.empty:
        print("\nincomplete control set, separation not evaluable")
        return

    margin = agonists.min() - antagonists.max()
    print(f"\nworst agonist:      {agonists.min():.3f}")
    print(f"best antagonist:    {antagonists.max():.3f}")
    print(f"separation margin:  {margin:+.3f}")

    if margin > 0:
        print("\nPASS: every agonist outranks every antagonist.")
    else:
        print("\nFAIL: at least one antagonist outranks an agonist.")
        print("Before changing weights, sanity-check the two agonists'")
        print("pocket_contact_recovery. If it is near zero for both, the")
        print("predictor did not place the N-terminus in the pocket at all,")
        print("which is a prediction failure rather than a scoring failure.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--complexes", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=None)
    args = parser.parse_args()

    reference_path = args.reference or geometry.fetch_reference()
    reference = geometry.load_structure(Path(reference_path))
    frame = geometry.build_pocket_frame(reference)

    measured = collect(args.complexes, frame, reference)
    report(measured)

    if not measured.empty:
        out = args.complexes.parent / "control_metrics.csv"
        measured.to_csv(out, index=False)


if __name__ == "__main__":
    main()
