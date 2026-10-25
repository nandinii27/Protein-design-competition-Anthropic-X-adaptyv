"""
Fused Gromov-Wasserstein comparison of design interfaces against the native agonist.

The question this answers is not "is this a confident interface" - co-folding
confidence already answers that, and the control set shows it answers it
identically for agonists and antagonists. The question is "does this design
engage the receptor the way the native agonist does", which is a comparison
between two attributed graphs with no correspondence between them.

There is no alignment available. The native hormone is 30 residues, a design
may be 45, and no residue in one corresponds to any particular residue in the
other. Fused Gromov-Wasserstein is built for exactly this: it couples a Gromov
term, which matches intra-graph geometry without needing shared coordinates,
to a Wasserstein term, which matches node attributes. The result is a single
distance between a design's interface and the native agonist's.

Each interface is represented as:

  C   pairwise alpha-carbon distances among the binder residues that contact
      the receptor, normalised by the graph's own maximum so that graphs of
      different size and extent remain comparable
  F   per-residue features: normalised depth along the pocket axis, hydropathy,
      and formal charge. Depth is the feature that carries the mechanism - it
      says where in the receptor that residue sits, which is the whole
      distinction between a groove binder and an agonist
  p   uniform mass over nodes

CALIBRATION WARNING. The alpha parameter trades structure against features and
there are four control peptides. Four points can be separated by almost any
metric under almost any parameter, so tuning alpha until the controls split is
fitting, not evidence. Leave it at the default, report the controls' distances
as a sanity reading rather than a validation, and treat this as an exploratory
ranking signal alongside insertion depth - not as a replacement for either.

    python fgw_interface.py --complexes outputs/complexa --out outputs/fgw.csv
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

import geometry
import score
import targets

try:
    import ot
except ImportError as exc:
    raise SystemExit("POT is required: pip install pot") from exc

CONTACT_CUTOFF_ANGSTROM = 5.0

# Trade-off between the structural (Gromov) and feature (Wasserstein) terms.
# In POT's convention the objective is (1 - alpha) * feature + alpha * structure,
# so 1.0 is pure geometry and 0.0 is pure attributes. Half weights both. Do not
# tune this on the control set.
DEFAULT_ALPHA = 0.5

FORMAL_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.5}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass
class InterfaceGraph:
    structure_matrix: np.ndarray   # normalised pairwise distances, n x n
    features: np.ndarray           # n x 3
    mass: np.ndarray               # uniform, length n
    size: int

    @property
    def is_empty(self) -> bool:
        return self.size < 2


def _residue_letter(residue) -> str:
    return THREE_TO_ONE.get(residue.get_resname().upper(), "X")


def build_interface_graph(structure,
                          receptor_chain_id: str,
                          binder_chain_id: str,
                          frame: geometry.PocketFrame) -> InterfaceGraph:
    """Represent the binder side of an interface as an attributed graph.

    Only binder residues that actually touch the receptor are kept. A design
    whose tail dangles into solvent should not have that tail counted as part
    of how it engages the receptor.
    """
    model = next(structure.get_models())
    receptor = model[receptor_chain_id]
    binder = model[binder_chain_id]

    receptor_atoms = np.asarray(
        [a.get_coord() for r in receptor for a in r if a.element != "H"]
    )

    contact_residues, contact_ca = [], []
    for residue in binder:
        if "CA" not in residue:
            continue
        atoms = np.asarray([a.get_coord() for a in residue if a.element != "H"])
        if atoms.size == 0:
            continue
        distances = np.linalg.norm(
            receptor_atoms[:, None, :] - atoms[None, :, :], axis=-1
        )
        if distances.min() < CONTACT_CUTOFF_ANGSTROM:
            contact_residues.append(residue)
            contact_ca.append(residue["CA"].get_coord())

    if len(contact_residues) < 2:
        return InterfaceGraph(np.zeros((0, 0)), np.zeros((0, 3)),
                              np.zeros(0), len(contact_residues))

    coordinates = np.asarray(contact_ca)

    structure_matrix = cdist(coordinates, coordinates)
    largest = structure_matrix.max()
    if largest > 0:
        structure_matrix = structure_matrix / largest

    features = []
    for residue, coordinate in zip(contact_residues, coordinates):
        letter = _residue_letter(residue)
        features.append([
            frame.depth(coordinate),
            score.HYDROPATHY.get(letter, 0.0) / 4.5,
            FORMAL_CHARGE.get(letter, 0.0),
        ])
    features = np.asarray(features, dtype=float)

    n = len(contact_residues)
    return InterfaceGraph(structure_matrix, features, np.full(n, 1.0 / n), n)


def fgw_distance(a: InterfaceGraph, b: InterfaceGraph,
                 alpha: float = DEFAULT_ALPHA) -> float:
    """Fused Gromov-Wasserstein distance between two interface graphs.

    Returns NaN rather than a number when either graph is too small to compare,
    so that a design with essentially no interface cannot be ranked as if it
    had been measured.
    """
    if a.is_empty or b.is_empty:
        return float("nan")

    feature_cost = cdist(a.features, b.features)

    return float(ot.gromov.fused_gromov_wasserstein2(
        feature_cost,
        a.structure_matrix,
        b.structure_matrix,
        a.mass,
        b.mass,
        loss_fun="square_loss",
        alpha=alpha,
    ))


def native_graph(reference, frame: geometry.PocketFrame) -> InterfaceGraph:
    """The reference point: the native agonist's own interface."""
    return build_interface_graph(
        reference, targets.RECEPTOR_CHAIN, targets.REFERENCE_PEPTIDE_CHAIN, frame
    )


def guess_binder_chain(structure) -> str:
    """The binder is the shorter chain. Design outputs do not use fixed labels."""
    model = next(structure.get_models())
    chains = [(c.id, sum(1 for r in c if "CA" in r)) for c in model]
    chains = [c for c in chains if c[1] > 0]
    if len(chains) < 2:
        raise ValueError("expected at least two chains in the complex")
    return min(chains, key=lambda c: c[1])[0]


def guess_receptor_chain(structure) -> str:
    model = next(structure.get_models())
    chains = [(c.id, sum(1 for r in c if "CA" in r)) for c in model]
    return max(chains, key=lambda c: c[1])[0]


def score_directory(directory: Path,
                    frame: geometry.PocketFrame,
                    native: InterfaceGraph,
                    alpha: float) -> pd.DataFrame:
    rows = []
    for path in sorted(directory.rglob("*.pdb")):
        if path.stem.endswith("_aligned"):
            continue
        try:
            complex_structure = geometry.load_structure(path)
            receptor_chain = guess_receptor_chain(complex_structure)
            binder_chain = guess_binder_chain(complex_structure)
            graph = build_interface_graph(
                complex_structure, receptor_chain, binder_chain, frame
            )
            distance = fgw_distance(graph, native, alpha)
        except Exception as error:
            print(f"skipped {path.name}: {error}")
            continue

        rows.append({
            "design_id": path.stem,
            "fgw_distance": distance,
            "interface_size": graph.size,
            "mean_interface_depth": (
                float(np.mean(graph.features[:, 0])) if graph.size else float("nan")
            ),
        })
        print(f"{path.stem}  fgw={distance:.4f}  n={graph.size}")

    frame_out = pd.DataFrame(rows)
    if not frame_out.empty:
        # Invert into a similarity in [0, 1] so it can join the composite score,
        # where every other term is higher-is-better.
        finite = frame_out["fgw_distance"].replace([np.inf, -np.inf], np.nan)
        largest = finite.max()
        if largest and largest > 0:
            frame_out["fgw_similarity"] = 1.0 - finite / largest
    return frame_out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--complexes", type=Path, required=True,
                        help="directory of design-receptor complex PDBs")
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("--out", type=Path, default=Path("outputs/fgw.csv"))
    args = parser.parse_args()

    reference_path = args.reference or geometry.fetch_reference()
    reference = geometry.load_structure(Path(reference_path))
    frame = geometry.build_pocket_frame(reference)
    native = native_graph(reference, frame)

    print(f"native interface: {native.size} contacting residues")
    if native.is_empty:
        raise SystemExit("native interface graph is empty; check chain ids")

    measured = score_directory(args.complexes, frame, native, args.alpha)
    if measured.empty:
        raise SystemExit("no complexes measured")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    measured.sort_values("fgw_distance").to_csv(args.out, index=False)
    print(f"\n{len(measured)} complexes measured")
    print("lowest distance is most similar to the native agonist interface")


if __name__ == "__main__":
    main()
