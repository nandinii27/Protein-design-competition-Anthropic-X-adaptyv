"""
Geometric measurement of pocket insertion for designed complexes.

This is the part of the pipeline that attempts to measure the second objective.
No hosted scoring tool outputs it. Co-folding confidence (ipTM, pAE) reports
whether a pose is plausible, not whether the design reaches the part of the
receptor that triggers activation.

The measurement rests on the two-domain mechanism of class B1 receptors: the
C-terminal half of the native peptide binds the extracellular domain groove and
supplies affinity, while the N-terminal residues insert into the transmembrane
bundle core and supply activation. A design that reproduces only the first
contact set is an antagonist. So the axis running from the native peptide's
N-terminal alpha carbon to its C-terminal alpha carbon is a one-dimensional
coordinate for "how deep into the receptor does this design reach", and native
agonist depth is the reference point.

Two quantities come out of it:

  insertion_depth    projection of the design's deepest alpha carbon onto the
                     pocket axis, normalised so 1.0 equals native depth and 0.0
                     equals the groove entrance. Values near 0 indicate an
                     ECD-only binder, which is the antagonist architecture.

  contact_recovery   fraction of the native peptide's receptor contacts that
                     the design reproduces, reported separately for the groove
                     contacts and the pocket contacts. Pocket recovery is the
                     quantity that matters for the activation objective.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBList
from Bio.PDB.Structure import Structure

import targets

CONTACT_CUTOFF_ANGSTROM = 5.0

# Residues of the native peptide treated as the activation module. In GLP-1
# numbering these are His7 through Thr13, the segment that penetrates the
# transmembrane bundle. Indices here are zero-based into the bound peptide
# chain as it appears in the reference structure, which may be truncated at the
# N-terminus relative to the full hormone. Check the first modelled residue
# before trusting these.
ACTIVATION_SPAN = slice(0, 7)


@dataclass
class InsertionMetrics:
    design_id: str
    insertion_depth: float
    nterm_depth: float
    pocket_contact_recovery: float
    groove_contact_recovery: float
    n_pocket_contacts: int
    n_groove_contacts: int

    def to_dict(self) -> dict:
        return asdict(self)


def fetch_reference(pdb_id: str = targets.REFERENCE_PDB_ID,
                    directory: Path = Path("structures")) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    downloader = PDBList()
    path = downloader.retrieve_pdb_file(
        pdb_id, pdir=str(directory), file_format="mmCif"
    )
    return Path(path)


def load_structure(path: Path) -> Structure:
    suffix = path.suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def _ca_coordinates(chain) -> tuple[np.ndarray, list]:
    coordinates, residues = [], []
    for residue in chain:
        if "CA" in residue:
            coordinates.append(residue["CA"].get_coord())
            residues.append(residue)
    return np.asarray(coordinates), residues


def _all_atom_coordinates(chain) -> tuple[np.ndarray, list]:
    coordinates, owners = [], []
    for residue in chain:
        for atom in residue:
            if atom.element == "H":
                continue
            coordinates.append(atom.get_coord())
            owners.append(residue)
    return np.asarray(coordinates), owners


@dataclass
class PocketFrame:
    """The coordinate frame derived from the native agonist pose."""

    origin: np.ndarray          # native peptide N-terminal CA, the deepest point
    axis: np.ndarray            # unit vector pointing outwards towards the ECD
    native_span: float          # distance from pocket floor to groove entrance
    pocket_contacts: set
    groove_contacts: set

    def depth(self, coordinate: np.ndarray) -> float:
        """Normalised depth: 1.0 at the native N-terminal position, 0.0 at the
        groove entrance, negative outside the receptor."""
        projection = float(np.dot(coordinate - self.origin, self.axis))
        return 1.0 - projection / self.native_span


def build_pocket_frame(structure: Structure,
                       receptor_chain_id: str = targets.RECEPTOR_CHAIN,
                       peptide_chain_id: str = targets.REFERENCE_PEPTIDE_CHAIN
                       ) -> PocketFrame:
    model = next(structure.get_models())
    receptor = model[receptor_chain_id]
    peptide = model[peptide_chain_id]

    peptide_ca, peptide_residues = _ca_coordinates(peptide)
    if len(peptide_ca) < 4:
        raise ValueError(
            f"peptide chain {peptide_chain_id} has too few residues; check the chain id"
        )

    origin = peptide_ca[0]
    outward = peptide_ca[-1] - origin
    native_span = float(np.linalg.norm(outward))
    axis = outward / native_span

    receptor_atoms, receptor_owners = _all_atom_coordinates(receptor)
    peptide_atoms, peptide_owners = _all_atom_coordinates(peptide)

    distances = np.linalg.norm(
        receptor_atoms[:, None, :] - peptide_atoms[None, :, :], axis=-1
    )
    close = np.argwhere(distances < CONTACT_CUTOFF_ANGSTROM)

    activation_residues = set(id(r) for r in peptide_residues[ACTIVATION_SPAN])

    pocket_contacts, groove_contacts = set(), set()
    for receptor_index, peptide_index in close:
        receptor_residue = receptor_owners[receptor_index]
        peptide_residue = peptide_owners[peptide_index]
        key = receptor_residue.get_id()[1]
        if id(peptide_residue) in activation_residues:
            pocket_contacts.add(key)
        else:
            groove_contacts.add(key)

    # A receptor residue contacted by both segments is counted as pocket, since
    # pocket contact is the scarcer and more informative signal.
    groove_contacts -= pocket_contacts

    return PocketFrame(
        origin=origin,
        axis=axis,
        native_span=native_span,
        pocket_contacts=pocket_contacts,
        groove_contacts=groove_contacts,
    )


def reference_contacts(frame: PocketFrame) -> dict[str, list[int]]:
    """Receptor residue numbers to use as hotspots in the generation tools."""
    return {
        "pocket": sorted(frame.pocket_contacts),
        "groove": sorted(frame.groove_contacts),
    }


def measure(design_id: str,
            complex_path: Path,
            frame: PocketFrame,
            receptor_chain_id: str,
            binder_chain_id: str) -> InsertionMetrics:
    """Measure a predicted design-receptor complex against the native frame.

    The predicted complex must be superposed on the reference structure before
    this is called, otherwise the frame does not apply. align_to_reference()
    below does that.
    """
    structure = load_structure(complex_path)
    model = next(structure.get_models())
    receptor = model[receptor_chain_id]
    binder = model[binder_chain_id]

    binder_ca, _ = _ca_coordinates(binder)
    depths = np.array([frame.depth(coordinate) for coordinate in binder_ca])

    receptor_atoms, receptor_owners = _all_atom_coordinates(receptor)
    binder_atoms, _ = _all_atom_coordinates(binder)
    distances = np.linalg.norm(
        receptor_atoms[:, None, :] - binder_atoms[None, :, :], axis=-1
    )
    contacted = {
        receptor_owners[index].get_id()[1]
        for index in np.unique(np.argwhere(distances < CONTACT_CUTOFF_ANGSTROM)[:, 0])
    }

    pocket_recovery = (
        len(contacted & frame.pocket_contacts) / len(frame.pocket_contacts)
        if frame.pocket_contacts else 0.0
    )
    groove_recovery = (
        len(contacted & frame.groove_contacts) / len(frame.groove_contacts)
        if frame.groove_contacts else 0.0
    )

    return InsertionMetrics(
        design_id=design_id,
        insertion_depth=float(depths.max()),
        nterm_depth=float(depths[0]),
        pocket_contact_recovery=pocket_recovery,
        groove_contact_recovery=groove_recovery,
        n_pocket_contacts=len(contacted & frame.pocket_contacts),
        n_groove_contacts=len(contacted & frame.groove_contacts),
    )


def align_to_reference(mobile: Structure, reference: Structure,
                       mobile_chain: str, reference_chain: str) -> Structure:
    """Superpose a predicted complex onto the reference by receptor backbone.

    Only residues present in both chains at the same sequence position are used,
    so a predicted construct that is truncated relative to the reference still
    aligns.
    """
    from Bio.PDB import Superimposer

    mobile_model = next(mobile.get_models())
    reference_model = next(reference.get_models())

    mobile_map = {
        residue.get_id()[1]: residue["CA"]
        for residue in mobile_model[mobile_chain]
        if "CA" in residue
    }
    reference_map = {
        residue.get_id()[1]: residue["CA"]
        for residue in reference_model[reference_chain]
        if "CA" in residue
    }
    shared = sorted(set(mobile_map) & set(reference_map))
    if len(shared) < 20:
        raise ValueError(
            f"only {len(shared)} shared residues between receptor chains; "
            "check chain identifiers and residue numbering"
        )

    superimposer = Superimposer()
    superimposer.set_atoms(
        [reference_map[i] for i in shared], [mobile_map[i] for i in shared]
    )
    superimposer.apply(list(mobile.get_atoms()))
    return mobile
