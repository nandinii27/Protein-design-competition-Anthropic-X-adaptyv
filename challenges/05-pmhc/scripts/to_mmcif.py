"""
Convert the design targets from PDB to annotated mmCIF.

The hosted generative endpoint takes an annotated mmCIF inline in the request
body, so entity and entity_poly records have to be present. Biopython's mmCIF
writer emits a minimal block without them; gemmi derives the entities from the
chains and writes a complete one.

Outputs
    targets/design_g12d.cif
    targets/design_wt.cif
    targets/design_apo.cif

Install once, then run from the repository root:
    .\\.venv\\Scripts\\python.exe -m pip install gemmi
    .\\.venv\\Scripts\\python.exe challenges\\05-pmhc\\to_mmcif.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import gemmi

HERE = Path(__file__).resolve().parent
TARGETS = HERE / "targets"
NAMES = ["design_g12d", "design_wt", "design_apo"]


def convert(stem: str) -> None:
    source = TARGETS / f"{stem}.pdb"
    if not source.exists():
        sys.exit(f"missing {source.name} - run make_design_targets.py first")

    structure = gemmi.read_structure(str(source))
    structure.setup_entities()
    structure.assign_label_seq_id()

    out = TARGETS / f"{stem}.cif"
    structure.make_mmcif_document().write_file(str(out))

    model = structure[0]
    chains = {chain.name: len(chain) for chain in model}
    size = out.stat().st_size / 1024
    print(f"{out.name:<20} {chains}  {size:.0f} KB")


def main() -> None:
    for stem in NAMES:
        convert(stem)

    print()
    print("peptide is chain C, residues 1-10")
    print("hotspots: C4 C5 C6 C7 C8")
    print("binder length: 55-90")


if __name__ == "__main__":
    main()
