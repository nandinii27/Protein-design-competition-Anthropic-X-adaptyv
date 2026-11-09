"""
Build the design and counter-screen targets from a full pMHC structure.

Produces, with chains relabelled A / B / C:

    target_g12d.pdb   alpha1-alpha2 (182) + beta2m (99) + VVVGADGVGK (10)
    target_wt.pdb     the same with Asp6 -> Gly
    target_apo.pdb    alpha1-alpha2 + beta2m, no peptide

BETA2-MICROGLOBULIN IS NOT OPTIONAL. An earlier version of this script
trimmed the target to alpha1-alpha2 plus peptide, on the reasoning that beta2m
sits far from the groove and only adds atoms. It does, and the empty socket it
leaves behind is a binding site: a generative run against that target produced
142 designs with zero peptide contact, 50 of which were more than 80%
identical to beta2m. Keep beta2m in the target and declare it off-limits in
the design spec (`not_binding: "all"`), which turns an open socket into a
negative-design constraint.

ALPHA-3 IS DROPPED. Residues 183 onward of the heavy chain form the Ig-like
alpha3 domain, which sits below the platform, contacts CD8 rather than a
groove binder, and roughly doubles the atom count for no interface.

THE WILD TYPE IS BUILT BY DELETION, NOT PREDICTION. Glycine's atoms are a
strict subset of aspartate's, so removing CB, CG, OD1 and OD2 converts the
residue while leaving every backbone coordinate identical. Re-predicting or
relaxing the wild-type structure would confound a sequence difference with a
conformational one, and the comparison between variants is the whole point.

Chains are identified by content rather than by letter, because chain naming
is not consistent between the PDB entry, the processed files and the
predicted complexes.

    python build_free_targets.py g12d_pmhc.pdb targets/
"""
from __future__ import annotations

import sys
from pathlib import Path

import gemmi

def amino_acids(chain):
    """Standard amino-acid residues of a chain, in order."""
    return [r for r in chain
            if gemmi.find_tabulated_residue(r.name)
            and gemmi.find_tabulated_residue(r.name).is_amino_acid()]


A1A2_LENGTH = 182         # alpha1 + alpha2, the peptide-binding platform
B2M_START = "IQRTPKIQVY"
HEAVY_START = "GSHSMRY"
MUT_RESI = 6
ASP_SIDECHAIN = {"CB", "CG", "OD1", "OD2"}


def classify(model):
    """Return (heavy, b2m, peptide) chains, identified by their content."""
    heavy = b2m = peptide = None
    for chain in model:
        residues = amino_acids(chain)
        if not residues:
            continue
        seq = gemmi.one_letter_code([r.name for r in residues])
        if 8 <= len(residues) <= 11:
            peptide = chain
        elif seq.startswith(B2M_START[:6]) or (95 <= len(residues) <= 103):
            b2m = chain
        elif seq.startswith(HEAVY_START[:5]) or seq[1:].startswith(HEAVY_START[:5]):
            if heavy is None or len(residues) > len(list(heavy)):
                heavy = chain
    missing = [n for n, c in (("heavy", heavy), ("b2m", b2m), ("peptide", peptide))
               if c is None]
    if missing:
        sys.exit(f"could not identify {missing}; chains present: "
                 f"{[(c.name, len(c)) for c in model]}")
    return heavy, b2m, peptide


def emit(path: Path, parts, drop_sidechain=False):
    """parts: list of (new_chain_id, source_chain). Writes a PDB."""
    st = gemmi.Structure()
    st.add_model(gemmi.Model("1"))
    for new_id, source in parts:
        chain = gemmi.Chain(new_id)
        for n, res in enumerate(amino_acids(source), start=1):
            copy = gemmi.Residue()
            copy.name = res.name
            copy.seqid = gemmi.SeqId(n, " ")
            for atom in res:
                if (drop_sidechain and new_id == "C" and n == MUT_RESI
                        and atom.name in ASP_SIDECHAIN):
                    continue
                copy.add_atom(atom)
            if drop_sidechain and new_id == "C" and n == MUT_RESI:
                copy.name = "GLY"
            chain.add_residue(copy)
        st[0].add_chain(chain)
    st.setup_entities()
    st.write_pdb(str(path))
    counts = [(c.name, len(c)) for c in st[0]]
    print(f"  {path.name:<20} {counts}")


def main() -> None:
    src = Path(sys.argv[1])
    out = Path(sys.argv[2] if len(sys.argv) > 2 else ".")
    out.mkdir(parents=True, exist_ok=True)

    st = gemmi.read_structure(str(src))
    st.remove_ligands_and_waters()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    heavy, b2m, peptide = classify(st[0])

    residues = amino_acids(heavy)
    print(f"heavy chain {heavy.name}: {len(residues)} residues "
          f"-> trimming to {A1A2_LENGTH} (alpha1-alpha2)")
    platform = gemmi.Chain(heavy.name)
    for r in residues[:A1A2_LENGTH]:
        platform.add_residue(r)

    pep_seq = gemmi.one_letter_code([r.name for r in amino_acids(peptide)])
    print(f"peptide {peptide.name}: {pep_seq}")

    emit(out / "target_g12d.pdb", [("A", platform), ("B", b2m), ("C", peptide)])
    emit(out / "target_wt.pdb", [("A", platform), ("B", b2m), ("C", peptide)],
         drop_sidechain=True)
    emit(out / "target_apo.pdb", [("A", platform), ("B", b2m)])


if __name__ == "__main__":
    main()
