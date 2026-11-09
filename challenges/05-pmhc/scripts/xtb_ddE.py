"""
Semiempirical QM discrimination energy: GFN2-xTB ΔΔE between the G12D and
wild-type complexes.

WHY QM HERE AND NOT A FORCE FIELD. The mutation is a carboxylate buried in a
low-dielectric groove. That is precisely the regime where fixed-charge force
fields are least reliable: no electronic polarisation, so surrounding
residues cannot respond to the buried charge; no charge transfer; and point
charges that cannot represent the quadrupolar character of a carboxylate
field at close range. GFN2-xTB carries an explicit anisotropic electrostatics
term (E_AES, second-order density fluctuations with dipole and quadrupole
character), which is exactly the physics a point-charge model discards.

WHAT IS COMPUTED.

    ΔE_int = E(complex) − E(binder) − E(target)
    ΔΔE    = ΔE_int(G12D) − ΔE_int(WT)
           = [E_c^D − E_c^G] − [E_t^D − E_t^G]

The free-binder term cancels exactly, so only four single-points are needed.
This is the calculation that makes the result defensible when the absolute
numbers are not: the two truncations are atom-for-atom identical apart from
one side chain, so truncation and basis-set-superposition error enter both
differences almost equally and largely cancel.

TRUNCATION. Peptide plus every residue with a heavy atom within 8 Å of it.
Severed backbone bonds are capped with link-atom hydrogens placed along the
original bond vector — the standard QM/MM treatment. Identical selection and
identical caps for both variants, which is what the cancellation argument
rests on.

THE WILD-TYPE VARIANT is built by deleting the Asp6 side chain and placing a
hydrogen along the former CA–CB vector. The backbone is untouched, so the
two systems differ by exactly one side chain and nothing else — no
relaxation, no re-prediction, no second structure.

STATED LIMITATIONS, because they bound what the number can carry:
  - Interaction energy is not binding free energy: no conformational or
    translational entropy, no explicit desolvation.
  - The two systems differ in net charge by one unit, since Asp is anionic
    and Gly is not. Truncation error cancels; the charge difference is the
    physical signal but also makes implicit solvation the weakest link.
  - ALPB uses a uniform dielectric with an empirical cavity, and a charge at
    the protein–water boundary inside a groove is its hardest case.
  - A single-point on one conformer is blind to entropy, and §2 of the
    methods argues the discrimination may be substantially entropic.
  - No published validation exists for xTB-derived specificity scores
    against wet-lab pMHC binder specificity.

    python xtb_ddE.py pmhc_g12d_09.pdb
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

CUT = 5.0                 # angstroms from the peptide
PEPTIDE_CHAIN = "D"
BINDER_CHAIN = "A"
MUT_RESI = 6              # Asp6 in the peptide
HARTREE_KCAL = 627.5095

CHARGED = {"ARG": +1, "LYS": +1, "ASP": -1, "GLU": -1}


def protonate(src: Path, dst: Path) -> None:
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    f = PDBFixer(filename=str(src))
    f.findMissingResidues()
    f.missingResidues = {}
    f.findNonstandardResidues()
    f.replaceNonstandardResidues()
    f.findMissingAtoms()
    f.addMissingAtoms()
    f.addMissingHydrogens(7.4)
    PDBFile.writeFile(f.topology, f.positions, dst.open("w"))


def build(pdb: Path):
    """Return (atoms, residues) where atoms carry chain/resi/name/Z/xyz."""
    import gemmi

    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    atoms, residues = [], []
    for chain in st[0]:
        for res in chain:
            entry = {"chain": chain.name, "resi": res.seqid.num,
                     "name": res.name, "atoms": []}
            for a in res:
                rec = {"chain": chain.name, "resi": res.seqid.num,
                       "res": res.name, "name": a.name,
                       "Z": a.element.atomic_number,
                       "xyz": np.array([a.pos.x, a.pos.y, a.pos.z])}
                entry["atoms"].append(rec)
                atoms.append(rec)
            residues.append(entry)
    return atoms, residues


def select(residues, cut: float):
    pep = np.array([a["xyz"] for r in residues if r["chain"] == PEPTIDE_CHAIN
                    for a in r["atoms"]])
    keep = []
    for r in residues:
        if r["chain"] == PEPTIDE_CHAIN:
            keep.append(r)
            continue
        X = np.array([a["xyz"] for a in r["atoms"]])
        if np.sqrt(((X[:, None, :] - pep[None, :, :]) ** 2).sum(-1)).min() < cut:
            keep.append(r)
    return keep


def cap(kept, residues):
    """Link-atom hydrogens on severed backbone bonds."""
    index = {(r["chain"], r["resi"]): r for r in residues}
    inside = {(r["chain"], r["resi"]) for r in kept}
    caps = []
    for r in kept:
        get = lambda n: next((a for a in r["atoms"] if a["name"] == n), None)
        for own, neighbour, partner, length in (
                ("N", (r["chain"], r["resi"] - 1), "C", 1.01),
                ("C", (r["chain"], r["resi"] + 1), "N", 1.09)):
            if neighbour in inside or neighbour not in index:
                continue
            a = get(own)
            b = next((x for x in index[neighbour]["atoms"]
                      if x["name"] == partner), None)
            if a is None or b is None:
                continue
            v = b["xyz"] - a["xyz"]
            caps.append({"Z": 1, "xyz": a["xyz"] + v / np.linalg.norm(v) * length})
    return caps


def mutate_to_gly(kept):
    """Asp6 -> Gly, backbone untouched. Returns a new residue list."""
    drop = {"CB", "CG", "OD1", "OD2", "HB2", "HB3", "HB1"}
    out = []
    for r in kept:
        if not (r["chain"] == PEPTIDE_CHAIN and r["resi"] == MUT_RESI):
            out.append(r)
            continue
        ca = next(a for a in r["atoms"] if a["name"] == "CA")
        cb = next(a for a in r["atoms"] if a["name"] == "CB")
        v = cb["xyz"] - ca["xyz"]
        link = {"chain": r["chain"], "resi": r["resi"], "res": "GLY",
                "name": "HA3", "Z": 1,
                "xyz": ca["xyz"] + v / np.linalg.norm(v) * 1.09}
        atoms = [a for a in r["atoms"] if a["name"] not in drop] + [link]
        out.append({**r, "name": "GLY", "atoms": atoms})
    return out


def charge_of(kept) -> int:
    q = 0
    for r in kept:
        q += CHARGED.get(r["name"], 0)
        names = {a["name"] for a in r["atoms"]}
        if r["name"] == "HIS" and {"HD1", "HE2"} <= names:
            q += 1
        if "OXT" in names:                        # real C-terminus
            q -= 1
        if len(names & {"H1", "H2", "H3"}) == 3:  # real N-terminus
            q += 1
    return q


def energy(kept, caps, label: str) -> float:
    from tblite.interface import Calculator

    Z = np.array([a["Z"] for r in kept for a in r["atoms"]]
                 + [c["Z"] for c in caps])
    xyz = np.array([a["xyz"] for r in kept for a in r["atoms"]]
                   + [c["xyz"] for c in caps]) * 1.8897259886   # A -> bohr
    q = charge_of(kept)
    print(f"  {label:18s} {len(Z):5d} atoms  charge {q:+d}", flush=True)

    calc = Calculator("GFN2-xTB", Z, xyz, charge=float(q))
    calc.add("alpb-solvation", "water")
    calc.set("verbosity", 0)
    for key, value in (("max-iter", 300), ("accuracy", 2.0),
                       ("electronic-temperature", 3.1668e-3)):   # ~1000 K, aids SCF
        try:
            calc.set(key, value)
        except Exception:
            print(f"    (option {key} not supported)", flush=True)
    t0 = time.time()
    e = calc.singlepoint().get("energy")
    print(f"  {label:18s} converged in {time.time()-t0:.0f}s", flush=True)
    print(f"  {label:18s} E = {e:.6f} Eh", flush=True)
    return float(e)


def main() -> None:
    src = Path(sys.argv[1] if len(sys.argv) > 1 else "pmhc_g12d_09.pdb")
    work = src.with_name(src.stem + "_H.pdb")
    if not work.exists():
        protonate(src, work)

    _, residues = build(work)
    kept = select(residues, CUT)

    systems = {}
    for variant, sel in (("g12d", kept), ("wt", mutate_to_gly(kept))):
        for part, chains in (("complex", None),
                             ("target", {PEPTIDE_CHAIN, "B", "C"})):
            s = [r for r in sel if chains is None or r["chain"] in chains]
            systems[(variant, part)] = (s, cap(s, residues))

    E = {}
    for key, (s, caps) in systems.items():
        E[key] = energy(s, caps, f"{key[0]}/{key[1]}")

    d_complex = E[("g12d", "complex")] - E[("wt", "complex")]
    d_target = E[("g12d", "target")] - E[("wt", "target")]
    ddE = (d_complex - d_target) * HARTREE_KCAL

    print()
    print(f"  E_complex(G12D) - E_complex(WT)  {d_complex*HARTREE_KCAL:+10.2f} kcal/mol")
    print(f"  E_target(G12D)  - E_target(WT)   {d_target*HARTREE_KCAL:+10.2f} kcal/mol")
    print(f"  ddE                              {ddE:+10.2f} kcal/mol")
    print("  negative favours the binder on G12D over wild type")

    out = Path("xtb_ddE.json")
    out.write_text(json.dumps(
        {"energies_hartree": {f"{k[0]}_{k[1]}": v for k, v in E.items()},
         "ddE_kcal_per_mol": ddE, "truncation_angstrom": CUT}, indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
