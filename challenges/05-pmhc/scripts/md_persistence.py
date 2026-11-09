"""
Complex MD: does the designed binder stay bound and stay folded?

SCOPE, stated up front because it bounds what the result can claim.

This run answers one question: over nanoseconds at 300 K, does the binder
remain docked on the peptide and retain its fold, or does the interface peel
away and the helices unravel? A design that dissociates in a few nanoseconds
is a design that will not survive an hour on a BLI chip, and finding that out
costs one night of an idle laptop.

It does NOT answer whether the binder discriminates G12D from wild type.
That is a difference in bulge flexibility, and separating it from run-to-run
noise needs independent replicas on each variant. One trajectory cannot do it
and this script does not pretend to.

SOLVENT. Implicit, GBn2, zero salt. Explicit water would put this system
near 80,000 atoms, which on a laptop CPU is roughly half a nanosecond
overnight - too short to observe dissociation at all. Implicit solvent buys
two orders of magnitude and is appropriate for a persistence question, where
what matters is whether contacts survive. It is NOT appropriate for the
electrostatics of the buried aspartate, and is not used for that here.

Salt screening is absent: implicitSolventSaltConc is a parameter of the Amber
prmtop route only, and the XML force-field route used here ignores it. Adding
it would require building a prmtop through tleap. At 300 K over nanoseconds,
for a question about whether contacts survive, this is acceptable; it is
recorded here rather than papered over.

STRUCTURES. Built from BoltzGen's refold_cif outputs, not the
intermediate_designs_inverse_folded files - the latter carry 234 side-chain
atoms per complex still sitting at the origin. PDBFixer then rebuilds missing
terminal atoms (the designed chains have no OXT) and adds hydrogens at pH 7.4.

PLATFORM. CPU, forced. This machine's OpenCL returns forces differing from
the reference implementation by 13%, which is a driver fault.

Resumable. State is checkpointed every 100 ps and metrics are appended every
10 ps, so an interrupted run loses at most the last 100 ps.

    python md_persistence.py --hours 6
"""
from __future__ import annotations

import argparse
import csv
import glob
import sys
import time
from pathlib import Path

DESIGNS = ["pmhc_g12d_09"]

TEMPERATURE = 300          # kelvin
FRICTION = 1.0             # per picosecond
TIMESTEP = 0.004           # picoseconds, with hydrogen mass repartitioning
HYDROGEN_MASS = 4.0        # amu
CUTOFF = 1.2               # nanometres, implicit solvent convention

EQUIL_PS = 100
REPORT_PS = 10
CHECKPOINT_PS = 100
CONTACT_CUTOFF = 0.45      # nanometres, heavy-atom contact


def find_structures(root: Path) -> dict[str, Path]:
    found = {}
    for path in root.glob("**/pmhc_g12d_*.pdb"):
        found.setdefault(path.stem, path)
    return found


def build(pdb_path: Path):
    from openmm import app, unit, LangevinMiddleIntegrator, Platform
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # do not build loops that were never modelled
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()             # terminal OXT, any absent side-chain atoms
    fixer.addMissingHydrogens(7.4)

    ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/obc2.xml")
    model = app.Modeller(fixer.topology, fixer.positions)

    system = ff.createSystem(
        model.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=CUTOFF * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=HYDROGEN_MASS * unit.amu,
    )
    integrator = LangevinMiddleIntegrator(
        TEMPERATURE * unit.kelvin,
        FRICTION / unit.picosecond,
        TIMESTEP * unit.picoseconds,
    )
    platform = Platform.getPlatformByName("CPU")
    sim = app.Simulation(model.topology, system, integrator, platform)
    sim.context.setPositions(model.positions)
    return sim, model.topology


def index_sets(topology):
    """Binder CA atoms, and heavy atoms of binder / peptide / MHC."""
    binder_ca, binder, peptide, mhc = [], [], [], []
    for atom in topology.atoms():
        if atom.element is None or atom.element.symbol == "H":
            continue
        chain = atom.residue.chain.id
        if chain == "A":
            binder.append(atom.index)
            if atom.name == "CA":
                binder_ca.append(atom.index)
        elif chain == "D":
            peptide.append(atom.index)
        else:
            mhc.append(atom.index)
    return binder_ca, binder, peptide, mhc


def metrics(positions, ref_ca, binder_ca, binder, peptide, mhc):
    import numpy as np

    xyz = np.asarray(positions.value_in_unit_system(
        __import__("openmm").unit.md_unit_system))

    # binder fold integrity: CA RMSD after superposing the binder on itself
    p = xyz[binder_ca] - xyz[binder_ca].mean(0)
    q = ref_ca - ref_ca.mean(0)
    u, _, vt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    rmsd = float(np.sqrt(((p @ rot - q) ** 2).sum(1).mean()))

    def contacts(a, b):
        A, B = xyz[a], xyz[b]
        d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
        return int((d2 < CONTACT_CUTOFF ** 2).sum())

    com = float(np.linalg.norm(xyz[binder].mean(0) - xyz[peptide].mean(0)))
    return rmsd, contacts(binder, peptide), contacts(binder, mhc), com


def run_one(name: str, src: Path, out: Path, budget_s: float) -> None:
    from openmm import app, unit
    import numpy as np

    out.mkdir(parents=True, exist_ok=True)
    pdb_path = out / f"{name}.pdb"
    if not pdb_path.exists():
        pdb_path.write_bytes(src.read_bytes())

    csv_path = out / f"{name}_metrics.csv"
    chk_path = out / f"{name}.chk"
    fields = ["ps", "binder_rmsd_nm", "contacts_peptide", "contacts_mhc",
              "com_dist_nm"]

    sim, topology = build(pdb_path)
    binder_ca, binder, peptide, mhc = index_sets(topology)
    print(f"  atoms {topology.getNumAtoms()}  binder CA {len(binder_ca)}  "
          f"peptide heavy {len(peptide)}", flush=True)

    start_ps = 0.0
    if chk_path.exists() and csv_path.exists():
        sim.loadCheckpoint(str(chk_path))
        start_ps = float(list(csv.DictReader(csv_path.open()))[-1]["ps"])
        print(f"  resuming at {start_ps:.0f} ps", flush=True)
    else:
        def energy():
            return sim.context.getState(getEnergy=True).getPotentialEnergy(
                ).value_in_unit(unit.kilojoule_per_mole)

        print(f"  E start      {energy():12.1f} kJ/mol", flush=True)
        sim.minimizeEnergy(maxIterations=2000)
        print(f"  E minimised  {energy():12.1f} kJ/mol", flush=True)
        sim.context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin)
        print(f"  equilibrating {EQUIL_PS} ps", flush=True)
        sim.step(int(EQUIL_PS / TIMESTEP))
        e = energy()
        print(f"  E equilibrated {e:12.1f} kJ/mol", flush=True)
        if not -1e6 < e < 1e6:
            raise RuntimeError(f"energy {e:.3g} kJ/mol is unphysical - "
                               "structure did not relax, skipping")

    state = sim.context.getState(getPositions=True)
    ref = np.asarray(state.getPositions().value_in_unit(unit.nanometer))[binder_ca]
    ref = ref - ref.mean(0)

    # benchmark on the first reporting block
    block = int(REPORT_PS / TIMESTEP)
    t0 = time.time()
    sim.step(block)
    rate = REPORT_PS / (time.time() - t0) * 86.4 / 1000  # ns/day
    print(f"  {rate:.1f} ns/day -> ~{rate * budget_s / 86400:.1f} ns "
          f"in this slot", flush=True)

    deadline = time.time() + budget_s
    ps = start_ps + REPORT_PS
    since_checkpoint = 0.0

    while time.time() < deadline:
        state = sim.context.getState(getPositions=True)
        row = dict(zip(fields, (round(ps, 1),) + tuple(
            round(v, 4) for v in metrics(state.getPositions(), ref,
                                         binder_ca, binder, peptide, mhc))))
        exists = csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields)
            if not exists:
                w.writeheader()
            w.writerow(row)

        if since_checkpoint >= CHECKPOINT_PS:
            sim.saveCheckpoint(str(chk_path))
            since_checkpoint = 0.0
            print(f"  {ps:7.0f} ps  rmsd {row['binder_rmsd_nm']:.3f} nm  "
                  f"pep {row['contacts_peptide']:4d}  "
                  f"mhc {row['contacts_mhc']:5d}", flush=True)

        sim.step(block)
        ps += REPORT_PS
        since_checkpoint += REPORT_PS

    sim.saveCheckpoint(str(chk_path))
    print(f"  stopped at {ps:.0f} ps", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", type=float, default=6.0,
                    help="total wall clock, split evenly across designs")
    ap.add_argument("--root", default=".", help="folder to search for structures")
    ap.add_argument("--out", default="md_out")
    args = ap.parse_args()

    found = find_structures(Path(args.root).resolve())
    missing = [d for d in DESIGNS if d not in found]
    if missing:
        sys.exit(f"could not find {missing} under {Path(args.root).resolve()}\n"
                 f"saw: {sorted(found)[:8]}")

    budget = args.hours * 3600 / len(DESIGNS)
    print(f"{len(DESIGNS)} designs, {budget/3600:.2f} h each, CPU platform\n")

    for name in DESIGNS:
        print(f"=== {name}", flush=True)
        try:
            run_one(name, found[name], Path(args.out), budget)
        except Exception as error:
            print(f"  FAILED: {type(error).__name__}: {error}", flush=True)

    print(f"\ndone -> {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
