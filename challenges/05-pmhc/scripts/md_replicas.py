"""
Replica MD: does the bound binder sense the rigidity difference?

The free-state ensemble showed the wild-type p4-p7 bulge to be 2.15x more
conformationally variable than G12D, with anchor residues identical between
variants. That is the discriminating signal physically available on this
target, because the mutated aspartate itself is buried and uncontactable.

This asks whether the difference survives in the BOUND complex. Three
independent seeds on the G12D complex and three on the wild-type complex,
same binder, same starting coordinates, differing by one side chain.

WHY REPLICAS. A single trajectory gives one sample of a noisy quantity. The
claim "the wild-type bulge is more mobile" is only meaningful if the gap
between variants exceeds the spread between seeds of the same variant. Three
replicas per variant is the minimum that produces an error bar at all, and
the comparison is reported against that spread rather than as a bare ratio.

The wild type is built by deleting the Asp6 side chain from the identical
structure, so the two systems share backbone coordinates exactly and nothing
is re-predicted or relaxed between them.

RMSF is measured on peptide CA atoms after superposing each frame on the
alpha1-alpha2 groove CA atoms, so it reports peptide motion relative to the
groove rather than tumbling of the whole assembly. The anchor residues
(p1-p3, p9-p10) are reported alongside the bulge as an internal control: they
should be equally rigid in both variants, and if they are not, the
measurement is picking up something other than bulge flexibility.

Resumable, and writes one CSV row per completed replica.

    python md_replicas.py --hours 0.9 --root <dir with pmhc_g12d_09.pdb>
"""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

import numpy as np

SEEDS = (1, 2, 3)
TEMPERATURE = 300
FRICTION = 1.0
TIMESTEP = 0.004
HYDROGEN_MASS = 4.0
CUTOFF = 1.6
EQUIL_PS = 50
SAMPLE_PS = 5
CONTACT = 0.45

BULGE = (4, 5, 6, 7)
ANCHOR = (1, 2, 3, 9, 10)
MUT_RESI = 6


def make_wt(src: Path, dst: Path) -> None:
    """Asp6 -> Gly by deleting the side chain. Backbone untouched.

    Done on the PDB text so that exactly the intended atoms are dropped and
    every other coordinate in the file is byte-identical to the G12D input.
    PDBFixer then adds the glycine HA3 when it protonates.
    """
    keep = {"N", "CA", "C", "O", "OXT"}
    lines = []
    for line in src.read_text().splitlines(True):
        if line.startswith(("ATOM", "HETATM")):
            name, chain, resi = line[12:16].strip(), line[21], line[22:26].strip()
            if chain == "D" and resi == str(MUT_RESI):
                if name not in keep or line[76:78].strip() == "H":
                    continue
                line = line[:17] + "GLY" + line[20:]
        lines.append(line)
    dst.write_text("".join(lines))


def build(pdb_path: Path, seed: int):
    from openmm import LangevinMiddleIntegrator, Platform, app, unit
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/obc2.xml")
    model = app.Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(model.topology,
                             nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=CUTOFF * unit.nanometer,
                             constraints=app.HBonds,
                             hydrogenMass=HYDROGEN_MASS * unit.amu)
    integrator = LangevinMiddleIntegrator(TEMPERATURE * unit.kelvin,
                                          FRICTION / unit.picosecond,
                                          TIMESTEP * unit.picoseconds)
    integrator.setRandomNumberSeed(seed)
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            platform = Platform.getPlatformByName(name)
            break
        except Exception:
            continue
    sim = app.Simulation(model.topology, system, integrator, platform)
    sim.context.setPositions(model.positions)
    return sim, model.topology, platform.getName()


def indices(topology):
    groove, peptide, binder = [], {}, []
    for atom in topology.atoms():
        if atom.name != "CA":
            continue
        chain = atom.residue.chain.id
        if chain == "B":
            groove.append(atom.index)
        elif chain == "D":
            peptide[atom.residue.id if isinstance(atom.residue.id, int)
                    else len(peptide) + 1] = atom.index
        elif chain == "A":
            binder.append(atom.index)
    order = sorted(peptide)
    return groove, [peptide[k] for k in order], binder


def superpose(mobile, ref):
    p = mobile - mobile.mean(0)
    q = ref - ref.mean(0)
    u, _, vt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, d]) @ vt


def run(pdb_path: Path, seed: int, budget_s: float, out: Path, variant: str):
    from openmm import unit

    sim, topology, platform = build(pdb_path, seed)
    groove, peptide, binder = indices(topology)
    print(f"  {variant} seed {seed} on {platform}: "
          f"{len(groove)} groove CA, {len(peptide)} peptide CA", flush=True)

    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(TEMPERATURE * unit.kelvin, seed)
    sim.step(int(EQUIL_PS / TIMESTEP))

    state = sim.context.getState(getPositions=True)
    xyz = np.asarray(state.getPositions().value_in_unit(unit.nanometer))
    ref_groove = xyz[groove] - xyz[groove].mean(0)

    frames = []
    block = int(SAMPLE_PS / TIMESTEP)
    deadline = time.time() + budget_s
    while time.time() < deadline:
        sim.step(block)
        xyz = np.asarray(sim.context.getState(getPositions=True)
                         .getPositions().value_in_unit(unit.nanometer))
        R = superpose(xyz[groove], ref_groove)
        frames.append(((xyz[peptide] - xyz[groove].mean(0)) @ R) * 10.0)

    traj = np.asarray(frames)                       # frames x 10 x 3, angstrom
    rmsf = np.sqrt(((traj - traj.mean(0)) ** 2).sum(-1).mean(0))
    np.save(out / f"{variant}_s{seed}_peptideCA.npy", traj)

    rec = {"variant": variant, "seed": seed, "frames": len(traj),
           "ns": round(len(traj) * SAMPLE_PS / 1000, 3)}
    for i, r in enumerate(rmsf, start=1):
        rec[f"rmsf_p{i}"] = round(float(r), 4)
    rec["bulge_rmsf"] = round(float(np.mean([rmsf[i - 1] for i in BULGE])), 4)
    rec["anchor_rmsf"] = round(float(np.mean([rmsf[i - 1] for i in ANCHOR])), 4)
    rec["bulge_over_anchor"] = round(rec["bulge_rmsf"] / rec["anchor_rmsf"], 3)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", type=float, default=0.9)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", default="/kaggle/working/replicas")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    src = next(Path(args.root).glob("**/pmhc_g12d_09.pdb"))
    g12d, wt = out / "g12d.pdb", out / "wt.pdb"
    g12d.write_bytes(src.read_bytes())
    make_wt(src, wt)
    print(f"source {src}\nbuilt {g12d.name} and {wt.name}", flush=True)

    jobs = [(v, s) for s in SEEDS for v in ("g12d", "wt")]   # interleaved, so an
    budget = args.hours * 3600 / len(jobs)                   # early stop still
    print(f"{len(jobs)} replicas, {budget/60:.1f} min each\n", flush=True)

    csv_path = out / "replicas.csv"
    done = set()
    if csv_path.exists():
        done = {(r["variant"], int(r["seed"]))
                for r in csv.DictReader(csv_path.open())}

    for variant, seed in jobs:                 # leaves balanced variants
        if (variant, seed) in done:
            continue
        try:
            rec = run(out / f"{variant}.pdb", seed, budget, out, variant)
        except Exception as error:     # one bad replica must not lose the set
            print(f"  {variant} seed {seed} FAILED: "
                  f"{type(error).__name__}: {error}", flush=True)
            continue
        exists = csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=list(rec))
            if not exists:
                w.writeheader()
            w.writerow(rec)
        print("CSV," + ",".join(str(v) for v in rec.values()), flush=True)
        print(f"    {rec['ns']} ns  bulge {rec['bulge_rmsf']:.3f} A  "
              f"anchor {rec['anchor_rmsf']:.3f} A  "
              f"ratio {rec['bulge_over_anchor']:.2f}", flush=True)

    print(f"\ndone -> {csv_path}")


if __name__ == "__main__":
    main()
