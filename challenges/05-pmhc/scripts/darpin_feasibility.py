"""
Can a DARPin reach the peptide?

The question §9 leaves open is geometric, not thermodynamic. A DARPin's
paratope is an elongated, shallow concave surface formed by the randomised
positions of its internal repeats. The pMHC epitope is a peptide lying in a
canyon between two tall helices. Whether that shape can reach that target is
a steric question with a numerical answer, and it is worth having before
committing a design campaign to the scaffold.

METHOD. Rigid-body scan. The DARPin paratope is oriented toward the peptide
and swept through spin and tilt, pushed along the approach axis until first
hard contact, and scored by peptide heavy-atom contacts. The best placements
are then re-scored by buried solvent-accessible surface on the peptide, which
is the same quantity reported for the submitted designs, so the two are
directly comparable.

The paratope is derived from the structure rather than assumed: the three
internal repeats are located by their conserved TPLHLAA motif, aligned, and
the positions that differ between them are the randomised library positions.

CONTROL. The identical scan is run with a submitted de novo binder as the
ligand. If the scan cannot recover a good placement for a binder that is
known to bind, the scan is wrong and its DARPin answer means nothing.

WHAT THIS CANNOT SAY. It is a steric feasibility test on one fixed scaffold
conformation with no side-chain or loop relaxation, and no energy model. It
answers "can this shape reach", not "would it bind". Real DARPin selection
remodels the paratope against the target.

    python darpin_feasibility.py 1mj0.pdb pmhc_g12d_09.pdb
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial import cKDTree

CLASH = 3.3          # angstroms, heavy-heavy, first hard contact
CONTACT = 4.5        # angstroms, heavy-heavy, counted as contact
SPINS = 24           # 15 degree steps about the approach axis
TILTS = (-30, -15, 0, 15, 30)
OFFSETS = (-4.0, 0.0, 4.0)   # lateral displacement off the groove axis; a
                             # binder need not approach through the peptide
                             # centroid, and the control fails badly without
                             # this freedom


def load(path: Path):
    st = gemmi.read_structure(str(path))
    st.remove_ligands_and_waters()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    return st[0]


RADIUS = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}


def coords(model, chains=None, residues=None, sidechain_only=False, radii=False):
    out, rad = [], []
    backbone = {"N", "CA", "C", "O"}
    for ch in model:
        if chains and ch.name not in chains:
            continue
        for res in ch:
            if residues and res.seqid.num not in residues:
                continue
            for a in res:
                if sidechain_only and a.name in backbone:
                    continue
                out.append([a.pos.x, a.pos.y, a.pos.z])
                rad.append(RADIUS.get(a.element.name, 1.70))
    return (np.asarray(out), np.asarray(rad)) if radii else np.asarray(out)


def sasa(xyz, radii, n_sphere=200, probe=1.4):
    """Shrake-Rupley. Returns per-atom accessible surface in square angstroms."""
    i = np.arange(n_sphere) + 0.5
    phi = np.arccos(1 - 2 * i / n_sphere)
    theta = np.pi * (1 + 5 ** 0.5) * i
    unit = np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)

    r = radii + probe
    tree = cKDTree(xyz)
    out = np.zeros(len(xyz))
    for k in range(len(xyz)):
        pts = xyz[k] + unit * r[k]
        near = [j for j in tree.query_ball_point(xyz[k], r[k] + r.max())
                if j != k]
        if near:
            d = np.linalg.norm(pts[:, None, :] - xyz[near][None, :, :], axis=-1)
            free = (d >= r[near][None, :]).all(axis=1)
        else:
            free = np.ones(n_sphere, bool)
        out[k] = 4 * np.pi * r[k] ** 2 * free.mean()
    return out


def paratope_residues(model) -> list[int]:
    ch = next(iter(model))
    seq = gemmi.one_letter_code([r.name for r in ch])
    nums = [r.seqid.num for r in ch]
    starts = [m.start() - 6 for m in re.finditer("TPLHLAA", seq)]
    reps = [seq[s:s + 33] for s in starts]
    var = [i for i in range(33) if len({r[i] for r in reps}) > 1]
    return sorted({nums[s + i] for s in starts for i in var})


def rotation(axis, angle):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def align_to(a, b):
    """Rotation taking unit vector a onto unit vector b."""
    v = np.cross(a, b)
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if a @ b > 0 else -np.eye(3)
    return rotation(v, np.arccos(np.clip(a @ b, -1, 1)))


def scan(ligand, face_idx, target_tree, peptide_tree, axis, origin):
    """Sweep orientations; push each to first contact; score peptide contacts."""
    centre = ligand[face_idx].mean(0)
    normal = centre - ligand.mean(0)
    normal /= np.linalg.norm(normal)
    base = align_to(normal, -axis) @ (ligand - centre).T

    perp = np.cross(axis, [0, 0, 1.0])
    perp = perp / np.linalg.norm(perp)
    perp2 = np.cross(axis, perp)

    best = None
    lateral = [perp * a + perp2 * b for a in OFFSETS for b in OFFSETS]
    for spin in range(SPINS):
        Rs = rotation(axis, 2 * np.pi * spin / SPINS)
        for t1 in TILTS:
            for t2 in TILTS:
                R = Rs @ rotation(perp, np.radians(t1)) @ rotation(perp2, np.radians(t2))
                moved = (R @ base).T
                for shift in lateral:
                    here = origin + shift
                    lo, hi = 0.0, 40.0      # push along the axis to first contact
                    for _ in range(20):
                        mid = (lo + hi) / 2
                        if target_tree.query(moved + here + axis * mid,
                                             k=1)[0].min() < CLASH:
                            lo = mid
                        else:
                            hi = mid
                    pts = moved + here + axis * hi
                    n = sum(len(x) for x in
                            peptide_tree.query_ball_point(pts, CONTACT))
                    if best is None or n > best[0]:
                        best = (n, pts)
    return best


def buried_peptide_surface(target_model, ligand_pts, ligand_radii) -> float:
    """Peptide surface lost when the ligand is placed: free minus bound."""
    pep, pep_r = coords(target_model, chains={"D"}, radii=True)
    rest, rest_r = coords(target_model, chains={"B", "C"}, radii=True)
    n = len(pep)

    free = sasa(np.vstack([pep, rest]), np.concatenate([pep_r, rest_r]))[:n].sum()
    bound = sasa(np.vstack([pep, rest, ligand_pts]),
                 np.concatenate([pep_r, rest_r, ligand_radii]))[:n].sum()
    return max(0.0, free - bound)


def main() -> None:
    darpin_path = Path(sys.argv[1])
    target_path = Path(sys.argv[2])

    darpin = load(darpin_path)
    complex_model = load(target_path)

    para = paratope_residues(darpin)
    print(f"paratope: {len(para)} residues  {para}")

    dar_all, dar_r = coords(darpin, radii=True)
    dar_face = coords(darpin, residues=set(para), sidechain_only=True)
    face_idx = cKDTree(dar_all).query(dar_face, k=1)[1].tolist()

    target = coords(complex_model, chains={"B", "C"})
    peptide = coords(complex_model, chains={"D"})
    b2m = coords(complex_model, chains={"C"})
    binder, binder_r = coords(complex_model, chains={"A"}, radii=True)

    axis = peptide.mean(0) - b2m.mean(0)
    axis /= np.linalg.norm(axis)
    origin = peptide.mean(0)

    target_tree = cKDTree(np.vstack([target, peptide]))
    peptide_tree = cKDTree(peptide)

    for label, lig, rad, idx in (
            ("DARPin E3_5", dar_all, dar_r, face_idx),
            ("de novo binder (control)", binder, binder_r, None)):
        if idx is None:   # control: its face is the surface it actually uses
            hits = cKDTree(lig).query_ball_point(peptide, 5.0)
            idx = sorted({j for h in hits for j in h})
            print(f"  control interface atoms: {len(idx)}")
        n, pts = scan(lig, idx, target_tree, peptide_tree, axis, origin)
        bsa = buried_peptide_surface(complex_model, pts, rad)
        print(f"{label:28s} best placement: {n:4d} peptide contacts, "
              f"{bsa:6.1f} A^2 peptide BSA", flush=True)

    print(f"\nreference: the submitted pose of this binder buries "
          f"{buried_peptide_surface(complex_model, binder, binder_r):.1f} A^2")


if __name__ == "__main__":
    main()
