"""
Geometric triage of BoltzGen designs - CPU only, no GPU required.

Answers the question that decides whether a design can possibly be specific:
does its interface land on the peptide, or on the MHC?

The alpha1/alpha2 helices are identical across every peptide HLA-A*11:01
presents, and b2-microglobulin is not part of the discriminating surface at
all. A binder that reads either will bind the wild-type complex and the empty
MHC just as well as the target, failing two of the three scored objectives.
Those designs are not worth GPU time.

Chains are identified by content, not by letter, so this works on runs with
or without b2m in the target:
    peptide   10 residues
    alpha1-a2 182 residues
    b2m       99 residues matching the human b2m sequence
    binder    whatever remains

Columns
    pep_bsa_frac   buried target surface on the peptide, as a fraction of all
                   buried target surface. The headline number; aim >= 0.30.
    contacts_4_8   binder heavy atoms within 5 A of peptide positions 4-8.
    hotspot_cover  how many of positions 4-8 are contacted at all.
    b2m_identity   percent identity of the binder to human b2-microglobulin.
                   High values mean the model regenerated the natural partner
                   rather than designing anything.

Usage:
    python geometric_filter.py <design_dir> <output_csv>
"""

from __future__ import annotations

import copy
import csv
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import three_to_index, index_to_one

B2M = ("IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLY"
       "YTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM")
A1A2_START = "GSHSMRY"
PEPTIDE_LEN = 10
HOTSPOTS = (4, 5, 6, 7, 8)
CUTOFF = 5.0

sr = ShrakeRupley()
parser = MMCIFParser(QUIET=True)


def one(residue) -> str:
    try:
        return index_to_one(three_to_index(residue.resname))
    except Exception:
        return "X"


def identity(query: str, reference: str) -> float:
    q = query.lstrip("M")
    n = min(len(q), len(reference))
    if n == 0:
        return 0.0
    return 100.0 * sum(a == b for a, b in zip(q[:n], reference[:n])) / n


def classify(model) -> dict | None:
    """Map chain ids to roles by content."""
    chains = {}
    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        if residues:
            chains[chain.id] = "".join(one(r) for r in residues)

    roles, remaining = {}, dict(chains)
    for cid, seq in list(remaining.items()):
        if len(seq) == PEPTIDE_LEN:
            roles["peptide"] = cid
            del remaining[cid]
            break
    for cid, seq in list(remaining.items()):
        if seq.startswith(A1A2_START):
            roles["mhc"] = cid
            del remaining[cid]
            break
    for cid, seq in list(remaining.items()):
        if len(seq) == len(B2M) and identity(seq, B2M) > 95:
            roles["b2m"] = cid
            del remaining[cid]
            break
    if len(remaining) != 1 or "peptide" not in roles or "mhc" not in roles:
        return None
    roles["binder"] = next(iter(remaining))
    roles["_seqs"] = chains
    return roles


def sasa_by_residue(model, keep: set[str]) -> dict:
    sub = copy.deepcopy(model)
    for chain in list(sub):
        if chain.id not in keep:
            sub.detach_child(chain.id)
    sr.compute(sub, level="R")
    return {(c.id, r.id[1]): r.sasa for c in sub for r in c if r.id[0] == " "}


def analyse(path: Path) -> dict | None:
    try:
        model = parser.get_structure(path.stem, str(path))[0]
    except Exception as error:
        print(f"  {path.name}: parse failed - {error}")
        return None

    roles = classify(model)
    if roles is None:
        print(f"  {path.name}: could not classify chains")
        return None

    binder_id, mhc_id, pep_id = roles["binder"], roles["mhc"], roles["peptide"]
    b2m_id = roles.get("b2m")
    target_ids = {mhc_id, pep_id} | ({b2m_id} if b2m_id else set())
    sequence = roles["_seqs"][binder_id]

    free = sasa_by_residue(model, target_ids)
    bound = sasa_by_residue(model, target_ids | {binder_id})
    buried = lambda cid: sum(
        max(0.0, free[k] - bound[k]) for k in free if k[0] == cid
    )
    pep_bsa, mhc_bsa = buried(pep_id), buried(mhc_id)
    b2m_bsa = buried(b2m_id) if b2m_id else 0.0
    total = pep_bsa + mhc_bsa + b2m_bsa

    peptide = [r for r in model[pep_id] if r.id[0] == " "]
    binder_atoms = np.array(
        [a.coord for r in model[binder_id] if r.id[0] == " "
         for a in r if a.element != "H"]
    )
    contacts, covered = 0, set()
    for index, residue in enumerate(peptide, start=1):
        if index not in HOTSPOTS:
            continue
        coords = np.array([a.coord for a in residue if a.element != "H"])
        n = int((np.linalg.norm(
            binder_atoms[:, None, :] - coords[None, :, :], axis=-1) < CUTOFF).sum())
        contacts += n
        if n:
            covered.add(index)

    return {
        "design": path.stem,
        "pep_bsa_frac": round(pep_bsa / total, 3) if total else 0.0,
        "contacts_4_8": contacts,
        "hotspot_cover": len(covered),
        "pep_bsa": round(pep_bsa, 1),
        "mhc_bsa": round(mhc_bsa, 1),
        "b2m_bsa": round(b2m_bsa, 1),
        "length": len(sequence),
        "n_cys": sequence.count("C"),
        "b2m_identity": round(identity(sequence, B2M), 1),
        "hotspots_hit": ",".join(str(i) for i in sorted(covered)),
        "binder_chain": binder_id,
        "sequence": sequence,
    }


def main() -> None:
    if len(sys.argv) < 3:
        sys.exit("usage: geometric_filter.py <design_dir> <output_csv>")
    design_dir, out_csv = Path(sys.argv[1]), Path(sys.argv[2])

    files = sorted(design_dir.glob("*.cif"))
    if not files:
        sys.exit(f"no .cif files in {design_dir}")
    print(f"{len(files)} designs in {design_dir}")

    probe = parser.get_structure("probe", str(files[0]))[0]
    roles = classify(probe)
    if roles is None:
        inv = {c.id: sum(1 for r in c if r.id[0] == " ") for c in probe}
        sys.exit(f"cannot classify chains in {files[0].name}: {inv}")
    print(f"roles: binder={roles['binder']} mhc={roles['mhc']} "
          f"peptide={roles['peptide']} b2m={roles.get('b2m', '-')}\n")

    rows = []
    for n, path in enumerate(files, start=1):
        row = analyse(path)
        if row:
            rows.append(row)
        if n % 25 == 0:
            print(f"  {n}/{len(files)}")

    if not rows:
        sys.exit("nothing analysed")

    rows.sort(key=lambda r: (r["pep_bsa_frac"], r["contacts_4_8"]), reverse=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    frac = np.array([r["pep_bsa_frac"] for r in rows])
    print(f"\nanalysed {len(rows)} designs -> {out_csv}")
    print(f"peptide BSA fraction: median {np.median(frac):.3f}  max {frac.max():.3f}")
    for cut in (0.10, 0.20, 0.30, 0.40):
        print(f"  >= {cut:.2f}: {int((frac >= cut).sum())}")
    print(f"  touches any hotspot : {sum(1 for r in rows if r['hotspot_cover'] > 0)}")
    print(f"  all 5 hotspots      : {sum(1 for r in rows if r['hotspot_cover'] == 5)}")
    print(f"  cysteine-free       : {sum(1 for r in rows if r['n_cys'] == 0)}")
    print(f"  b2m-like (>80% id)  : {sum(1 for r in rows if r['b2m_identity'] > 80)}")

    print(f"\n{'design':<20} {'frac':>6} {'pep':>7} {'mhc':>8} {'cov':>4} "
          f"{'len':>4} {'cys':>4} {'b2m%':>6}")
    for r in rows[:20]:
        print(f"{r['design']:<20} {r['pep_bsa_frac']:>6.3f} {r['pep_bsa']:>7.1f} "
              f"{r['mhc_bsa']:>8.1f} {r['hotspot_cover']:>4} {r['length']:>4} "
              f"{r['n_cys']:>4} {r['b2m_identity']:>6.1f}")


if __name__ == "__main__":
    main()
