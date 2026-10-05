"""pH-selectivity scoring by thermodynamic linkage on predicted pKa shifts.

No current structure predictor takes pH as an input, so folding a complex at
pH 6.5 and at 7.4 returns the identical structure and the highest-weighted
objective gets no signal from the forward model. This module supplies that
signal from an empirical pKa predictor instead.

For a single titratable site, the pH dependence of the binding free energy is
set by how much the site's pKa shifts on complex formation:

    ddG(pH) = -RT * ln[ (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH)) ]

Summed over interface ionisable sites. The selectivity quantity is the
difference between the two assay conditions:

    pH_selectivity = ddG(pH_low) - ddG(pH_high)

A histidine whose pKa is raised on binding stabilises its protonated form in
the complex, so binding is favoured at the lower pH. That is the wanted switch.
A pKa that drops on binding gives the inverse behaviour and the design is
rejected rather than merely down-ranked.

Caveat to carry into the methods document: the pKa predictor is an empirical
model parameterised largely on natural proteins. Its accuracy on de novo
designed interfaces is unknown and this score is unvalidated.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

RT_KCAL = 0.001987 * 298.15  # kcal/mol at 25 C

IONISABLE = {"HIS", "ASP", "GLU", "LYS", "ARG", "CYS", "TYR"}

HIS_NITROGENS = ("ND1", "NE2")
CARBOXYLATE_OXYGENS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}

# propka writes a fixed-width summary block; residue, number, chain, pKa.
PKA_LINE = re.compile(r"^\s*([A-Z]{2,3})\s+(\d+)\s+([A-Za-z0-9])\s+(-?\d+\.\d+)")


@dataclass(frozen=True)
class Site:
    resname: str
    resnum: int
    chain: str

    def __str__(self) -> str:
        return f"{self.resname}{self.resnum}_{self.chain}"


class _ChainSelect(Select):
    def __init__(self, keep: set[str]):
        self.keep = keep

    def accept_chain(self, chain):
        return chain.get_id() in self.keep

    def accept_residue(self, residue):
        return residue.get_id()[0] == " "


def write_chains(structure_path: str | Path, chains: set[str], out_path: str | Path) -> Path:
    """Extract a chain subset into its own PDB, hetero records dropped."""
    structure = PDBParser(QUIET=True).get_structure("s", str(structure_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), _ChainSelect(chains))
    return Path(out_path)


def _find_pka(stem: str, *directories: Path) -> Path | None:
    """Locate a .pka file by stem across the directories it may have landed in.

    The predictor resolves its output path differently across versions and entry
    points: some write next to the input, others into the process working
    directory. Searching both is cheaper than depending on which.
    """
    seen: set[Path] = set()
    for directory in directories:
        directory = Path(directory).resolve()
        if directory in seen or not directory.is_dir():
            continue
        seen.add(directory)
        candidate = directory / f"{stem}.pka"
        if candidate.exists():
            return candidate
    return None


def run_propka(pdb_path: str | Path, workdir: str | Path | None = None) -> dict[Site, float]:
    """Run the pKa predictor on one structure and parse the summary block.

    The structure is staged into its own directory and the predictor is invoked
    with that directory as the process working directory, so concurrent designs
    cannot overwrite each other's output.
    """
    pdb_path = Path(pdb_path).resolve()
    workdir = Path(workdir).resolve() if workdir else pdb_path.parent
    workdir.mkdir(parents=True, exist_ok=True)
    staged = workdir / pdb_path.name
    if staged != pdb_path:
        shutil.copy(pdb_path, staged)

    origin = Path.cwd()
    errors: list[str] = []

    try:
        os.chdir(workdir)
        try:
            from propka.run import single

            single(staged.name, optargs=["--quiet"], write_pka=True)
        except Exception as exc:
            errors.append(f"python api: {exc}")
            binary = shutil.which("propka3") or shutil.which("propka")
            if binary is None:
                errors.append("no propka executable on PATH")
            else:
                completed = subprocess.run(
                    [binary, "--quiet", staged.name],
                    capture_output=True,
                    text=True,
                )
                if completed.returncode != 0:
                    errors.append(f"cli: {completed.stderr.strip()[:200]}")
    finally:
        os.chdir(origin)

    pka_file = _find_pka(staged.stem, workdir, origin, pdb_path.parent)
    if pka_file is None:
        detail = "; ".join(errors) if errors else "no error reported"
        raise FileNotFoundError(f"no pKa output for {staged.name} ({detail})")

    if pka_file.parent != workdir:
        shutil.move(str(pka_file), str(workdir / pka_file.name))
        pka_file = workdir / pka_file.name

    return parse_pka(pka_file)


def parse_pka(pka_path: str | Path) -> dict[Site, float]:
    """Parse the SUMMARY section of a .pka file into site -> predicted pKa."""
    values: dict[Site, float] = {}
    in_summary = False
    for line in Path(pka_path).read_text().splitlines():
        if line.strip().startswith("SUMMARY OF THIS PREDICTION"):
            in_summary = True
            continue
        if in_summary:
            if line.strip().startswith("-") or line.strip().startswith("Free energy"):
                break
            match = PKA_LINE.match(line)
            if match:
                resname, resnum, chain, pka = match.groups()
                if resname in IONISABLE:
                    values[Site(resname, int(resnum), chain)] = float(pka)
    return values


def linkage_ddg(pka_free: float, pka_bound: float, ph: float) -> float:
    """Protonation contribution to binding free energy at one pH, kcal/mol.

    Negative values favour binding.
    """
    numerator = 1.0 + 10.0 ** (pka_bound - ph)
    denominator = 1.0 + 10.0 ** (pka_free - ph)
    return -RT_KCAL * math.log(numerator / denominator)


def interface_sites(
    complex_path: str | Path,
    binder_chain: str,
    target_chain: str,
    cutoff: float = 6.0,
) -> set[Site]:
    """Ionisable binder sites with any heavy atom within cutoff of the target."""
    structure = PDBParser(QUIET=True).get_structure("c", str(complex_path))[0]
    target_atoms = np.array(
        [a.get_coord() for a in structure[target_chain].get_atoms() if a.element != "H"]
    )
    sites: set[Site] = set()
    for residue in structure[binder_chain]:
        if residue.get_id()[0] != " " or residue.get_resname() not in IONISABLE:
            continue
        coords = np.array([a.get_coord() for a in residue if a.element != "H"])
        distances = np.linalg.norm(coords[:, None, :] - target_atoms[None, :, :], axis=-1)
        if distances.min() <= cutoff:
            sites.add(Site(residue.get_resname(), residue.get_id()[1], binder_chain))
    return sites


def his_carboxylate_contacts(
    complex_path: str | Path,
    binder_chain: str,
    target_chain: str,
    cutoff: float = 4.0,
) -> list[dict]:
    """Binder histidine nitrogen to target carboxylate oxygen pairs within cutoff.

    This is the hard geometric filter. A design with no such pair has no
    structural mechanism for pH switching regardless of what the linkage score
    says, and is discarded before ranking.
    """
    model = PDBParser(QUIET=True).get_structure("c", str(complex_path))[0]
    acidic = [r for r in model[target_chain] if r.get_resname() in CARBOXYLATE_OXYGENS]
    contacts = []
    for residue in model[binder_chain]:
        if residue.get_resname() != "HIS":
            continue
        for nitrogen in HIS_NITROGENS:
            if nitrogen not in residue:
                continue
            n_coord = residue[nitrogen].get_coord()
            for partner in acidic:
                for oxygen in CARBOXYLATE_OXYGENS[partner.get_resname()]:
                    if oxygen not in partner:
                        continue
                    distance = float(np.linalg.norm(n_coord - partner[oxygen].get_coord()))
                    if distance <= cutoff:
                        contacts.append(
                            {
                                "his_resnum": residue.get_id()[1],
                                "his_atom": nitrogen,
                                "target_residue": f"{partner.get_resname()}{partner.get_id()[1]}",
                                "target_atom": oxygen,
                                "distance": round(distance, 2),
                            }
                        )
    return contacts


def score_design(
    complex_path: str | Path,
    binder_chain: str,
    target_chain: str,
    ph_low: float = 6.5,
    ph_high: float = 7.4,
    interface_cutoff: float = 6.0,
    contact_cutoff: float = 4.0,
    min_his_contacts: int = 2,
    workdir: str | Path | None = None,
) -> dict:
    """Full pH assessment of one predicted complex.

    Runs the pKa predictor twice, once on the isolated binder and once on the
    complex, restricts to interface sites, and sums the linkage contributions.
    Returns the selectivity score, the per-site shifts, the geometric contacts
    and a pass flag combining the hard filter with the sign of the score.
    """
    complex_path = Path(complex_path)
    temp_root = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="phscore_"))
    temp_root.mkdir(parents=True, exist_ok=True)

    free_pdb = write_chains(complex_path, {binder_chain}, temp_root / f"{complex_path.stem}_free.pdb")
    bound_pdb = write_chains(
        complex_path, {binder_chain, target_chain}, temp_root / f"{complex_path.stem}_bound.pdb"
    )

    pka_free = run_propka(free_pdb, temp_root)
    pka_bound = run_propka(bound_pdb, temp_root)

    sites = interface_sites(complex_path, binder_chain, target_chain, interface_cutoff)

    per_site = []
    ddg_low = ddg_high = 0.0
    for site in sorted(sites, key=lambda s: s.resnum):
        free = pka_free.get(site)
        bound = pka_bound.get(site)
        if free is None or bound is None:
            continue
        low = linkage_ddg(free, bound, ph_low)
        high = linkage_ddg(free, bound, ph_high)
        ddg_low += low
        ddg_high += high
        per_site.append(
            {
                "site": str(site),
                "pka_free": round(free, 2),
                "pka_bound": round(bound, 2),
                "pka_shift": round(bound - free, 2),
                "ddg_low": round(low, 3),
                "ddg_high": round(high, 3),
            }
        )

    contacts = his_carboxylate_contacts(complex_path, binder_chain, target_chain, contact_cutoff)
    distinct_his = {c["his_resnum"] for c in contacts}

    selectivity = ddg_low - ddg_high
    his_shifts = [s["pka_shift"] for s in per_site if s["site"].startswith("HIS")]

    return {
        "design": complex_path.stem,
        "ph_selectivity_kcal": round(selectivity, 3),
        "ddg_ph_low": round(ddg_low, 3),
        "ddg_ph_high": round(ddg_high, 3),
        "n_interface_ionisable": len(per_site),
        "n_his_carboxylate_contacts": len(contacts),
        "n_distinct_his_in_contact": len(distinct_his),
        "max_his_pka_shift": round(max(his_shifts), 2) if his_shifts else None,
        "geometry_pass": len(distinct_his) >= min_his_contacts,
        "sign_pass": selectivity < 0.0,
        "ph_pass": len(distinct_his) >= min_his_contacts and selectivity < 0.0,
        "per_site": per_site,
        "contacts": contacts,
    }