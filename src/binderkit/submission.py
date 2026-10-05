"""Submission CSV writer and validator.

Required columns are name, sequence and molecule_class. Rows are ordered
best-first by the submitter's own ranking. Extra metric columns are kept and
written after the required ones: the selection stage reads them, so more
legible information is strictly better than less.

Multi-chain formats carry their chains in one sequence field separated by a
colon, heavy chain first.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ["name", "sequence", "molecule_class"]

MOLECULE_CLASSES = {"protein", "nanobody", "scfv", "fab_kappa", "fab_lambda"}

MULTI_CHAIN_CLASSES = {"scfv", "fab_kappa", "fab_lambda"}

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

MAX_ROWS = 20

# Columns worth carrying through when present, in the order they are written.
METRIC_COLUMNS = [
    "ph_selectivity_kcal",
    "ddg_ph_low",
    "ddg_ph_high",
    "max_his_pka_shift",
    "n_his_carboxylate_contacts",
    "iptm_human",
    "pae_interaction_human",
    "iptm_mouse",
    "pae_interaction_mouse",
    "self_consistency_rmsd",
    "plddt_binder",
    "length",
    "n_cys",
    "net_charge_ph7",
    "surface_hydrophobic_patch",
    "max_identity_to_reference",
    "epitope_patch",
    "pareto_rank",
    "min_distance_to_selected",
]


@dataclass
class ValidationResult:
    ok: bool
    errors: list
    warnings: list

    def report(self) -> str:
        lines = []
        for error in self.errors:
            lines.append(f"ERROR {error}")
        for warning in self.warnings:
            lines.append(f"WARN  {warning}")
        if not lines:
            lines.append("OK")
        return "\n".join(lines)


def _chains(sequence: str) -> list[str]:
    return [part.strip() for part in sequence.split(":")]


def validate(df: pd.DataFrame, max_rows: int = MAX_ROWS) -> ValidationResult:
    """Structural validation of a submission table.

    Errors block submission. Warnings are judgement calls left to the submitter.
    """
    errors: list[str] = []
    warnings: list[str] = []

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        return ValidationResult(False, [f"missing required columns: {missing}"], [])

    if len(df) == 0:
        errors.append("submission is empty")
    if len(df) > max_rows:
        errors.append(f"{len(df)} rows exceeds the maximum of {max_rows}")

    names = df["name"].astype(str)
    if names.duplicated().any():
        errors.append(f"duplicate names: {sorted(names[names.duplicated()].unique())}")
    if names.str.strip().eq("").any():
        errors.append("blank name field")

    sequences = df["sequence"].astype(str)
    if sequences.duplicated().any():
        errors.append(f"duplicate sequences at rows {list(sequences[sequences.duplicated()].index)}")

    for index, row in df.iterrows():
        molecule_class = str(row["molecule_class"]).strip().lower()
        sequence = str(row["sequence"]).strip().upper()

        if molecule_class not in MOLECULE_CLASSES:
            errors.append(f"row {index}: unknown molecule_class {molecule_class!r}")

        chains = _chains(sequence)
        if molecule_class in MULTI_CHAIN_CLASSES and len(chains) != 2:
            errors.append(f"row {index}: {molecule_class} needs two colon separated chains")
        if molecule_class not in MULTI_CHAIN_CLASSES and len(chains) != 1:
            errors.append(f"row {index}: {molecule_class} must be a single chain")

        for chain_index, chain in enumerate(chains):
            invalid = sorted(set(chain) - VALID_AA)
            if invalid:
                errors.append(f"row {index} chain {chain_index}: non-standard residues {invalid}")
            if not chain:
                errors.append(f"row {index} chain {chain_index}: empty sequence")
            elif len(chain) < 20:
                warnings.append(f"row {index} chain {chain_index}: unusually short, {len(chain)} aa")

        if "H" in sequence and sequence.count("H") < 2:
            warnings.append(f"row {index}: fewer than two histidines, weak pH mechanism")
        if "H" not in sequence:
            warnings.append(f"row {index}: no histidine, no pH switch mechanism")

    return ValidationResult(not errors, errors, warnings)


def build(
    df: pd.DataFrame,
    name_column: str = "name",
    sequence_column: str = "sequence",
    molecule_class: str = "protein",
    name_prefix: str | None = None,
    extra_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Assemble the submission table in required-column-first order.

    Row order is preserved: the caller has already ranked best-first.
    """
    out = pd.DataFrame()
    if name_prefix:
        out["name"] = [f"{name_prefix}-{i + 1:03d}" for i in range(len(df))]
    else:
        out["name"] = df[name_column].astype(str).values

    out["sequence"] = df[sequence_column].astype(str).str.upper().str.strip().values

    if "molecule_class" in df.columns:
        out["molecule_class"] = df["molecule_class"].astype(str).str.lower().values
    else:
        out["molecule_class"] = molecule_class

    carried = extra_columns if extra_columns is not None else METRIC_COLUMNS
    for column in carried:
        if column in df.columns:
            out[column] = df[column].values
    return out


def write(df: pd.DataFrame, path: str | Path, validate_first: bool = True) -> ValidationResult:
    """Write the CSV, refusing on validation errors unless checking is disabled."""
    result = validate(df) if validate_first else ValidationResult(True, [], [])
    if validate_first and not result.ok:
        return result
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return result


def read(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"name": str, "sequence": str, "molecule_class": str})
