# Challenge 04 — TfR1 conditional binder

Anthropic × Adaptyv protein design competition, Proteinbase.
Deadline: Sunday 1 November 2026, 23:59 AoE. **Submitted.**

Design a binder that engages human TfR1 at pH 7.4, releases at pH 6.0, and
cross-reacts with mouse TfR1. Full rationale in `methods.md`.

## Outcome

10 designs submitted, generated through NVIDIA's hosted RFD4-Proteina endpoint.
Unfiltered, no pH engineering applied. The intended route — Proteina-Complexa on
a local GPU — did not produce designs in the time available; the install fixes
are in the notebook and documented in `methods.md`.

The work worth reading here is the epitope analysis, not the designs.

## Epitope

Both candidate patches are conserved, solvent-exposed, spatially contiguous
regions of the apical domain, outside the transferrin footprint. Numbering is
UniProt / PDB numbering of mature human TfR1.

- **EP1, conservation-first:** `K240, K241, E244, I262, E266, S273`
- **EP2, pH handle (used):** `T263, F264, E266, A269, H300, L303` — centred on a
  conserved, solvent-exposed histidine that titrates between pH 7.4 and 6.0 on
  the *target* side, so the binder does not have to carry the switch itself.

The selection independently recovered the apical species-specificity loop
(`R208N`, the `V210` deletion, `Y211D`, `L212P`) from conservation and solvent
exposure alone — the same patch the viral-entry literature identifies.

## Layout

```
04-tfr1/
  methods.md                        epitope rationale, install notes, limitations
  proteina_complexa_colab.ipynb     Colab notebook with the full install fix
  targets/
    human_TfR1_apical.pdb           apical domain 189-383, chain A
    mouse_TfR1_apical.pdb           threaded mouse model, same frame
    human_TfR1_ectodomain_chainA.pdb  full monomer 121-760
    tfr1_residue_table.csv          per-residue SASA, conservation, domain
    epitope_summary.json            selection output
  raw/                              inputs for reproducing targets/
  outputs/                          design outputs
```

`../../src/tfr1_epitope.py` regenerates everything in `targets/` from `raw/`.

## Reproducing the target preparation

```powershell
python -I src\tfr1_epitope.py challenges\04-tfr1\raw challenges\04-tfr1\targets
```

Inputs: PDB `11QC` (RCSB), UniProt `P02786` (human), `Q62351` (mouse, fetched via
the EBI Proteins API — the UniProt REST endpoint was returning 500 for this
accession).

The script validates every PDB it writes. That check exists because an earlier
hand-written writer produced a column-misaligned file that parsed as plausible
garbage and killed a design run minutes in.

## Running Proteina-Complexa

`proteina_complexa_colab.ipynb` — needs an A100 or H100, mirrors checkpoints and
outputs to Google Drive so a runtime reset costs only the environment rebuild.

Five install defects are patched in it, all stemming from the repo assuming its
own Docker image. The one most likely to catch others: `atomworks` installs
unversioned behind a `|| echo "Warning:"`, silently pulls 3.x, and 3.x dropped
`AtomSelectionStack`. **Pin 2.2.1.** See `methods.md` for the full list.

Outstanding: on a 16 GB card, JAX preallocates most of the VRAM before PyTorch
gets any. `XLA_PYTHON_CLIENT_PREALLOCATE=false` is the untested first thing to
try.

## Not in this folder

`submission.csv` and the RFD4-Proteina generation script — add from the session
where the designs were produced.
