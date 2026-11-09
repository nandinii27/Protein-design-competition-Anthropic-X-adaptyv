# Protein Design — Anthropic × Adaptyv 2026

Design and scoring code for all five challenges of the 2026 protein design competition on Proteinbase. `binderkit` is the shared library; each challenge directory holds its own run code, methods document and submission.

**Every one of the five asks for a binder, and in none of them does affinity decide the ranking.** Four are conditional-binding problems — bind here, release there, or bind without activating — and the fifth asks a binder to tell two peptides apart that differ by one buried side chain. No generation stack produces a signal for any of those. Each challenge directory therefore begins by identifying what the scoring actually rewards, and then builds the missing signal, because the part everyone's pipeline handles is the part that was never in question.

## The five challenges

| directory | target | what the scoring actually asks |
|---|---|---|
| `01-egfr` | EGFR domain III | Bind at pH 6.5, release at pH 7.4, cross-react with the mouse orthologue, and land on the functional epitope — in that priority order. No structure predictor takes pH as an input, so folding the complex at either pH returns the same answer. |
| `02-tnfa` | TNF-α soluble trimer (1TNF) | The same switch running the other way: bind at pH 7.4 with no detectable binding at pH 6.0, plus mouse cross-reactivity. Histidine is the only side chain whose pKa falls inside that window, which is why grafting it is the standard move and why the standard move is not enough. |
| `03-glp1r` | GLP-1R | Binding is the easy half. The second objective is **activation**, and the extracellular-domain groove the challenge page recommends is precisely the architecture that yields antagonists — the same outcome the 2025 Baker-lab GPCR miniprotein work reached at GLP1R, GIPR, GCGR and CGRPR. |
| `04-tfr1` | TfR1 apical domain | Bind at pH 7.4, release at pH 6.0, cross-react with mouse TfR1, and stay outside the transferrin footprint. |
| `05-pmhc` | HLA-A\*11:01 / KRAS G12D | Bind the neoantigen decamer `VVVGADGVGK`, discriminate it from wild type, and depend on the peptide rather than the groove. The mutated aspartate presents 8.8 Å² and the reference TCR contacts none of it. |

Each directory documents its own reasoning in `methods.md` or `README.md`. This file covers what is shared.

## Where the missing signal comes from

**pH, in `01` and `02`.** The signal is injected twice. Upstream, the epitope is restricted to the intersection of conserved, acidic and solvent-exposed surface — protonated histidines need carboxylate partners, and conservation gates the cross-species objective. That intersection is small and determines everything downstream, so it is computed first. Downstream, predicted complexes are scored with a thermodynamic linkage function over predicted pKa shifts, gated by a geometric hard filter on histidine-to-carboxylate contacts. **The scorer was validated against known systems before being trusted on designs, where it recovered the FcRn switch histidines blind.** It remains unvalidated on designed interfaces.

**Target-side chemistry, in `04`.** Rather than requiring the binder to carry the switch, the selected epitope is centred on a conserved, solvent-exposed histidine on TfR1 that titrates across the assay window by itself. The same selection independently recovered the apical species-specificity loop — `R208N`, the `V210` deletion, `Y211D`, `L212P` — from conservation and solvent exposure alone, which is the patch the viral-entry literature identifies.

**Architecture, in `03`.** Hotspots span both receptor domains rather than the extracellular groove alone, because the groove supplies the affinity and the transmembrane insertion supplies the activation. Selecting on the recommended surface optimises the objective that was never in doubt.

**Conformation, in `05`.** Discrimination cannot come from contact, so it is sought in dynamics: replica MD on both peptide variants puts position 6 at 2.66× the mobility in wild type, with anchor residues unchanged. Filtering is built on peptide-specific buried surface, target-free monomer refold, a ProteinMPNN readout test with apo controls, and a three-target counter-screen.

## `binderkit`

```
src/binderkit/
  epitope.py      conservation, SASA, acidic intersection, patch ranking
  phscore.py      pKa prediction free vs bound, linkage score, His geometry
  liability.py    cysteines, hydrophobic patch, motifs, length, novelty
  select.py       Pareto front over objectives, diversity pick
  submission.py   CSV builder and validator
```

PROPKA3 is pulled in as a dependency and runs on CPU in seconds per structure. Backbone generation, sequence design and structure prediction are external and not wrapped here; this package handles epitope definition and everything after folding.

```
pip install -e .
```

Challenges carry their own extras where they need them — `03-glp1r` has a `requirements.txt`, and `05-pmhc` installs `boltzgen`, `chai_lab`, `openmm`, `tblite` and `prody` per stage inside its notebook, since the heavy ones are only useful on a GPU session.

## The general shape of a run

```
read the scoring rubric          identify the objective that decides the
        |                        ranking, which is never affinity alone
epitope or target definition     restrict the surface before generating
        |                        anything; this determines everything after
backbone generation              external, aimed at the selected hotspots
        |
sequence design                  external
        |
structure prediction             external
        |
cheap geometric filter           discard the bulk before spending GPU on
        |                        co-folding
objective-specific scoring       the filter built for this challenge's real
        |                        objective, with a control
liability and novelty            hard filters, including sequence identity
        |                        against the natural proteins in the target
Pareto front                     no weighted sum; the objectives are not
        |                        commensurable
diversity pick                   max-min sequence distance over the front
        |
submission                       validated CSV, best row first
```

**Every filter carries a control, and the controls have earned their place.** The pH scorer was calibrated on FcRn and barnase–barstar before being applied to a design. The ProteinMPNN readout test in `05` is run twice per design, with and without the binder, because the groove already prefers the mutant residue and only the shift attributable to the binder means anything. A rigid-body docking scan in the same challenge was reported as inconclusive because its known-positive control recovered only 42% of a pose it should have found. Two results in this repository reversed once a control was added.

## Conventions

- Predicted complexes: binder as chain `B`, target as chain `A`, one PDB per design, named by design identifier. `05-pmhc` departs from this because its target is a trimer — see its README.
- Free energies in kcal/mol. Negative `ph_selectivity_kcal` means binding is more favourable at the lower pH, which is the wanted direction for `01` and `04` and the wrong one for `02`.
- Thresholds in function signatures are defaults, not published constants. Recalibrate against whatever reference data is available and record what was used in the methods document.
- Prior competition data is used for threshold calibration only, never upstream of generation. The de novo rule forbids deriving a design from an existing binder, and earlier rounds optimised affinity alone, so their successful epitopes are actively wrong for objectives that are not affinity.
- **Submitted sequences are checked by local alignment against every natural protein in the target assembly.** In `05-pmhc`, seven of twenty-five designs came back at 85–93% identity to β2-microglobulin; they pass every structural filter, because they are a real protein, and only an explicit identity check removes them.
