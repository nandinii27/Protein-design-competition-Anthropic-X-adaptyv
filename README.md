# Protein Design

Design and scoring code for the 2026 protein design competition challenges.

`binderkit` is the shared library; each challenge directory holds its own run
notebook, methods document and submission.

## The problem this code exists to solve

Challenge 01 asks for a de novo binder that engages its target at pH 6.5 and
releases at pH 7.4, cross-reacts with the mouse orthologue, and binds the
functional epitope — in that priority order.

No available structure predictor takes pH as an input. Folding a complex at 6.5
and at 7.4 returns the identical structure, so a conventional generation stack
produces a signal for affinity, a signal for cross-reactivity, and nothing at
all on the highest-weighted objective.

This package injects the missing signal in two places:

- **Upstream**, by restricting the epitope to the intersection of conserved,
  acidic and solvent-exposed surface. Protonated histidines need carboxylate
  partners; conservation gates the cross-species objective. The intersection is
  small and it determines everything downstream, so it is computed first.
- **Downstream**, by scoring predicted complexes with a thermodynamic linkage
  function over predicted pKa shifts, gated by a geometric hard filter on
  histidine to carboxylate contacts.

The linkage score is unvalidated on designed interfaces. 

## Layout

    src/binderkit/
      epitope.py      conservation, SASA, acidic intersection, patch ranking
      phscore.py      pKa prediction free vs bound, linkage score, His geometry
      liability.py    cysteines, hydrophobic patch, motifs, length, novelty
      select.py       Pareto front over objectives, diversity pick
      submission.py   CSV builder and validator

    challenges/01-egfr/
      run.ipynb       exploration and run log
      methods.md      submitted methods description
      submission.csv  final table, best row first
      outputs/        intermediate tables and predicted structures

## Install

    pip install -e .

PROPKA3 is pulled in as a dependency and runs on CPU in seconds per structure.
Backbone generation, sequence design and structure prediction are external and
not wrapped here; this package handles epitope definition and everything after
folding.

## Pipeline

    epitope definition        conserved AND acidic AND exposed, on the
                              functional surface
            |
    backbone generation       external, aimed at the selected hotspots
            |
    sequence design           external, histidine fixed facing carboxylates
            |
    structure prediction      external, complex with the human target
            |
    cross-species refold      survivors re-folded against the mouse sequence
            |
    pH scoring                pKa shift free vs bound, linkage sum, geometry
                              filter
            |
    liability and novelty     hard filters
            |
    Pareto front              no weighted sum, objectives are not commensurable
            |
    diversity pick            max-min sequence distance over the front
            |
    submission                validated CSV, best row first

Prior competition data is used for threshold calibration only, never upstream
of generation: the de novo rule forbids deriving designs from an existing
binder, and the prior round optimised affinity alone, so its successful
epitopes are actively wrong for the cross-species objective.

## Conventions

- Predicted complexes: binder as chain `B`, target as chain `A`, one PDB per
  design, named by design identifier.
- Free energies in kcal/mol. Negative `ph_selectivity_kcal` means binding is
  more favourable at the lower pH, which is the wanted direction.
- Thresholds in function signatures are defaults, not published constants.
  Recalibrate them against whatever reference data is available and record what
  was used in the methods document.
