# Experimental — burial-aware multi-histidine grafting

Not run yet due to time constraints. 

## Why

The submitted pipeline grafts histidines selected on reach and orientation
only, one or two at a time. Both choices cap the achievable pH switch, and the
cap is quantifiable.

For a single titratable site the linkage relation gives

    ddG(pH) = -RT ln[(1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH))]

Take a generous case, pKa 6.0 free shifting to 8.0 on binding:

    pH 6.5:  ddG = -0.592 * ln(24.8) = -1.90 kcal/mol
    pH 7.4:  ddG = -0.592 * ln(4.79) = -0.93 kcal/mol
    selectivity = -0.97 kcal/mol, about five-fold

So one histidine with an unusually large shift buys roughly five-fold. Ten- to
hundred-fold needs 1.4-2.7 kcal/mol, which means two or three histidines acting
together. The best design in the submission reached -1.21 kcal/mol.

The second limit is burial. A pKa shift comes from the environment changing
around the residue. A histidine that stays solvent-exposed in both states sees
almost no change, so its pKa barely moves. Large shifts need the side chain to
become buried on binding, where desolvation and a nearby carboxylate stabilise
the protonated form. The submitted filter tested distance and orientation and
said nothing about burial, which is the most likely reason the observed shifts
were small.

## What `graft_v2.py` does differently

- **Requires burial.** A position must be exposed in the free binder
  (>20 A^2) and lose at least 15 A^2 of accessibility on complex formation.
- **Allows up to three simultaneous substitutions**, scored as a set, with a
  bonus for engaging distinct carboxylates. Independent salt bridges contribute
  more protonation-linked free energy than several contacts to one.
- **Builds the side chain explicitly.** Each candidate gets 18 rotamers placed
  from ideal internal coordinates by NeRF, scored for clash against the rest of
  the complex and for imidazole-nitrogen to carboxylate-oxygen distance near
  3.0 A. The best rotamer's coordinates go into the structure the pKa predictor
  reads, rather than relying on whatever a fold predictor happened to produce.
- **Stages the search by cost.** Geometric scoring is instant and ranks all
  combinations; the pKa predictor runs only on a shortlist; re-prediction runs
  only on what survives that.
- Runs on all accepted designs rather than only those passing developability
  filters, which also widens backbone coverage.

## How to run

From the challenge directory, with the project environment active:

    python graft_v2.py

Outputs land in `outputs/`: candidate positions with burial and rotamer
quality, every combination scored, the shortlist with sequences, and a
ColabFold-format CSV for re-prediction.

Then re-predict `colabfold_input_v2.csv` (AlphaFold2-multimer, MMseqs2 MSA,
1 model, 3 recycles), score the results, and feed them through the existing
selection.

## What would count as success

Three-histidine combinations reaching past -2 kcal/mol. That crosses from
roughly four-fold to twenty-fold, which is a different claim about the design.

If everything still clusters near -1 kcal/mol, the ceiling is structural rather
than a search problem, and the next lever is pH-aware generation rather than
post-hoc placement: a differentiable surrogate mapping local structural
environment to predicted pKa shift, trained on experimental pKa data, used as a
guidance term during backbone generation. PROPKA is rule-based and not
differentiable, which is why that cannot be bolted on.

## Known limitations of this approach

- Rotamers are placed on a fixed backbone with no repacking of neighbours; a
  real substitution would perturb the surrounding side chains.
- Clash scoring is a hard distance cutoff, not an energy function.
- The pKa predictor remains unvalidated on de novo designed interfaces, so the
  objective being optimised is itself uncalibrated.
- Negative design against pH 7.4 is discrete here, by enumeration and scoring.
  Gradient-based negative design is not possible without the differentiable
  surrogate above.
