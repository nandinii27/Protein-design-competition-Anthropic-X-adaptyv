# Methods — Challenge 01, conditional EGFR binder

De novo single-chain minibinders against domain III of the human EGFR
extracellular region, designed for preferential binding at pH 6.5 over pH 7.4
and for retained binding to mouse EGFR. No existing binder was used as a seed,
template or starting point at any stage.

## The problem this pipeline is built around

No available structure predictor takes pH as an input. Folding a binder-target
complex at pH 6.5 and at pH 7.4 returns the identical structure, so a
conventional generation stack produces a signal for affinity, a signal for
cross-species binding, and nothing at all on the highest-weighted objective.

The approach injects pH awareness in two places. Upstream, by restricting the
epitope to the intersection of conserved, acidic and solvent-exposed surface,
so that a protonation-dependent salt bridge has somewhere to form and does so
on a surface present in both orthologues. Downstream, by grafting histidines at
geometrically identified positions, re-predicting the mutants, and scoring them
with a thermodynamic linkage function over predicted pKa shifts.

## 1. Epitope selection

Three constraints were required simultaneously:

- **conserved** between human (UniProt P00533) and mouse (UniProt Q01279) EGFR
- **acidic**, providing a carboxylate partner for a protonated histidine
- **solvent exposed**, so a designed binder can reach it

Human and mouse full-length sequences were aligned pairwise (BLOSUM62, gap open
−11, extend −1). The structure chain was aligned separately to the human
sequence so that conservation could be mapped onto author residue numbering.
Relative solvent accessibility was computed with Shrake-Rupley on the full
assembly and normalised by residue-type maximum accessible area, with an
exposure cutoff of 0.25.

Structure: PDB 6ARU, chain A, the EGFR extracellular region. Chains B and C of
that entry are a cetuximab Fab and were removed before any accessibility
calculation; leaving them in would have reported the antibody-covered surface
as buried.

Within the domain III window (residues 310–480, mature-protein numbering), 171
residues were evaluated. Seven satisfied all three constraints:

    Glu320, Asp323, Asp364, Glu367, Glu400, Glu431, Glu472

Candidate patches were built around each anchor using a 12 Å CA radius and
ranked by anchor count, then by the conserved fraction of the patch. Patches
with fewer than two anchors were discarded: a single salt bridge contributes
too little protonation-linked free energy to switch binding.

Three viable patches emerged, all fully or near-fully conserved:

| centre | anchors | residues | conserved fraction |
|---|---|---|---|
| Glu320 | Glu320, Asp323 | 11 | 1.000 |
| Glu400 | Glu400, Glu431 | 6 | 1.000 |
| Glu367 | Asp364, Glu367 | 11 | 0.818 |

Glu320 and Glu367 patches share residues 330, 333 and 334 — they are adjacent
regions of one continuous conserved acidic stretch, not independent sites. The
CA–CA distance between Glu320 and Glu367 is 18.4 Å, within the span of a
65–100 residue minibinder, so the union was targeted to give the designer four
carboxylates rather than two.

Final epitope, 19 residues:

    314, 317, 318, 319, 320, 322, 323, 325, 330, 333, 334,
    336, 362, 364, 365, 366, 367, 369, 372

### Relationship to the cetuximab epitope

Cetuximab contacts (any heavy atom within 4.5 Å in 6ARU) are residues 349–473.
None of the three conserved acidic patches overlaps that footprint. This is not
coincidence: the cetuximab surface differs between human and mouse, which is
why cetuximab fails on mouse EGFR, so conserved residues are necessarily
elsewhere. The conservation and chemistry constraints push the epitope off the
canonical blocking surface by construction.

Distance from each patch centre to the nearest cetuximab contact residue:
Glu320 18.5 Å, Glu367 14.6 Å, Asp364 10.2 Å. A minibinder of this size buries
roughly 800–1000 Å² and spans 15–20 Å, so partial occlusion of the ligand site
is plausible at these offsets but is not established. Functional blocking is
therefore a weaker claim for these designs than for a cetuximab-epitope binder,
and this was accepted as the cost of satisfying the cross-species objective.

## 2. Generation

Backbone generation and sequence design: BindCraft, four-stage protocol,
default design, prediction, interface and template protocols, relaxed filters.
Target: the domain III window trimmed to residues 300–490 (191 residues), which
reduces the cost of backpropagating through the network and avoids spending
capacity on domains the binder never contacts. Binder lengths sampled from
65–100. Hotspots as listed above.

18 designs were accepted across 10 distinct backbones (generator seeds).
Hotspot RMSD for accepted designs ranged 1.2–3.0 Å, confirming the binders
engage the intended patch rather than drifting elsewhere on the target.

Relaxed rather than default filters were used deliberately. Default filters are
calibrated for experimental binding success, which is the lowest-weighted
objective here, and two further orthogonal filters are applied downstream that
the generator knows nothing about. Pool size was prioritised over per-design
binding quality, consistent with the stated ranking.

## 3. pH scoring

For each complex, pKa values were predicted twice with PROPKA3: once on the
isolated binder chain, once on the two-chain complex. Sites were restricted to
ionisable binder residues with any heavy atom within 6.0 Å of the target.

For each such site, the protonation contribution to binding free energy at a
given pH follows the standard linkage relation:

    ddG(pH) = -RT * ln[ (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH)) ]

Contributions are summed over interface sites, and the reported quantity is

    pH_selectivity = ddG(6.5) - ddG(7.4)

Negative values indicate binding more favourable at the lower pH.

Two filters are applied before any ranking:

1. **Geometric.** At least one binder histidine must have a side-chain nitrogen
   (ND1 or NE2) within 4.0 Å of a target carboxylate oxygen in the predicted
   complex. A design without this has no structural mechanism for pH switching
   regardless of its score.
2. **Sign.** Designs with a non-negative selectivity score are discarded rather
   than down-ranked, since a positive score indicates the inverse switch.

### Result on the generated designs

Of the 18 accepted designs: 6 had negative selectivity, 3 had a histidine
contacting any target carboxylate, and 1 had a histidine contacting one of the
four conserved anchors. The intersection of correct sign and histidine-anchor
contact was **empty**.

This is the expected outcome of a generator with no pH objective, and it is
reported. Magnitudes were also small — the largest was −0.88 kcal/mol,  
roughly a four-fold affinity difference, where a useful
tumour-selective switch would want ten- to hundred-fold.

One incidental observation: ProteinMPNN placed up to five histidines at some
interfaces without being asked to. Selecting an epitope on conserved acidic
surface appears to enrich interface histidine as a side effect, since histidine
is a natural hydrogen-bond partner for Asp and Glu. Those histidines were not,
however, positioned against the chosen anchors.

## 4. Histidine grafting

Since no generated design combined the correct sign with a mechanism, histidine
was placed deliberately.

For each design passing the liability filters (8 of 18), every binder position
was evaluated against the four conserved anchors:

- CB within 4.5–10.0 Å of an anchor carboxylate oxygen. The window follows
  histidine geometry: CB to NE2 is about 3.5 Å and an imidazole-carboxylate
  hydrogen bond wants N···O near 3 Å. Closer and the ring clashes; further and
  it cannot reach under any rotamer.
- The CA→CB vector pointing toward that oxygen, rejecting positions whose side
  chain projects into the binder core.
- Glycine, proline, cysteine and existing histidine excluded.

Positions were scored as cosine divided by distance. Single mutants were
generated at the best positions and pairs at two positions, giving 35 unique
mutants from 8 designs (21 single, 14 double).

Mutants were re-predicted as complexes with ColabFold (AlphaFold2-multimer v3,
MMseqs2 UniRef+Environmental MSA, 1 model, 3 recycles, no relaxation). Mutating
the sequence and re-folding, rather than editing a side chain in place, gives a
structure predicted under the same procedure as the parent; an unrelaxed
in-place mutation would produce geometry that the pKa prediction cannot be
trusted on.

### Result on the mutants

Of 35 re-predicted mutants: 14 improved on their parent's selectivity score, 11
had negative selectivity, 2 placed a histidine on a conserved anchor, and 10
cleared geometry, sign and liability filters together.

The delta against the parent is the quantity that isolates the graft. An
absolute score mixes the grafted histidine with whatever incidental protonation
effects the parent already carried.

Only 2 of the 10 contact a *conserved* anchor; the remaining 8 pass the
geometric filter via a carboxylate elsewhere on the target. Those support the pH
objective but not necessarily the cross-species objective, since the contacted
residue may not be conserved. This distinction is preserved in the submission
metadata as `n_his_on_conserved_anchor`.

## 5. Developability

Applied as hard filters, not ranking signals: length 60–110, zero free
cysteines, maximum contiguous hydrophobic run of 5, largest exposed hydrophobic
cluster below 600 Å² (clustered on exposed hydrophobic residues within 8 Å and
summed over SASA), and absence of N-glycosylation sequons and NG, DG and DP
motifs.

## 6. Selection

Three objectives were ranked by non-dominated sorting, not combined into a
weighted sum: the brief states a priority order, and a priority order does not
define exchange rates between objectives.

- pH selectivity, minimised
- interface confidence (ipTM), maximised
- histidines contacting a conserved anchor, maximised

The third is a proxy for cross-species binding. A histidine on a conserved
carboxylate supports low-pH binding on a surface present in both orthologues.

A global cap of two designs per backbone was applied. Variants of one backbone
differ by a residue or two and measure nearly the same thing, so without the cap
the submission would have spent nine slots on one protein.

The submission contains 18 designs across all 10 backbones:

- **5 switching designs** (3 grafted, 2 parents) with negative predicted
  selectivity, ranging −1.21 to −0.02 kcal/mol
- **13 affinity controls**, spanning predicted selectivity −0.88 to +0.58

The controls are deliberate. No prior data establishes whether a PROPKA-derived
linkage score tracks measured pH dependence on de novo designed interfaces. A
submission of switching designs alone cannot distinguish "the method works" from
"the method has no discriminating power", because there is nothing to compare
against. Designs spanning the predicted range turn the assay into a calibration
of the score itself, so that a null result on the switching designs is still
informative. Each row is labelled by role in the `design_role` column.

## 7. Limitations

- PROPKA is an empirical model parameterised largely on natural proteins. Its
  accuracy on de novo designed interfaces is unknown. The linkage score derived
  from it has no prior experimental validation in this setting and should be
  read as a hypothesis-generating ranking, not a calibrated prediction.
- The linkage relation as applied treats titratable sites as independent.
  Coupled titration between neighbouring ionisable residues is not modelled.
- Predicted selectivity magnitudes are small. The best design at −1.21 kcal/mol
  corresponds to roughly an eight-fold affinity difference between the two pH
  conditions, below what would be needed for tumour selectivity in vivo.
- Each complex is a single static structure from a single model with no
  relaxation. Grafted histidine side-chain geometry is AlphaFold's unrefined
  prediction, and the pKa calculation depends directly on it.
- Protonation-dependent conformational change is not represented anywhere.
- Mouse cross-reactivity is argued structurally — the epitope is conserved by
  construction — and not predicted. No mouse complex was folded.
- Functional blocking of ligand binding is inferred from proximity to the
  cetuximab epitope (10–18 Å), not demonstrated.
- Self-consistency between designed and re-predicted backbones was assessed only
  through the generator's own filters, not independently for the mutants.

## 8. Reproducibility

Code: see repository. Epitope selection, pH linkage scoring, liability
filtering and selection are implemented in `src/binderkit`; the per-stage driver
scripts are in `challenges/01-egfr`.

Inputs: PDB 6ARU; UniProt P00533 (human EGFR) and Q01279 (mouse EGFR).

Generation: BindCraft, four-stage protocol, relaxed filters, settings JSON in
the repository. Re-prediction: ColabFold AlphaFold2-multimer v3, single model,
3 recycles, MMseqs2 UniRef+Environmental, no relaxation.

Intermediate tables (epitope patches, parent pH scores, graft candidates,
mutant metrics, final selection) are included in the repository.
