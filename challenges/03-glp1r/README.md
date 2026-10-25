# Challenge 03 — GLP1R agonist

Design a binder to human GLP1R that also **activates** it. Binding is scored, but a
design that binds without activating does not meet the second objective.

Submitted: 20 de novo designs, 38–59 residues.

---

## Why this challenge is not a binder problem

GLP1R is a class B1 GPCR with two domains: a small extracellular domain (ECD) sitting
on a seven-helix transmembrane bundle (TMD). The natural hormone is one long helix
doing two separate jobs:

- its **C-terminal half** docks in the ECD groove — this supplies essentially all the
  affinity;
- its **N-terminal residues** insert into the TMD core — this is what activates, by
  opening the intracellular face so Gs couples and cAMP rises.

Remove the N-terminus and you get exendin(9-39) or GLP-1(9-36): tight groove binders
that are pure antagonists. The challenge page recommends targeting the ECD groove,
and following that literally produces the antagonist architecture. The 2025 Baker-lab
de novo GPCR miniprotein work got antagonists at GLP1R, GIPR, GCGR and CGRPR for this
reason.

Everything below exists to select for the second objective, because the first takes
care of itself.

## Approach

**Generation.** Hotspots were chosen to span both domains — five TMD-pocket residues
(190, 234, 241, 310, 364) and three ECD-groove residues (68, 89, 128), taken from the
native contact set in PDB 6X18. A design satisfying all eight has to bridge the full
distance, so the constraint itself encodes the activation mechanism. An ECD-only
binder cannot satisfy the first five.

Two tracks were run: all-atom structure-conditioned design (96 designs) and a
sequence-only masked language model conditioned on the ECD (120 peptides).

**Scoring.** Three signals, in descending order of how much they are trusted:

1. `pocket_contact_recovery` — fraction of the native agonist's TMD-pocket contacts
   the design reproduces. Validated on controls (below).
2. `insertion_depth` — projection of the design's deepest Cα onto the pocket axis,
   normalised so 1.0 is native depth.
3. `fgw_distance` — Fused Gromov-Wasserstein distance between the design's interface
   contact graph and the native agonist's. Nodes are contacting binder residues,
   edges are normalised Cα distances, node features are pocket depth, hydropathy and
   charge. FGW is used because no residue correspondence exists between a 30-residue
   hormone and a 50-residue design. **Exploratory, not validated.**

Hard filters run first and disqualify regardless of score: length bounds, GRAVY
(the receptor is supplied in DDM/CHS micelles, so exposed hydrophobic surface
partitions into detergent), long hydrophobic runs, free cysteines (the ECD carries
three native disulfides). A novelty filter rejects anything above 40% identity to
GLP-1 or exendin-4.

## Result 1 — interface confidence cannot detect agonism

Four controls with established pharmacology were co-folded against the receptor and
ranked by interface confidence, the metric the standard pipeline selects on:

| control | ipTM | truth |
|---|---|---|
| GLP-1(9-36) | 0.962 | antagonist |
| GLP-1(7-36) | 0.957 | agonist |
| exendin(9-39) | 0.941 | antagonist |
| exendin-4 | 0.929 | agonist |

Perfectly interleaved, antagonist on top. All four score above 0.92, which reads as a
confident, well-predicted interface. Ranking designs by ipTM selects antagonists with
high confidence.

## Result 2 — pocket contact recovery does

Same four molecules, measured geometrically:

| control | truth | pocket recovery | groove recovery | insertion depth |
|---|---|---|---|---|
| GLP-1(7-36) | agonist | **0.93** | 0.68 | 1.01 |
| exendin-4 | agonist | **0.89** | 0.76 | 1.01 |
| GLP-1(9-36) | antagonist | 0.64 | 0.74 | 0.92 |
| exendin(9-39) | antagonist | **0.39** | 0.71 | 0.75 |

Separation margin +0.125, every agonist above every antagonist.

Groove recovery is flat across all four (0.68–0.76) with the antagonists scoring
*higher* than GLP-1(7-36) — every one of these grips the groove equally well, which
is why they all co-fold confidently. Pocket recovery is what separates them.

**Limitation: n = 4.** This licenses the metric; it does not validate it.
GLP-1(9-36) at 0.64 sits closer to the agonists than it should, which would likely
become a real failure at larger n.

## Result 3 — sequence-only generation failed on this target

| track | passed hard filters | rate |
|---|---|---|
| structure-conditioned | 93 / 96 | 97% |
| masked LM, 30-mer | 9 / 60 | 15% |
| masked LM, 39-mer | 0 / 60 | 0% |

The 39-mers collapsed into low-complexity repeats with C-terminal cysteine clusters
(`...SSCLCPC`, `...CCGCCC`); 13 were rejected for free cysteines alone. The 30-mers
died mostly on `X` tokens — unresolved mask positions, not amino acids.

This contradicts the expectation set by the blind cyclic-peptide agonist work at this
receptor (Li, Wiita, Helleday & Bryant 2025), which succeeded sequence-only. The
likely reason is that those peptides were cyclic and chemically synthesised, with no
expression or detergent constraints. Those constraints do real work here.

No PepMLM designs were submitted.

## Selection

Ranked by pocket recovery, with FGW distance as a continuous tie-break, then greedy
diversity selection at a 70% identity ceiling.

Pocket recovery **saturates**: 24 of 93 designs sit exactly at 0.929 and another 37 at
0.893. The hotspot constraint worked well enough to compress the top of the
distribution, so pocket recovery cannot rank within its own ceiling. FGW spread that
tied group from 0.132 to 0.216 and is what orders the final set — which means the
ordering inside the tied group rests on an unvalidated metric.

All signals correlate strongly (pocket recovery with FGW −0.84, with insertion depth
0.82, with ipTM 0.72). They measure overlapping things; they diverge at the top, which
is where selection happens.

**Overlap with a pure ipTM top-20: 4 of 20.** Sixteen of the submitted designs would
not have been chosen by the default pipeline, which would instead have included
designs at 0.750 and 0.786 pocket recovery.

Selected set: mean ipTM 0.720, mean pocket recovery 0.929, mean FGW 0.162.

## What this does and does not show

It shows, on this target, that interface confidence does not distinguish agonists from
antagonists, and that a mechanism-derived geometric metric does on a four-point
control set. That is a well-posed question with preliminary evidence.

It does not show that the submitted designs activate GLP1R. Nothing computational here
measures efficacy. The wet-lab result is the test.

## Known gaps

- The bivalent architecture (ECD-anchored domain plus a flexible pocket-reaching arm)
  was designed but never run — the monthly job quota was exhausted.
- Affinity is measured on a detergent-solubilised construct with undisclosed
  thermostabilising mutations and no G protein; activation is measured on cell-surface
  receptor. Optimising hard against the first drifts away from the active-state pocket
  the second depends on.
- No molecular dynamics; all poses are single static predictions.

## Files

| file | role |
|---|---|
| `targets.py` | sequences, control set, constraints |
| `geometry.py` | pocket frame, insertion depth, contact recovery |
| `fgw_interface.py` | FGW interface comparison |
| `score.py` | hard filters and composite score |
| `novelty_filter.py` | identity guard against natural peptides |
| `measure_controls.py` | the agonist/antagonist gate |
| `analyze_designs.py` | full measurement pass over all designs |
| `designs.py` | generated sequences with reported metrics |
| `design_metrics.csv` | all 96 designs, every measurement |
| `submission.csv` | the 20 submitted |

## Reading

- Muratspahić et al. 2025 — de novo miniprotein agonists and antagonists targeting
  GPCRs. The GLP1R antagonist result is the most relevant prior here.
- Li, Wiita, Helleday & Bryant 2025 — blind de novo dual cyclic peptide agonists for
  GCGR and GLP1R. bioRxiv 2025.06.06.658268.
- Chen et al. 2025 — target sequence-conditioned peptide binder design. Nat Biotechnol.
- Sexton & Wootten (Monash) on class B activation; Tate (MRC-LMB) on thermostabilised
  constructs.
