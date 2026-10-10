# Challenge 02 — conditional TNF-alpha binder

Anthropic x Adaptyv protein design competition, round 2.

Target: human TNF-alpha soluble trimer (PDB 1TNF). Objectives, in the order the
challenge ranks them: pH-selective binding (bind at pH 7.4, no detectable
binding at pH 6.0), mouse cross-reactivity, affinity.

**Outcome: no design passed in-silico validation.** AF2 folds the binders but
does not dock them, so the pH switch was never tested on a design. The scoring
function built for it *was* validated against known systems (section 4), where
it recovered the FcRn switch histidines blind. 20 designs were submitted anyway,
with paired controls, and the numbers and post-mortem are below.

---

## 1. Why not histidine grafting

The standard approach is: design a binder, look at the interface, mutate chosen
residues to histidine, hope. Histidine is used because its side-chain pKa
(around 6.0 to 6.5) is the only one that falls inside the 6.0 to 7.4 assay
window. Every other titratable side chain is either fully charged or fully
neutral across that entire range and contributes nothing.

| Residue | Side chain | pKa |
|---|---|---|
| Asp | carboxylate | ~3.9 |
| Glu | carboxylate | ~4.3 |
| **His** | **imidazole** | **~6.0-6.5** |
| Cys | thiol | ~8.3 |
| Tyr | phenol | ~10.1 |
| Lys | amine | ~10.5 |
| Arg | guanidinium | ~12.5 |

The problem with grafting is not the chemistry. It is the energy budget.

From Wyman linkage, over a window of width dpH the maximum achievable change in
binding is

    d(log Ka) <= d(nu_H+) * dpH

where d(nu_H+) is the number of protons released on binding. The assay window is
1.4 pH units, so each ideally placed titratable site contributes at most 1.4 log
units — about 25-fold, and that is a ceiling, not a realistic value.

| Sites | Ceiling |
|---|---|
| 1 | ~25x |
| 2 | ~630x |
| 3 | ~16,000x |

A therapeutic switch needs 100 to 1000-fold. **One histidine cannot reach that
even in the ideal limit. Three coupled sites is a floor, not a preference.**

The second failure mode is stability. Burying a charge in a low-dielectric
protein interior costs desolvation energy of several kcal/mol, comparable to the
entire folding stability of a designed miniprotein (typically 3-8 kcal/mol). A
histidine mutated into a finished interface, which protonates at pH 6.0 with no
charge partner, can cost more than the protein has.

---

## 2. Approach: anchor the switch on the target, not the binder

Titratable and charged residues already present on the target are free. They are
already folded and already stable. The target absorbs the electrostatic cost; the
binder is designed around them rather than mutated afterwards.

### 2.1 Map the reachable epitope

`src/target_switch_anchors.py`

In a trimer, the A/B contact set is dominated by the large **buried**
protomer-protomer packing interface, which no binder can reach. A first pass
using the centroid of all A/B contacts returned core residues at relative SASA
0.00 to 0.09 — unreachable.

Corrected: filter contacts to solvent-exposed (relative SASA >= 0.15), then
measure distance to the **nearest** exposed boundary residue rather than to a
centroid, because the exposed boundary wraps around the protomer junction and
is not a compact spot.

Human (UniProt P01375, residues 77-233) was aligned to mouse (P06804) and
conservation mapped onto 1TNF numbering.

Result: 14 titratable candidates within 12 A of the exposed groove.

### 2.2 Score conservation on chemistry class, not residue identity

Strict identity scoring called every usable anchor non-conserved, because
several positions are human Glu / mouse Asp. Both are carboxylates with nearly
identical intrinsic pKa — the switch chemistry survives the substitution.
Histidine has no substitute in this window, so His to Tyr destroys it.

Scoring on chemistry class took the usable anchor count from **0 to 7**.

The one histidine on the exposed epitope, **A73, is His in human and Tyr in
mouse**. A textbook design built around the target's own histidine would have
failed objective 2 silently — the design would look fine and simply not
cross-react.

### 2.3 The switch has to run on cations, not carboxylates

The epitope's titratable character is almost entirely carboxylate. Those sit at
pKa ~4, fully deprotonated at both 6.0 and 7.4. They never switch — they are a
static negative patch.

And the obvious pairing runs backwards. A binder histidine facing a target
carboxylate:

- pH 7.4: His neutral, hydrogen bond, moderate
- pH 6.0: His+ ... Glu- salt bridge, **stronger**

That binds harder in the endosome, which is the opposite of the requirement.

Linkage says what is needed: d(nu_H+) < 0, net proton **release** on binding,
which requires binding to **lower** the pKa of some group. A pKa drops when the
protonated form is destabilised by its environment. The protonated form of His
is +1. So the switch needs **target-side positive charge**.

### 2.4 The cation triad

Conserved positive residues on the exposed boundary, in 1TNF numbering:

| Target | rel. SASA | Mouse | Role |
|---|---|---|---|
| A103 ARG | 0.33 | Lys — charge kept | switch cation |
| B112 LYS | 0.23 | identical | switch cation |
| B103 ARG | 0.60 | Lys — charge kept | switch cation |

Separations: A103-B112 8.6 A, A103-B103 12.2 A, B112-B103 20.2 A.
Maximum span **20.2 A**, coverable by one face of a 70-95 residue binder.

A103 is on one protomer; B103 and B112 on the other. The site is a genuine
**composite inter-protomer epitope**, which is the epitope class the challenge
asks for. It fell out of the charge analysis rather than being imposed.

**Mechanism.** The binder presents three histidines, one facing each cation.

- pH 7.4: His neutral, hydrogen bonding and packing, binder binds
- pH 6.0: each His protonates to +1, now facing a +1 target residue at 4-6 A,
  three like-charge repulsions, binder releases

d(nu_H+) = -3, which clears the 100-1000x requirement with headroom.

The carboxylate anchors (A53, A104, B104, B10, B110, B107, A127) are **affinity
anchors** for binding contacts. The switch rides on the cation triad. Separating
those two roles is what distinguishes this from grafting.

### 2.5 Switch geometry

`src/make_switch_motif.py`

For each cation, the outward solvent normal was computed and an ideal histidine
ring centroid placed along it at a 5.0 A standoff — close enough for a neutral
imidazole to hydrogen bond, close enough for a protonated one to carry real
repulsion.

All three sites cleared the protein surface at the nominal standoff with no
adjustment: clearances 4.57, 5.00, 3.79 A.

**Known weakness.** Net charge within 8 A of A103 and B112 is 0, and -2 to -3
within 12 A. Carboxylates in the wider neighbourhood partially screen the
cationic patch. The local 5 A pair is strongly positive; the medium-range field
is not. A protonated histidine feels both. This is the most likely source of
switch attenuation and it was known before the run.

---

## 3. Pipeline and results

    1TNF  ->  anchor analysis  ->  switch geometry  ->  RFdiffusion
          ->  switch coverage filter  ->  ProteinMPNN + AF2  ->  linkage scoring

### 3.1 Backbone generation

RFdiffusion, binder design mode, target `target_AB.pdb` (chains A and B,
residues 6-157 each, contiguous).

    contigmap.contigs=[A6-157 B6-157 85-85]
    ppi.hotspot_res=[A103,B112,B103]
    diffuser.T=50

Two arms:

| Arm | Checkpoint | Designs | Rate |
|---|---|---|---|
| base | Complex_base | 4 (probe) | 1.9 min/design |
| beta | Complex_beta | 32 | 4.1 min/design |

A cropped target was built and **discarded** — cropping fragmented the chains
into 20 disconnected segments, which RFdiffusion cannot scaffold against. Full
contiguous chains are required.

### 3.2 Switch coverage filter

`src/switch_positions.py`

Identifies which binder positions face each cation, requires distinct positions
per site, and discards backbones that cannot reach at least two sites — before
any sequence design runs. This is the main compute saving in the pipeline.

| Arm | 3 sites | 2 sites | Pass rate |
|---|---|---|---|
| beta (32) | 3 | 9 | 38% |
| base (4) | 3 | 1 | 100% |

Measured helix content of the base-arm binders was 89, 93, 94 and 99 percent —
near-pure helical bundles. The beta checkpoint was enabled for the main run to
get loops and sheets capable of holding three side chains at defined positions.
On this sample it did **not** improve switch coverage. Four designs is not a
meaningful sample, but it does not support the hypothesis either.

### 3.3 Sequence design and validation

ProteinMPNN plus AlphaFold2, 12 surviving backbones, 4 sequences each = **48
sequences**. `initial_guess=True`, `num_recycles=3`, `use_solubleMPNN=True`,
cysteine excluded.

| Metric | Median | Best | Target |
|---|---|---|---|
| i_pTM | 0.130 | 0.161 | > 0.5 |
| i_PAE | 27.9 A | 27.2 A | < 10 A |
| self-consistency RMSD | 36.7 A | 13.4 A | < 2 A |
| pLDDT | 0.735 | 0.892 | > 0.8 |

Verified independently of the scores: in the AF2 models the binder sits near the
target but makes **0 to 12 atom contacts under 5 A**. A real interface has 100 to
300.

**The binders fold. They do not dock.**

`initial_guess` soft-initialises AF2 with the designed coordinates, biasing the
prediction toward the intended pose. It is the standard recommendation for
binder design, and it means these RMSD numbers are more permissive than a blind
refold. They still failed.

---

## 4. Validation of the pH-selectivity scorer

No design docked, so the switch was never scored on a real binder. The scorer
itself was instead characterised against known systems.

`src/linkage_score.py` runs PROPKA on the complex, the free binder and the free
target, computes d(nu_H+) at each pH, integrates the linkage relation across the
window, and reports fold-selectivity plus the titratable groups responsible.

### 4.1 Positive control — FcRn / IgG Fc (PDB 1I1A)

The canonical histidine-driven pH switch. FcRn binds IgG at pH 6.0 in the
endosome and releases it at pH 7.4 in plasma — the inverse direction to this
challenge, which also tests the sign convention.

Result: `delta_log_Ka` = -1.02 (about 10-fold tighter at pH 6.0), 3 effective
sites. Driving residues, found without any prior information:

| Residue | pKa free | pKa bound | d pKa |
|---|---|---|---|
| C436 HIS | 7.06 | 8.56 | +1.50 |
| C310 HIS | 6.42 | 7.06 | +0.64 |
| C435 HIS | 6.47 | 6.71 | +0.24 |

Chain C is the Fc. **His310 and His435 are the known FcRn-binding histidines**,
with His436 in the same cluster. Their pKa is upshifted on binding, so the
complex stabilises the protonated form, so binding is favoured in acid. The
mechanism was recovered correctly from a 600-residue structure.

A135 GLU on the FcRn side contributes in the opposite direction, consistent with
the acidic residues known to pair with those histidines.

Caveats: the real switch is >100-fold, so the magnitude is roughly two orders
conservative — expected from continuum electrostatics on a single crystal
structure with no conformational change. And the `B1 N+` contribution is the
free N-terminus of beta-2-microglobulin, an artifact of splitting a chain that
is never actually free.

### 4.2 Interface control — barnase / barstar, cognate pair (1BRS, A + D)

Chosen as a negative control and **it is not one**. Result: -0.52, 2 sites.

| Residue | d pKa | Contribution |
|---|---|---|
| A73 GLU | +1.53 | +0.501 |
| A102 HIS | -1.04 | -0.340 |
| D35 ASP | +1.48 | +0.200 |

Barnase Glu73 and His102 are its catalytic pair; barstar Asp35 is part of the
acidic loop that inserts into the barnase active site. This is an
electrostatically steered interface with titratable residues at its centre, so a
few-fold pH dependence is plausible rather than spurious. Bad control choice,
not a tool failure.

Worth noting: Glu73 (+0.501) and His102 (-0.340) pull in opposite directions and
partially cancel. That is the buffering effect — a charge network around a
titratable site damps the pH response instead of amplifying it.

### 4.3 True negative control — barnase / barstar, non-cognate pair (1BRS, A + E)

Same protein, same method, no interface.

Result: fold-selectivity **1.00, zero effective sites, empty contributor list.**

### 4.4 Characterisation

| System | fold | sites |
|---|---|---|
| FcRn / Fc | 0.1 | 3 |
| barnase / barstar, cognate | 0.3 | 2 |
| barnase / barstar, non-cognate | 1.0 | 0 |

Correct sign on a known switch, correct residues identified blind, no signal
where there is no interface, magnitude conservative by about two orders.

**Valid as a ranking signal. Not valid as an absolute predictor.** Which is how
it was specified to be used.

---

## 5. Submission

`src/build_submission.py`

20 designs, the Track 2/3 limit. 16 switch variants and 4 paired controls — the
same sequences without the histidine substitutions, so that if anything binds,
"the design binds" is separable from "the histidines did something". Without
controls a positive result would be uninterpretable.

Selection order:

1. switch-site count first. A design engaging two sites cannot reach the
   required selectivity regardless of how well it folds.
2. self-consistency RMSD second, as a weak tiebreak. Every value is a failure
   between 13 and 56 A, so it carries little information here.
3. at most 2 sequences per backbone, for diversity.
4. designs whose switch positions fall on the chain termini excluded — those
   reflect a disordered tail passing near the target, not a structured contact.

Composition: 4 designs with three histidines, 12 with two, across 9 distinct
backbones. All 85 residues, single chain, de novo.

**Known weakness of the submitted set.** 12 of 16 switch designs carry only two
histidines, which by the linkage ceiling caps them near 630-fold, below the
upper end of the 100-1000x requirement. The four three-site designs clear it on
paper but carry the worst fold metrics in the set (RMSD 37-56, pLDDT 0.53-0.63,
against 13-29 and 0.61-0.89 for the two-site designs). That tension is
unresolved and there is no data here to resolve it.

---

## 6. Post-mortem

**The sample was far too small.** Published RFdiffusion binder campaigns report
in-silico pass rates of a few percent and filter from tens of thousands of
sequences. This run generated 48. At a 2 percent pass rate the expected number
of successes from 48 is about one, and zero is the single most likely outcome.
Nothing about the result is evidence against the epitope or the switch design —
it is the expected outcome of undersampling.

The decision to generate 4 sequences per backbone instead of 32 was made to fit
a same-day schedule. That traded the one variable that determines whether
anything validates.

**The target is also hard.** A shallow composite groove between two protomers is
a harder docking problem than a flat domain face. That is what made it the
interesting part of the challenge and it is also why it failed.

**What the result does not tell us:** whether the pH switch works. No design got
far enough to be scored on it.

---

## 7. Environment notes

The RFdiffusion Colab notebook is pinned to 2023 dependencies and Colab now runs
Python 3.13. Three fixes were needed:

1. **DGL.** Modern DGL publishes no Python 3.13 wheels, so `pip install dgl`
   resolved backwards to **0.1.3 (2018)** — a different library generation, not
   an API rename. Shimming is not possible.
   Fix: RFdiffusion runs as a **subprocess** via its shebang, so it does not have
   to share the kernel's Python. Built a Python 3.11 environment with uv
   (dgl 2.5.0+cu121, torch 2.5.1+cu121) and rewrote the shebang of
   `run_inference.py` to point at it. The notebook kernel stays on 3.13.
2. **Kernel-side import.** The setup cell imports `inference.utils.parse_pdb`,
   which drags DGL into the 3.13 kernel and bypasses the shebang fix. Replaced
   with a local PDB parser.
3. **JAX.** `jnp.clip` dropped the `a_min`/`a_max` keyword names. ColabDesign's
   bundled AF2 still uses them. Injected a compatibility shim into
   `designability_test.py`, which also runs as a subprocess.

The subprocess/shebang trick generalises to any pinned research code that runs
its model as a script.

## 8. Tooling bugs caught

Two bugs in this repo's own analysis scripts would each have produced a
confidently wrong conclusion:

1. **Coordinate frame.** RFdiffusion recentres output (~67 A translation here)
   and renumbers chains to 1-N. Motif coordinates stored in the original frame
   matched nothing, so every design reported zero switch coverage. Fixed by
   superposing the design's target chains onto the reference and transforming
   the motif into that frame (alignment RMSD 0.13 A).
2. **Virtual CB.** RFdiffusion emits poly-glycine backbones. Glycine has no CB,
   so the side-chain position function returned None for every residue and
   coverage was again zero. Fixed by constructing CB from N, CA, C.

A third was an over-count: the same binder residue was credited with covering
two switch sites 8-20 A apart. Fixed by enforcing one-to-one assignment, which
took beta-arm coverage from 14 to 12 backbones.

---

## 9. What transfers

1. The linkage ceiling. 1.4 log units per switching proton over a 1.4 pH-unit
   window. Any pH-switch design can be sized against this in seconds, and a
   single-histidine design can be rejected before it is built.
2. A linkage-based scorer validated on FcRn/IgG recovers the known switch
   residues blind and returns nothing on a non-interface pair. It ranks; it does
   not predict absolute selectivity.
3. Switch chemistry can be anchored on the target rather than the binder. The
   target absorbs the desolvation cost.
4. Cross-species conservation of a switch must be scored on chemistry class, not
   residue identity. Strict identity called every usable anchor on this epitope
   non-conserved.
5. Filtering backbones on switch-site reachability before sequence design is
   cheap and removes most of the wasted compute.
6. Sample size is the variable that decides whether anything validates. Do not
   trade it for schedule.

## 10. Open problem

No structure-prediction tool in this stack models protonation. AF2 has no
concept of a titratable group; RFdiffusion's backbone frames carry no charge.
Proton-PottsMPNN (Jacobsen, Ovchinnikov et al., 2026) addresses the sequence
design side by representing protonated and deprotonated His, Asp and Glu as
distinct tokens, and recovered 237 and 288 pH-dependent binders experimentally —
but transition pH values landed between 4.0 and 5.8, with nothing at 6.2.

Producing a switch is solved. **Placing the transition at a chosen pH is not.**
That reduces to predicting pKa in designed rather than natural environments,
where existing predictors are weakest, and the dominant determinant is burial —
which solvent-accessible surface area describes poorly, since it cannot
distinguish a shallow dish from an enclosed cavity.

---

## 11. Files

| Path | Contents |
|---|---|
| `src/target_switch_anchors.py` | epitope mapping, conservation, anchor ranking |
| `src/make_switch_motif.py` | switch geometry from the cation triad |
| `src/switch_positions.py` | coverage filter and ProteinMPNN His bias |
| `src/linkage_score.py` | PROPKA-based pH-selectivity ranking |
| `src/build_submission.py` | submission assembly, switch placement, controls |
| `challenges/02-tnfa/data/` | 1TNF.cif, target_AB.pdb, human and mouse FASTA |
| `challenges/02-tnfa/motif/` | switch_motif.json, switch_motif_sites.pdb |
| `challenges/02-tnfa/results/` | anchors.csv, switch reports, bias file |
| `challenges/02-tnfa/results/control_*.csv` | scorer validation, per-residue contributions |
| `challenges/02-tnfa/submission20.csv` | the 20 submitted designs with metadata |

Reproduce:

    python src/target_switch_anchors.py --pdb challenges/02-tnfa/data/1TNF.cif \
      --human-fasta challenges/02-tnfa/data/human_tnf.fasta \
      --mouse-fasta challenges/02-tnfa/data/mouse_tnf.fasta \
      --shell 12.0 --out challenges/02-tnfa/results/anchors.csv

    python src/make_switch_motif.py --pdb challenges/02-tnfa/data/1TNF.cif \
      --out-prefix challenges/02-tnfa/motif/switch_motif

Expected: 14 candidates, 7 switch-conserved, A73 flagged SWITCH LOST; then three
sites with clearances 4.57 / 5.00 / 3.79 and maximum span 20.2 A.
