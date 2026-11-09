# Methods — de novo binders to HLA-A\*11:01 / KRAS G12D `VVVGADGVGK`

## 1. The discrimination problem, measured before anything was designed

The obvious design protocol for a neoantigen target is to place diffusion hotspots on the mutated residue. For this target that protocol is wrong, and the reason is measurable.

Per-residue solvent accessibility was computed on the peptide in 7OW6 (chains A–C, TCR removed). **Asp6 — the residue that is Gly in wild-type KRAS — presents 8.8 Å² of solvent-accessible surface. One carboxylate oxygen is entirely buried; the other exposes 6.4 Å².** For comparison, the therapeutic TCR in the same structure buries **zero** of that side chain across its entire interface.

The mutation is therefore not available for direct contact. A binder cannot read the aspartate the way an antibody reads an epitope. Discrimination has to come from **indirect readout**: the conformation the buried charge imposes on the peptide backbone, the electrostatic field it projects through the groove opening, or the difference in dynamics between the two variants.

Every subsequent decision follows from that measurement.

## 2. Conformational state: designing against the analyte that is actually in the assay

7OW6 is a ternary complex. The peptide conformation in it has been remodelled by TCR engagement, and the assay analyte is free pMHC trimer in solution. Designing against an induced-fit conformation optimises for a state that is not present when the measurement happens.

Five-model AlphaFold2 ensembles were generated for the free G12D trimer and the free wild-type trimer, and each was superposed on the α1-α2 groove Cα atoms only, with the peptide then measured **in place** rather than refitted.

Two results:

**Induced-fit remodelling is real but localised.** Groove superposition RMSD is 0.60 Å. Against that frame, D6 moves 2.2–3.0 Å and G7 moves 1.7–2.4 Å across all five models, while the anchor residues move 0.3–0.9 Å. The bulge is remodelled; the ends are not. Designing against the TCR-bound bulge would have aimed the interface at coordinates that shift by up to 3 Å in the free state.

**The variants differ in rigidity, not only in chemistry.** Inter-model spread across the five predictions, measured over the p4–p7 bulge:

| | bulge spread | anchor spread | bulge/anchor | mean pLDDT |
|---|---|---|---|---|
| G12D | 0.44 Å | 0.18 Å | 2.5 | 72.4 |
| wild type | 0.95 Å | 0.19 Å | 5.1 | 67.3 |

The anchors are identical between variants, which is the internal control that makes the bulge comparison meaningful. The wild-type bulge is **2.15× more variable**, and two of five wild-type models failed to dock the peptide at all against zero of five for G12D. Wild-type KRAS carries glycine at peptide positions 4, 6, 7 and 9; the G12D aspartate anchors a bulge that is otherwise unusually flexible.

This is the discriminating signal that is physically accessible: **the mutant presents a more ordered surface, and conformational selection against that order is a route to specificity that does not require touching the mutation.**

## 3. Design

Backbones and sequences were generated with **BoltzGen** against the complete free-state target: HLA-A\*11:01 α1-α2 (182 residues), β2-microglobulin (99 residues) and the G12D decamer, with binding directed to peptide positions 4–8.

**β2-microglobulin is included in the target and explicitly excluded as a binding partner.** This is not cosmetic. A trimmed α1-α2-plus-peptide target leaves the β2m socket empty, and a generative model will fill it: a run against the trimmed target produced 142 designs with zero peptide contact, 50 of which were more than 80% identical to β2m itself. The fix is to supply the full trimer and declare `not_binding: all` on the β2m chain, which converts an open socket into a negative-design constraint. A probe run after the correction put 5 of 5 designs on the declared hotspots.

Format discipline throughout, for Adaptyv's reducing cell-free expression system: single-chain miniproteins, 55–90 residues, cysteine-free preferred, no disulfide-dependent formats.

## 4. Five independent filters

No single predictor score was trusted. Each design was scored by measurements that fail in different ways.

**On-target interface confidence.** Boltz-2 co-folding of each design against the full trimer, scored with **ipSAE** rather than ipTM — ipSAE is interface-local and is not inflated by a well-predicted target. Reported alongside the minimum interaction PAE.

**Interface geometry.** Buried solvent-accessible surface recomputed from the refolded complexes, split between peptide and MHC. The reported quantity is Å² buried on the *peptide* and its fraction of total buried target surface. This filter alone removes designs that satisfy a cofolding score by lying along the α1/α2 helices.

**Monomer self-consistency.** Each binder was folded **alone**, target absent, with Chai-1 in single-sequence mode. A design that only holds its shape when the target templates it is both an expression risk and a false positive in the cofolding score. This filter is not redundant with ipSAE: it demoted two designs ranked 3rd and 5th on ipSAE, whose monomer pLDDT is 0.665 and 0.725 against 0.87–0.91 for the designs that survived.

**Mutation readout, with a control.** ProteinMPNN was given each designed complex with the peptide masked and asked for the conditional probability of each residue at each peptide position. The same calculation was run on the identical structure with the binder chain deleted. The reported quantity is the shift between them:

```
Δ = log[P(Asp)/P(Gly)]_complex − log[P(Asp)/P(Gly)]_apo
```

The control is what makes the measurement mean anything: the HLA groove already prefers aspartate at p6 on its own, so a high probability in the complex says nothing about the binder. Average P(Gly) at peptide positions 4, 7 and 9 — glycine in both variants — is reported as a structure-quality check on each measurement.

## 5. Result of the readout test

The apo baseline is stable across all 25 structures at log-odds **2.81 to 3.13**: the groove alone favours aspartate over glycine at p6 by roughly 20:1, consistently.

Against that baseline, **no design adds direct readout of the mutation.** Two sit at the zero line (+0.08, +0.02) and both have weak glycine controls (0.59, 0.76), making their backbones the least reliably read of the set. Every design with a glycine control above 0.85 returns a *negative* Δ between −0.35 and −1.40.

This is the outcome the burial measurement in §1 predicts. A binder seated on the TCR-facing surface has no physical route to a side chain presenting 8.8 Å². Reporting it is more useful than a cofolding score that implies otherwise: mutation specificity on this target is not achievable through contact, and any claim of it from a co-folding score alone is unsupported.

The route that remains open is the rigidity difference quantified in §2, which is a property of the ensemble rather than of a single structure, and is addressed in §10.

## 6. Three-target counter-screen

Each of the two leading designs was co-folded against three analytes — the G12D trimer, the wild-type trimer, and the empty groove — within a single predictor, with the on-target fold recomputed rather than carried over from the design stage, since a difference between two predictors is not a difference between two targets. Chai-1 with MSA-backed co-folding (HLA alignment depth 10,483; β2m 6,283), identical settings across all analytes. The leading design's triple was repeated under a second seed.

| design | seed | binder–MHC ipTM, G12D / apo / WT | peptide dependence | G12D − WT, binder–peptide ipTM |
|---|---|---|---|---|
| `_09` | 42 | 0.827 / 0.317 / 0.824 | **−62%** | +0.005 |
| `_06` | 42 | 0.715 / 0.585 / 0.758 | −18% | −0.012 |
| `_06` | 7 | 0.728 / 0.566 / 0.763 | −22% | −0.009 |

### Cross-panel selectivity

Both leading designs were co-folded against ten HLA-A\*11:01 ligands under identical settings. The on-target measurement for `_09` repeated at 0.1122 against 0.1107 in the counter-screen — two sessions, 1.3% apart — which is what makes the panel comparable.

| peptide | `_09` ipTM pep | % of on-target | `_06` ipTM pep | % of on-target |
|---|---|---|---|---|
| VVVGADGVGK on-target | 0.112 | 100 | 0.055 | 100 |
| VVVGACGVGK G12C | 0.095 | 84 | 0.094 | **172** |
| VVVGAVGVGK G12V | 0.083 | 74 | 0.068 | **124** |
| VVVGASGVGK G12S | 0.075 | 67 | 0.065 | **118** |
| VVVGAAGVGK G12A | 0.095 | 85 | 0.061 | **111** |
| VVVGAGDVGK G13D | 0.057 | 51 | 0.029 | 52 |
| AVFDRKSDAK unrelated (EBV) | 0.010 | 8 | 0.010 | 17 |
| VVGADGVGK 9-mer register | 0.014 | 12 | 0.005 | 10 |
| VVVGARGVGK G12R | 0.013 | 12 | 0.005 | 8 |
| IVTDFSVIK unrelated (EBV) | 0.020 | 18 | 0.004 | 7 |

**Neither design is a generic groove binder.** Two unrelated A\*11:01 ligands score at 7–18%, and binder–MHC ipTM against them collapses toward the apo value. Both are **register-specific** — the 9-mer form of the same neoantigen is rejected at 10–12%, so they read one peptide geometry rather than "a bulge" — and both **reject G12R outright** at 8–12%, where the substitution is bulky and cationic.

**The designs separate inside the mutation family, which is where it matters.** `_09` peaks on its own target, with the other codon-12 substitutions at 67–85%. `_06` does not: it is maximal on G12C at 172% of its on-target score and exceeds its design target on four of the five variants. `_06` is a KRAS codon-12 scaffold binder that happens to prefer cysteine at position 6; `_09` is the only design in the set whose highest-scoring ligand is the one it was designed against. G13D at ~52% for both shows that position matters as much as charge — displacing the aspartate by one residue halves the signal.

The accurate description of the lead design is therefore a **KRAS codon-12 neoepitope binder with correct register, correct HLA-ligand selectivity, and a modest on-target preference within the mutation family.** That is a quantified selectivity profile rather than a claim of specificity, and it is the honest ceiling of what indirect readout achieves here.

**Peptide dependence separates the designs sharply.** Removing the peptide costs `_09` 62% of its binder-to-MHC interface confidence; `_06` loses 18–22% and retains a complex ipTM of 0.87 against the empty groove. A binder that grips the apo groove is reading HLA, not neoantigen.

**Mutation specificity is absent in both**, and `_06` scores marginally *better* on wild type. The seed replicate on `_06` reproduces the specificity gap to within 0.005 and peptide dependence to within 0.03, so these are properties of the structures rather than of the sampler.

This is the third independent line of evidence converging on the same conclusion — the solvent-accessibility measurement of §1 predicted it, the ProteinMPNN test of §5 measured it, and the counter-screen confirms it under a different model and a different observable.

## 7. Semiempirical quantum electrostatics

The mutation is a carboxylate buried in a low-dielectric groove, which is the regime where fixed-charge force fields are least reliable: no electronic polarisation, so the surrounding residues cannot respond to the buried charge; no charge transfer; and point charges that cannot represent the quadrupolar character of a carboxylate field at close range. GFN2-xTB carries an explicit anisotropic electrostatics term — second-order density fluctuations with dipole and quadrupole character — which is the physics a point-charge model discards.

Four single-points with ALPB aqueous solvation on a 5 Å truncation around the peptide (1,020–1,251 atoms), link-atom hydrogens on severed backbone bonds, identical selection and identical caps for both variants. The wild-type system is built by deleting the Asp6 side chain and placing a hydrogen along the former CA–CB vector, so the backbone is untouched and the two systems differ by exactly one side chain.

```
ΔΔE = [E_complex(G12D) − E_complex(WT)] − [E_target(G12D) − E_target(WT)]
    = −8301.47 − (−8300.37)  =  −1.10 kcal/mol
```

The free-binder term cancels exactly, which is why only four single-points are needed. The two −8,300 kcal/mol terms are the intrinsic energy of an aspartate side chain against a hydrogen; they are meant to cancel and they do, but the residual is one part in eight thousand and the result should be read with that in mind.

**ΔΔE = −1.10 kcal/mol favours the binder on the mutant**, and sits inside the few-kcal/mol uncertainty that GFN2-xTB interaction energies carry on non-covalent complexes. It is reported as a tie-breaker, not as evidence of specificity. Its value is that it leans the same way as the only other number that did: the counter-screen gave this design +0.005 toward G12D on binder–peptide ipTM. Two weak signals from different levels of theory, on different observables, agreeing in sign.

Four limitations bound what it can carry. Interaction energy is not binding free energy — no conformational or translational entropy, no explicit desolvation. The two systems differ in net charge by one unit, which is the physical signal but also makes implicit solvation the load-bearing approximation, and ALPB's uniform dielectric is weakest exactly at a charge on the protein–water boundary. A single-point on one conformer is blind to the entropic contribution that §2 argues may carry most of the real discrimination. And no published validation exists for xTB-derived specificity scores against wet-lab pMHC binder specificity; the precedents cover affinity correlation on pMHC-II substitution panels and protein–ligand systems.

## 8. Complex molecular dynamics

6.63 ns of implicit-solvent MD on the lead design's complex — the full binder, α1-α2, β2-microglobulin and peptide assembly, 5,702 atoms with hydrogens, Amber ff14SB with OBC2 implicit solvent, hydrogen mass repartitioning at a 4 fs step, 300 K Langevin, on a T4 through OpenMM's OpenCL platform at 377 ns/day.

Discarding the first nanosecond:

| metric | mean | first ns → last ns | trend |
|---|---|---|---|
| binder Cα RMSD | 1.82 ± 0.17 Å | 1.56 → 1.84 | +0.01 Å/ns |
| binder–peptide contacts | 52.9 ± 10.9 | 52.7 → 53.5 | +1.9 /ns |
| binder–MHC contacts | 284 ± 26 | 274 → 283 | −0.74 /ns |
| binder–peptide centre-of-mass | 12.29 ± 0.24 Å | 12.14 → 12.42 | +0.05 Å/ns |

**The designed fold holds**: Cα RMSD plateaus below 2 Å and stops rising. **The binder does not migrate from the peptide onto the helices**, which is the specific failure this run was built to detect and which would appear as peptide contacts decaying while MHC contacts persist. Peptide contacts never approach zero across 663 frames — the minimum is 21 — and the centre-of-mass separation drifts 0.45 Å over the whole trajectory.

This rules out fast failure. It does not rule out slow failure, and a single trajectory on a single variant says nothing about mutant-versus-wild-type discrimination, which is addressed in §10.

### Replica simulations on both variants

Whether the rigidity difference of §2 survives into the bound complex was tested directly: three independent seeds on the G12D complex and three on the wild-type complex, identical binder, identical starting coordinates, the wild type built by deleting the Asp6 side chain from the same file so that every other atom is unchanged. 2.15 ns each. RMSF measured on peptide Cα atoms after superposing each frame on the α1-α2 groove, so it reports peptide motion relative to the groove rather than tumbling of the assembly.

| | G12D | wild type | ratio | |
|---|---|---|---|---|
| anchor RMSF (p1-p3, p9-p10) | 0.562 ± 0.029 Å | 0.570 ± 0.062 Å | 1.01× | overlapping |
| bulge RMSF (p4-p7) | 0.658 ± 0.039 Å | 0.994 ± 0.257 Å | 1.51× | no overlap |
| **position 6 RMSF** | **0.532 ± 0.030 Å** | **1.412 ± 0.264 Å** | **2.66×** | no overlap |
| bulge / anchor | 1.172 ± 0.054 | 1.725 ± 0.257 | 1.47× | no overlap |

**The anchors are the control and they pass**: identical between variants to within 1%, which rules out a difference in global tumbling or sampling quality. Against that baseline every wild-type replica exceeds every G12D replica on bulge RMSF, on position-6 RMSF, and on the bulge-to-anchor ratio. With three replicas per variant, complete separation is an exact one-sided rank test at p = 0.05; the direction was specified in advance by the free-state ensemble, which is what licenses the one-sided test.

**The effect is spatially localised on the mutation.** Per-position ratios run 0.97, 0.92, 0.97, 1.00 across p1-p4, rise to 1.27 at p5, peak at **2.66 at p6**, fall to 1.34 and 1.40 at p7-p8, and return to 1.09 by p9-p10. A global sampling artefact would not produce a spike centred on the substituted residue decaying over two residues either side.

**Mechanistically this closes the loop opened in §1.** The aspartate presents 8.8 Å² and cannot be contacted, but it pins its own backbone segment; in wild type that position is glycine, the least conformationally restrained residue. The difference a binder could read on this target is therefore conformational rather than chemical, and it persists in the bound state.

**What this does not show.** It establishes that the peptide differs between variants while the binder is bound. It does not establish that the binder *exploits* that difference — glycine's intrinsic backbone freedom accounts for the observation without reference to the design at all. Separating the two requires the same replica measurement on apo pMHC and a demonstration that the binder amplifies or depends on the gap. At 377 ns/day that is hours of GPU time, and it is the single most informative experiment remaining.

## 9. Selection

**De novo compliance is a gate, not a criterion.** Every sequence was aligned against β2-microglobulin and against HLA-A\*11:01 α1-α2 by Smith-Waterman with BLOSUM62, scored as identity over the query length. **Seven of twenty-five designs returned 85–93% identity to β2-microglobulin** and were excluded outright. These are residual β2m copies: the `not_binding` constraint of §3 stopped them dominating the run, as they had when the socket was left empty, but it did not stop the model producing them. They pass every structural filter — they contact the peptide region, they fold perfectly alone because they are a real protein — so nothing short of an explicit identity check removes them. Every surviving design sits at 3–13% to β2m and 4–12% to the heavy chain, and the figure is recorded per design in the submission.

The remaining eighteen were ordered under an explicit lexicographic key rather than a weighted score, because a composite would allow a strong cofolding number to mask the absence of peptide contact:

1. **Peptide contact** — more than 50 Å² of peptide surface buried in the refolded complex. A design on the helices cannot satisfy the specificity objectives.
2. **Cysteine-free** — the expression system is reducing.
3. **Folds without the target** — monomer pLDDT ≥ 0.80.
4. **On-target interface confidence** — ipSAE.

Where the counter-screen of §6 was run, its measured peptide dependence takes precedence over the ipSAE tie-break: `_09` is ranked first despite a lower ipSAE than `_06`, because it is three times more peptide-dependent and that is the objective being scored.

Ranks 1–9 pass all four. Ranks 10–11 contact the peptide well but do not fold independently. Ranks 12–13 contact the peptide and carry cysteines. Ranks 14–18 make no peptide contact and are declared as helix-surface controls.

**The epitope declared for each design is measured, not assumed** — the peptide residues within 4.5 Å of the binder in that design's own refolded complex, listed per design in the submission.

## 10. Direction

**Screening breadth.** The counter-screen and cross-panel of §6 were run on the leading designs; extending both across the whole submission turns a two-design comparison into a selectivity profile for the set. The panel itself widens naturally to a larger IEDB-derived collection of A\*11:01 ligands, and to the G12D decamer presented on other common alleles — the allele-level selectivity question any therapeutic eventually faces, and one this design's rejection of unrelated A\*11:01 ligands makes worth asking.

**Ensemble molecular dynamics.** The single trajectory in §8 establishes that the complex is stable; it cannot address discrimination. The rigidity difference in §2 is the accessible discriminating signal, and separating it from sampling noise requires independent replicas on both the mutant and wild-type complexes, with the bulge RMSF difference compared against between-replica variance. At the 377 ns/day measured here that is a day of GPU time, not a research programme.

**Semiempirical QM electrostatics.** A buried carboxylate in a low-dielectric groove is the regime where fixed-charge force fields are least reliable — no polarisation, no charge transfer, and point charges that cannot represent the anisotropy of a carboxylate field at close range. GFN2-xTB carries an explicit anisotropic electrostatics term, and the quantity of interest, ΔΔE between variants on atom-for-atom identical truncations, cancels truncation and basis-set superposition error to first order. Precedent exists for semiempirical QM on pMHC single-substitution panels.

**A DARPin arm.** Designed ankyrin repeat proteins are the scaffold best matched to this expression system: no cysteines, no disulfides, very high thermal stability, and reliable bacterial and cell-free expression — the failure mode that kills de novo designs for reasons unrelated to their interface. The scaffold is precedented against this exact target class, with DARPin binders to peptide-MHC complexes and DARPin T-cell engagers directed at MHC-complexed antigens developed out of the Plückthun laboratory and Molecular Partners. Its elongated concave paratope is a plausible geometric match for a peptide lying in a groove between two helices, which is the shape problem a compact miniprotein has to solve by arcing over the top.

A rigid-body feasibility scan was run to test whether the geometry works at all, since the elongated concave paratope of a repeat protein is a very different shape from a compact miniprotein arcing over a groove. The paratope was derived from the structure rather than assumed — the three internal repeats of E3_5 (PDB 1MJ0) located by their conserved TPLHLAA motif, aligned, and the positions differing between them taken as the randomised library positions, recovering 21 residues across the three repeats. The scaffold was then swept through spin, tilt and lateral offset and pushed to first contact against the pMHC target, scored by peptide contacts and buried peptide surface. **The identical scan was run with a submitted de novo binder as a control.**

The control is what makes the result reportable, and what bounds it. It recovers 82 Å² of buried peptide surface for a binder whose actual pose buries 194 Å² — 42% of the truth. Rigid-body search without side-chain or backbone relaxation therefore cannot measure a ceiling here, and the DARPin's 77 Å² is a floor set by the search, not a limit set by the geometry. What the comparison does support is that **a DARPin is not sterically excluded from this groove**: under the identical handicap it reaches 94% of what a known binder reaches. One further observation survives the control's failure because it is a ratio rather than an absolute — the DARPin makes 32 peptide contacts to the control's 15 while burying slightly less surface, the signature of a shallow extended paratope against a compact one.

Two constraints make it a deliberate choice rather than a default. The consensus framework is fixed and published, so a DARPin is de novo only at its randomised interface positions — a submission would need the compliance argument made on paratope identity rather than whole-sequence identity. And an N2C DARPin is roughly 96 residues against the 55–90 window adopted here from Adaptyv's own expression data. Both are worth accepting in a campaign with time to screen the scaffold against the miniprotein arm head to head.

**Surface divergence analysis.** Persistent homology and Fused Gromov-Wasserstein divergence over the free-state ensembles, to characterise the mutant-versus-wild-type surface difference as a distribution rather than a pairwise structural comparison.
