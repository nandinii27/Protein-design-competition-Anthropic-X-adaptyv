# Challenge 05 — de novo binders to HLA-A\*11:01 / KRAS G12D

Anthropic × Adaptyv 2026, final challenge. Design a protein that binds the KRAS G12D neoantigen decamer `VVVGADGVGK` presented on HLA-A\*11:01, and that discriminates it from the wild-type peptide and from the empty groove.

Full write-up in [`methods.md`](methods.md).

---

## The finding

**The mutation is not available for contact.** Per-residue solvent accessibility on the peptide in 7OW6 puts Asp6 — the residue that is Gly in wild-type KRAS — at **8.8 Å²**, with one carboxylate oxygen entirely buried. The therapeutic TCR in that same structure buries **zero** of it.

So the standard protocol for a neoantigen target, placing diffusion hotspots on the mutated residue, is wrong here. Discrimination has to come from indirect readout: the conformation the buried charge imposes on the peptide backbone, or the difference in dynamics between the variants.

**Four independent methods agree that direct readout was not achieved**, which is the result rather than a shortfall:

| method | observable | outcome |
|---|---|---|
| solvent accessibility | Å² exposed at p6 | 8.8 Å², one oxygen fully buried |
| ProteinMPNN with apo control | Δ log-odds P(Asp)/P(Gly) | −0.35 to −1.40 for every clean design |
| three-target counter-screen | ipTM, G12D vs WT | +0.005 (lead design), within noise |
| GFN2-xTB with ALPB | ΔΔE | −1.10 kcal/mol, within method error |

**What the lead design does achieve** is a quantified selectivity profile:

- **peptide-dependent** — removing the peptide costs 62% of its binder-to-MHC interface confidence (0.827 → 0.317)
- **register-specific** — the 9-mer form of the same neoantigen is rejected at 12%
- **HLA-ligand selective** — two unrelated A\*11:01 epitopes score at 8–18%
- **peaks on its own target** within the KRAS codon-12 family; the runner-up does not, scoring 172% on G12C
- **does not resolve** the small codon-12 substitutions G12A/C/V/S, which sit at 67–85%
- **the rigidity difference is real and localised** — replica MD on both variants (3 seeds each) puts peptide position 6 at 2.66× more mobile in wild type than in G12D, with anchor residues identical between variants and no overlap between replica sets
- **stable in simulation** — 6.63 ns of implicit-solvent MD leaves the binder fold intact (Cα RMSD 1.82 ± 0.17 Å, flat) and the peptide interface unchanged (52.9 ± 10.9 contacts, no decay, minimum 21 across 663 frames)

Accurately: a KRAS codon-12 neoepitope binder with correct register and correct ligand selectivity, which does not resolve the identity of residue 6 among small substitutions. That is the honest ceiling of indirect readout on this target.

---

## Pipeline

Each stage writes a CSV that the next one reads. Run in order.

**1 — Target preparation.** `prepare_targets.py`, `build_free_targets.py`
Pulls 7OW6, builds the G12D, wild-type and apo targets, computes per-residue peptide accessibility. The wild type is built by deleting the Asp side-chain atoms, so the backbone geometry is preserved exactly rather than re-predicted.

**2 — Free-state ensemble.** `freestate_ensemble.py`, `compare_ensembles.py`
Five-model AlphaFold2 ensembles of the free G12D and wild-type trimers, superposed on the groove Cα atoms only and measured in place. Establishes that induced-fit remodelling is localised (2–3 Å at the bulge, 0.3–0.9 Å at the anchors) and that the wild-type bulge is 2.15× more conformationally variable.

**3 — Design.** `pmhc_g12d.yaml` → BoltzGen
Full trimer as target, binding directed to peptide positions 4–8, and `not_binding: all` on β2-microglobulin. That last constraint is load-bearing: a trimmed α1-α2 target leaves the β2m socket empty and the model fills it — 142 designs, zero peptide contact, 50 of them >80% identical to β2m.

**4 — Geometry filter.** `geometric_filter.py`
Buried surface split between peptide and MHC, hotspot coverage, cysteine count, length. Identifies chains by content rather than letter, so it runs on any output layout.

**5 — Confidence.** `extract_fold.py`
Recovers ipSAE, ipTM and interaction PAE from BoltzGen's folding-step `.npz` outputs.

**6 — Monomer self-consistency.** `monomer_refold.py`
Folds each binder alone with Chai-1 single-sequence. Not redundant with ipSAE — it demoted two designs ranked 3rd and 5th on confidence whose monomer pLDDT is 0.665 and 0.725.

**7 — Mutation readout.** `mpnn_specificity_test.ipynb`
ProteinMPNN conditional probabilities at the masked peptide, run twice per design: once as the complex, once with the binder deleted. The apo control is the experiment — the groove already prefers Asp at p6, so only the shift attributable to the binder means anything.

**8 — Counter-screen.** `counterscreen3.py`
Each design against the G12D trimer, the wild-type trimer and the empty groove, within one predictor, MSA-backed, with a seed replicate on the leading design.

**9 — Cross-panel.** `crosspan.py`
Ten A\*11:01 ligands: five further codon-12 substitutions, G13D, the 9-mer register, and two unrelated EBV epitopes.

**10 — Quantum electrostatics.** `xtb_ddE.py`
GFN2-xTB ΔΔE on a 5 Å truncation with ALPB solvation and link-atom caps. The free-binder term cancels exactly; truncation and basis-set superposition error cancel to first order because both variants use atom-for-atom identical truncations.

**11 — Molecular dynamics.** `md_gpu.py`
Implicit-solvent MD of the complex, reporting binder RMSD, peptide contacts, MHC contacts and centre-of-mass separation. Benchmarks itself before committing and aborts below 4 ns/day. 6.63 ns on the lead design at 377 ns/day on a T4 via OpenMM's OpenCL platform; the interface holds and the binder does not migrate onto the helices.

**Scaffold feasibility.** `darpin_feasibility.py`
Rigid-body scan asking whether a DARPin paratope can reach a peptide in a groove between two helices, with a submitted binder as control. The control recovers only 42% of its own known pose, which bounds what the result can claim — reported as "not sterically excluded", not as a ceiling.

**12 — Replica MD on both variants.** `md_replicas.py`
Three seeds each on the G12D and wild-type complexes, same binder, wild type built by side-chain deletion from the identical file. Reports per-residue peptide RMSF with the anchors as an internal control.

**Figures.** `render_design.py`, `animate_design.py`
Matched-camera stills, and a camera flight over the complex breathing on its lowest-frequency normal mode.

---

## Results

| file | contents |
|---|---|
| `submission.csv` | 20 designs, ordered, each with a measured epitope |
| `master_metrics.csv` | every measurement, per design |
| `ontarget_metrics.csv` | Boltz-2 confidence and interface geometry |
| `monomer_refold.csv` | target-free fold quality |
| `mpnn_specificity.csv` | mutation readout with apo controls |
| `counterscreen3.csv` | G12D / wild type / apo |
| `crosspan.csv` | ten-ligand selectivity panel |
| `xtb_ddE.json` | GFN2-xTB discrimination energy |
| `md_metrics_09.csv` | 6.63 ns trajectory metrics, 10 ps resolution |
| `replicas.csv` | per-residue peptide RMSF, 3 replicas per variant |

## Selection

Twenty designs chosen under an explicit lexicographic key rather than a weighted score, because a composite would let a strong confidence number mask the absence of peptide contact:

1. contacts the peptide at all — more than 50 Å² buried
2. cysteine-free, for a reducing cell-free expression system
3. folds without the target — monomer pLDDT ≥ 0.80
4. on-target interface confidence (ipSAE)

Measured peptide dependence from the counter-screen overrides the confidence tie-break at the top, which is why the lead design is ranked above one with a higher ipSAE.

## Environment

Designs generated on Kaggle T4. Analysis, QM and figures on CPU. OpenMM's CUDA platform is absent from the PyPI wheel; its OpenCL platform runs on the T4 and agrees with the reference implementation to 7×10⁻⁶.
