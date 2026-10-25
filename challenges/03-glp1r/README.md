# Machine Learning for Rare Disease Research — MVA Hackathon 2026

Analysis for **Rare Disease, Real Kid: The MVA Hackathon 2026**, organised by
Sage Bionetworks with the MVA Society, Hugging Face and BEACON.

The case is a child with **Mosaic Variegated Aneuploidy (MVA)**, an ultra-rare
autosomal recessive chromosomal instability syndrome affecting fewer than 50
people worldwide. Cells across the body carry different chromosome counts
because the mitotic spindle assembly checkpoint is defective.

---

## Results

### Track 1 — Variant prediction

**Submitted. Scored 100 / 100 rank points, F-max 1.000, full match at rank 1.**

Causal genotype: compound heterozygous *BUB1B* (GRCh38).

| | Variant | HGVSc | HGVSp | Consequence | gnomAD AF |
|---|---|---|---|---|---|
| Allele 1 | chr15:40209701 T>G | c.2210T>G | p.Leu737Ter | stop gained, exon 17/23 | 1.0e-4 |
| Allele 2 | chr15:40220612 T>G | c.3006T>G | p.Asn1002Lys | missense, exon 23/23 | absent |

A null allele predicted to undergo nonsense-mediated decay, in trans with a
missense allele in the C-terminal pseudokinase domain.

### Track 2 — Drug repurposing

**Not submitted as an official entry.** The analysis is published here as
`track2_report.pdf` for anyone who wants to read the approach.

Its argument in one line: the lesion degrades a scaffolding interface rather
than an enzyme, which makes the tractable intervention a **rebalancing of the
Aurora B / PP2A-B56 set point** rather than cytotoxic clearance of aneuploid
cells — and the standard aneuploidy-exploitation strategies are contraindicated
here because the patient's own tissue shares the vulnerability they depend on.

---

## Contents

| File | What it is |
|---|---|
| `track2_report.pdf` | Full Track 2 analysis: mechanism, proposal, rejected alternatives, limitations |
| `voi_analysis.py` | Value-of-information model for choosing the next experiment |
| `mosaic_aneuploidy_burden.py` | Mosaic aneuploidy screen from the call set (BAF deviation + read depth) |
| `track1_report_v1.md` | Track 1 methods report as submitted |

---

## Method, Track 1

Pipeline, in order:

1. **Gene region definition.** Canonical GRCh38 coordinates for *BUB1B*,
   *CEP57* and *TRIP13* retrieved from the Ensembl REST lookup endpoint rather
   than hardcoded, padded 20 kb per flank to retain promoter and deep intronic
   territory. Gene body and flank scored separately — an unpartitioned region
   yields a misleadingly dense candidate set dominated by neighbouring haplotype
   structure.
2. **Panel extraction.** Single streaming pass over the VCF, retaining records
   inside the padded intervals. 105 variants; only 14 of the *BUB1B* records lie
   within the transcribed gene.
3. **Genotype and quality parsing.** GT, DP, GQ, PID, FILTER. Note the INFO `AF`
   field here is the within-sample allele fraction from the caller, **not** a
   population frequency, and was not used for rarity filtering. Non-PASS records
   downweighted rather than discarded.
4. **Annotation.** Ensembl VEP REST endpoint: canonical transcript consequence,
   HGVS, exon/intron numbering, SIFT, PolyPhen, gnomAD frequencies.
5. **Ranking.** Conjunction of population rarity, predicted functional
   consequence, and compatibility with biallelic recessive inheritance in a
   single gene. No model was fitted; at N=1 fitting one would not be meaningful.

**Exclusions.** *CEP57* carried no novel or rare variants. *TRIP13*'s two
gene-body candidates share a read-backed phasing identifier and are therefore in
cis, leaving one intact copy and excluding a recessive mechanism at that locus.

---

## Supporting analyses

### `voi_analysis.py`

Prioritisation pipelines return a ranked list but do not say what would resolve
the remaining ambiguity. This scores each candidate next experiment by expected
entropy reduction over a hypothesis set, reported in bits and normalised by cost
and turnaround.

Two non-obvious outputs: parental sequencing is less decisive than intuition
suggests, and a functional checkpoint assay has near-zero *diagnostic* value
while being the highest-value class for the *mechanistic* question. An
experiment's worth is defined against a specified uncertainty, not in the
abstract.

All priors are explicit elicited estimates, stated inline so they can be
challenged and substituted. They are not fitted and not claimed to be
calibrated.

```powershell
python voi_analysis.py
```

### `mosaic_aneuploidy_burden.py`

Estimates mosaic aneuploidy burden per chromosome directly from the call set,
using B-allele frequency deviation at heterozygous SNPs plus relative read
depth.

**Result: suggestive, not established.** Sixteen chromosomes cluster tightly at
the binomial noise floor; four are outliers implying 2-8% nominal mosaic
fraction. But the outlier ranking correlates with GC-rich chromosomes, the
acrocentrics and chr9 heterochromatin, which is characteristic of mappability
artefact rather than biology. The depth arm is unreliable at this resolution
because median depth is an integer. A validated mosaic chromosomal alteration
caller with GC correction would be needed to separate signal from artefact.

Reported as a screen, not a diagnostic.

```powershell
python mosaic_aneuploidy_burden.py
```

---

## Running the code

```powershell
pip install pandas requests matplotlib
```

The scripts expect the challenge VCF at `./data/WGS_EX2312012_HGWCNDSX7.vcf.gz`.
Access is gated; request it on the
[dataset page](https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data)
and download with your own token.

---

## Data handling

**This repository contains no patient sequence data or genotype-scale
derivative, and never will.** `data/` is gitignored.

Under the hackathon data transfer agreement the proband's VCF, phenotype
document and every genotype-scale derivative stay on the analysis machine and
are deleted at the conclusion of the hackathon. What is published here is code,
configuration and derived findings — the named candidate variants and
chromosome-level aggregate statistics — which the challenge terms explicitly
permit.

Participant submissions, code and reports are released under **CC BY 4.0**.

---

## AI assistance disclosure

Anthropic API, Claude, commercial terms, no training on customer content.

---

## Acknowledgements

This work was made possible through the Hackathon, organized by Sage Bionetworks
in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking,
Evaluation, and Assessment Consortium for Science), with prize sponsorship from
AWS and Anthropic. I am deeply grateful to the child and their family who
generously contributed their data and their story to advance research into this
rare disease.
