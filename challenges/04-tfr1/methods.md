# Challenge 04 — TfR1 conditional binder

## Target

Human transferrin receptor 1 (TfR1 / CD71), soluble ectodomain homodimer,
residues 121–760 (UniProt P02786). Structure reference: PDB 11QC, chains A–B.
Mouse ortholog: UniProt Q62351, construct 123–763.

Three objectives, ranked in this order by the challenge: pH-selective binding
(measurable K_D against human TfR1 at pH 7.4, no detectable binding at pH 6.0),
mouse cross-reactivity at pH 7.4, and affinity.

## Epitope selection

### The constraint the brief does not state

The organisers recommend the apical domain, and separately require mouse
cross-reactivity. These pull in opposite directions. The apical domain is the
receptor's species-divergence hotspot — it is the determinant New World
arenaviruses and parvovirus B19 read to distinguish human from rodent TfR1, and
it is why TfR1-targeted brain shuttles are built twice, once per species.

Quantified on this target, from a global alignment of the two ectodomains
(`src/tfr1_epitope.py`):

| | value |
|---|---|
| Ectodomain identity, human vs mouse | 78.1% |
| Apical residues modelled in 11QC (189–383) | 195 |
| Apical surface residues (SASA ≥ 25 Å²) | 100 |
| — conserved human/mouse | 59 |
| — divergent | 41 |

The divergent surface set contains the species-specificity loop itself:
`R208N, V210-(deletion), Y211D, L212P`, plus `D204Q, K205S, N215S, G217E`.
A binder whose interface lands there is a human-only binder by construction and
loses objective 2 before the assay runs. That result came out of conservation
and solvent exposure alone, with no reference to a viral complex.

So the epitope was not taken as "the apical domain" but as a conserved,
solvent-exposed, spatially contiguous sub-patch of it, under three filters:

1. **Conservation** — identical residue in human and mouse.
2. **Exposure** — residue SASA ≥ 25 Å², computed on the biological *dimer*, so
   A/B interface residues are correctly scored as unavailable. A monomer
   calculation passes several of them.
3. **Contiguity** — greedy patch growth within a 13 Å Cβ radius of the running
   centroid, so the hotspot list describes one face rather than scattered
   positions.

Residues in the transferrin footprint (helical domain 607–760 plus the reported
protease-like contacts) were excluded throughout, so the binder does not compete
with iron delivery and the receptor keeps recycling — the same constraint used
by transferrin-sparing brain shuttles.

### Epitope 1 — conservation-first

`K240, K241, E244, I262, E266, S273`

All six identical in mouse. Net formal charge zero, two lysines and two
glutamates available as electrostatic anchors. The maximum-exposure conserved
patch on the apical face.

### Epitope 2 — pH handle on the target side (used)

`T263, F264, E266, A269, H300, L303`

Centred on **His300**: solvent-exposed (42 Å²), apical, conserved in mouse.
Shares E266 with epitope 1, so the two patches are adjacent rather than
unrelated.

The reason for this patch is mechanistic. The standard route to pH selectivity
is histidine engineering on the *binder* — substituting His into the designed
interface so protonation below its pKa breaks the interaction. Its documented
failure modes: His substitutions usually cost affinity at both pH values rather
than creating a switch; reported switch ratios are typically single-digit while
this assay demands *no detectable binding* at pH 6.0; and a His buried at a
designed interface frequently has its pKa suppressed below 6.0, so it never
titrates inside the assay window.

Designing against a histidine already present on the target inverts that. His300
titrates between pH 7.4 and 6.0 regardless of what the binder does, so the switch
is encoded in the receptor and the binder does not sacrifice its own core to
carry it. Acidic residues cannot serve this role: Asp and Glu have pKa near 4 and
remain deprotonated at 6.0.

## What was actually run

**Intended:** Proteina-Complexa locally on a Colab GPU, hotspot-conditioned on
EP2. The install could not be made to produce designs in the time available; the
fixes are in `proteina_complexa_colab.ipynb` and documented below.

**Submitted:** 10 designs generated through NVIDIA's hosted RFD4-Proteina
endpoint, which needs no local GPU. Unfiltered. No pH engineering applied — the
His300 rationale above shaped hotspot choice only.

This is not a strong submission and is not presented as one. The defensible work
in this folder is the epitope analysis.

## Proteina-Complexa install notes

The public release assumes it is running inside its own Docker image. Five
separate consequences on a fresh host, all patched in the notebook:

1. `env/build_uv_env.sh` passes `--system` to every `uv pip install`, so uv
   ignores the venv it just created and installs into the host Python. Correct
   as root in a disposable container; wrong anywhere else. On Colab that means
   Python 3.13.
2. Downstream of (1): `dm-tree` 0.1.8 — hard-pinned by `proteinfoundation` —
   publishes no cp313 wheel, so it builds from source, and its bundled abseil
   (lts_20210324) does not compile against GCC 13.
3. `atomworks` is installed unversioned behind `|| echo "Warning: ..."`, so a
   wrong version installs silently. Current PyPI is 3.0.4, which removed
   `AtomSelectionStack` from `atomworks.io.utils.selection`; `pdb_utils.py`
   imports it, so `generate` dies on an ImportError. **Pin 2.2.1** — verified
   present there, absent from 3.0.0 onward.
4. `.env` ships with `LOCAL_CODE_PATH=/path/to/Proteina-Complexa` as a literal
   placeholder, and most other paths derive from it.
5. `complexa init` is two-phase and does not say so: `init` writes `.env`,
   `init uv` writes `env.sh`. Nothing works until both have run.

Also worth knowing: `complexa validate` rejects Hydra `++` overrides while
`complexa generate` accepts them, and `--verbose` does not actually suppress
file logging — run output still lands in `logs/`.

## Known limitations

**pH selectivity cannot be screened in silico here.** Neither AlphaFold-class
co-folding nor the design model takes protonation state as input. An ipTM or
ipSAE value computed "at pH 6.0" measures nothing. The pH strategy is structural
— place the interface over a titratable target-side histidine — and must be
validated by assay.

**Single-target generation.** RFD4-Proteina conditions on one target, so mouse
cross-reactivity is a post-hoc property rather than an optimisation objective.
EP2 makes this survivable only because all six hotspots are identical in human
and mouse.

**Threaded mouse model.** `mouse_TfR1_apical.pdb` is the mouse sequence threaded
onto the human 11QC backbone, with the 48 substituted positions in 189–383
truncated to Cβ so a design model rebuilds them rather than inheriting human
rotamers. Backbone shifts in divergent loops are not captured, which
systematically favours conserved-core epitopes. Read any mouse-side score as
optimistic.

## Reading

- Kariolis et al. 2020, *Sci Transl Med* — transferrin-sparing apical engagement
  for brain delivery.
- Yu et al. 2011, *Sci Transl Med* — lower anti-TfR1 affinity gives better brain
  uptake; high affinity drives lysosomal sorting.
- Lawrence et al. 1999, *Science* — ectodomain architecture and the domain
  boundaries used here.
- Radoshitzky et al. 2008, *PNAS* — apical-domain determinants of TfR1 species
  specificity.
- Igawa et al. 2010, *Nat Biotechnol* — histidine engineering for pH-dependent,
  recycling binders.
