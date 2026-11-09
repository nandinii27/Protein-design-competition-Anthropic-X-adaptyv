# Targets

USE THESE. Built by scripts/build_free_targets.py from g12d_pmhc.pdb.
Chains A = HLA-A*11:01 alpha1-alpha2 (182), B = beta2-microglobulin (99),
C = peptide (10).

    target_g12d.pdb   VVVGADGVGK, the design and on-target analyte
    target_wt.pdb     VVVGAGGVGK, built by side-chain deletion
    target_apo.pdb    empty groove

SUPERSEDED. The design_*.pdb and design_*.cif files are alpha1-alpha2 plus
peptide with beta2-microglobulin removed. Designing against them leaves the
beta2m socket empty and the generative model fills it: 142 designs, zero
peptide contact, 50 of them more than 80 percent identical to beta2m. They
are kept because methods.md section 3 describes that failure and these are
the evidence. Do not design against them.

SOURCE. ternary.pdb, g12d_pmhc.pdb, wt_pmhc.pdb and apo_mhc.pdb come from
7OW6 via scripts/prepare_targets.py.
