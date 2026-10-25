"""
Target definitions and constants for the GLP1R agonist design challenge.

All sequences here are from public databases. Verify before submission:
  - GLP1R:     UniProt P43220 (isoform P43220-1)
  - GLP-1:     UniProt P01275 (proglucagon, GLP-1(7-36)amide fragment)
  - Exendin-4: UniProt P26349 (Heloderma suspectum exendin-4)
"""

# ---------------------------------------------------------------------------
# Target receptor
# ---------------------------------------------------------------------------

# Human GLP1R, UniProt P43220-1, residues 24-432 (the assay construct range).
GLP1R_24_432 = (
    "RPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWA"
    "SSVPQGHVYRFCTAEGLWLQKDNSSLPWRDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVI"
    "ASAILLGFRHLHCTRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVF"
    "LLMQYCVAANYYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTR"
    "NSNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLLGTHEVIFA"
    "FVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQLEFRKSWERWRLEHLHIQRDSS"
)

# Offset of the first residue of GLP1R_24_432 in UniProt numbering.
UNIPROT_OFFSET = 24

# Extracellular domain. The ECD of class B1 receptors runs from the mature
# N-terminus to roughly the start of TM1. Residues 24-145 (UniProt numbering)
# cover the three-layer alpha-beta-beta-alpha fold plus the three conserved
# disulfides. This is the region to condition sequence-only generators on,
# because the rest of the construct is buried in detergent or membrane and is
# not a bindable surface.
ECD_START_UNIPROT = 24
ECD_END_UNIPROT = 145


def uniprot_slice(start: int, end: int) -> str:
    """Return the construct subsequence for an inclusive UniProt residue range."""
    if start < UNIPROT_OFFSET:
        raise ValueError(f"start {start} is before the construct start {UNIPROT_OFFSET}")
    i = start - UNIPROT_OFFSET
    j = end - UNIPROT_OFFSET + 1
    if j > len(GLP1R_24_432):
        raise ValueError(f"end {end} is past the construct end")
    return GLP1R_24_432[i:j]


GLP1R_ECD = uniprot_slice(ECD_START_UNIPROT, ECD_END_UNIPROT)

# ---------------------------------------------------------------------------
# Structural reference
# ---------------------------------------------------------------------------

# Active-state cryo-EM structure: human GLP1R bound to GLP-1, Gs-coupled.
REFERENCE_PDB_ID = "6X18"
RECEPTOR_CHAIN = "R"

# Chain holding the bound GLP-1 peptide in the reference structure. Confirm this
# against the actual mmCIF before running the geometry module; chain labelling
# varies between depositions.
REFERENCE_PEPTIDE_CHAIN = "P"

# ---------------------------------------------------------------------------
# Functional control set
# ---------------------------------------------------------------------------
# These are the ground truth for whether a scoring stack can separate agonism
# from binding. Agonists carry an intact N-terminus that inserts into the TMD
# core; antagonists are N-terminally truncated, bind the ECD, and do not
# activate. Any filter that ranks the antagonists above the agonists selects
# for the wrong objective.

CONTROL_PEPTIDES = {
    "GLP-1(7-36)": {
        "sequence": "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
        "label": "agonist",
        "note": "endogenous hormone, amidated C-terminus in vivo",
    },
    "exendin-4": {
        "sequence": "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
        "label": "agonist",
        "note": "full agonist, basis of exenatide",
    },
    "exendin(9-39)": {
        "sequence": "DLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
        "label": "antagonist",
        "note": "N-terminally truncated exendin-4, competitive antagonist",
    },
    "GLP-1(9-36)": {
        "sequence": "EGTFTSDVSSYLEGQAAKEFIAWLVKGR",
        "label": "antagonist",
        "note": "DPP-4 cleavage metabolite, essentially inactive at the receptor",
    },
}

# ---------------------------------------------------------------------------
# Submission constraints (from the challenge page)
# ---------------------------------------------------------------------------

MIN_LENGTH = 10
MAX_LENGTH = 250

# Design-class stratification used by the organisers.
MINIBINDER_RANGE = (40, 100)

# Length windows to sample for the agonist-architecture tracks. The natural
# agonists are 30-39 residues; a continuous helix spanning the ECD groove and
# reaching the TMD pocket needs roughly that span.
AGONIST_LENGTHS = (28, 32, 36, 40)

# Detergent-specific filter. The receptor is supplied in DDM/CHS micelles, so
# designs with large exposed hydrophobic surface partition into detergent and
# give aggregation rather than signal. GRAVY above this is a reject.
MAX_GRAVY = 0.0

# Free cysteines complicate expression and can scramble against the three
# native ECD disulfides.
MAX_FREE_CYSTEINES = 0
