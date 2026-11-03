"""Author an RFD4-Proteina annotated mmCIF for the TfR1 apical-domain binder task.

Format is reproduced from NVIDIA's own de_novo_binder.cif example:
  entity 0 = target, chain B, coordinates and sequence fixed, hotspots flagged per atom
  entity 1 = designed binder, chain A, ONE placeholder HETATM row whose label_comp_id
             is the length range, with expsegmin/expsegmax carrying min/max length

The _chem_comp and _chem_comp_bond blocks are copied verbatim from the example,
which already defines all 20 standard residues.
"""
import re, sys, json
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1
import warnings
warnings.filterwarnings("ignore")

EXAMPLE = Path(sys.argv[1])
TARGET_PDB = Path(sys.argv[2])
OUT_CIF = Path(sys.argv[3])

BINDER_MIN, BINDER_MAX = 55, 95
HOTSPOT_ATOMS = {206: ["OD1", "ND2"], 209: ["CD1", "CD2"], 374: ["CE", "NZ"]}

T2O = {k.upper(): v for k, v in protein_letters_3to1.items()}
ex = EXAMPLE.read_text().splitlines()

# ---- lift the _chem_comp and _chem_comp_bond loops verbatim ----
def block_from(header_last_tag):
    i = next(k for k, l in enumerate(ex) if l.strip() == header_last_tag)
    # walk back to the "loop_"
    j = i
    while ex[j].strip() != "loop_":
        j -= 1
    out = [ex[j]]
    k = j + 1
    while ex[k].strip().startswith("_"):
        out.append(ex[k]); k += 1
    while k < len(ex) and not ex[k].startswith("loop_") and ex[k].strip() != "#":
        out.append(ex[k]); k += 1
    return out

chem_comp = block_from("_chem_comp.formula_weight")
chem_bond = block_from("_chem_comp_bond.pdbx_aromatic_flag")

# ---- read the target ----
st = PDBParser(QUIET=True).get_structure("t", str(TARGET_PDB))[0]["A"]
residues = [r for r in st if r.id[0] == " " and r.get_resname() in T2O]
seq = "".join(T2O[r.get_resname()] for r in residues)
assert len(residues) == 195, len(residues)

def wrap(s, n=80):
    return "\n".join(s[i:i+n] for i in range(0, len(s), n))

L = []
L.append("data_tfr1_apical.cif")
L.append("#")
L.append("_entry.id       tfr1_apical.cif")
L.append("_entry.author   rfd4")
L.append("_entry.date     2026-11-02")
L.append("#")
L.append("loop_")
L.append("_entity.id ")
L.append("_entity.type ")
L.append("_entity.pdbx_description ")
L.append("0 polymer 'Transferrin receptor 1 apical domain (P02786 189-383)'")
L.append("1 polymer 'de novo binder'")
L.append("#")
L.append("loop_")
L.append("_entity_poly.entity_id ")
L.append("_entity_poly.type ")
L.append("_entity_poly.nstd_linkage ")
L.append("_entity_poly.nstd_monomer ")
L.append("_entity_poly.pdbx_seq_one_letter_code ")
L.append("_entity_poly.pdbx_seq_one_letter_code_can ")
L.append("_entity_poly.pdbx_strand_id ")
L.append("_entity_poly.pdbx_target_identifier ")
L.append("0 polypeptide(L) no no  ")
L.append(";" + wrap(seq))
L.append(";")
L.append(" ")
L.append(";" + wrap(seq))
L.append(";")
L.append(" B ?")
L.append("1 polypeptide(L) no yes X                                                                                                                         X                                                                                                                         A ?")
L.append("#")
L.append("loop_")
L.append("_entity_poly_seq.entity_id ")
L.append("_entity_poly_seq.hetero ")
L.append("_entity_poly_seq.mon_id ")
L.append("_entity_poly_seq.num ")
for i, r in enumerate(residues, start=1):
    L.append("0 no  %s %d" % (r.get_resname(), i))
L.append("1 no  %d-%d %d" % (BINDER_MIN, BINDER_MAX, 1))
L.append("#")
L += chem_comp
L.append("#")
L += chem_bond
L.append("#")

# ---- atom_site ----
COLS = ["group_PDB","type_symbol","label_atom_id","label_alt_id","label_comp_id",
        "label_asym_id","label_entity_id","label_seq_id","pdbx_PDB_ins_code",
        "auth_seq_id","auth_comp_id","auth_asym_id","auth_atom_id","B_iso_or_equiv",
        "occupancy","pdbx_formal_charge","annotation_expsegmax_1_residue",
        "mask_index_1_residue","mask_expsegmax_1_residue","mask_hotspot_1_atom",
        "mask_sequence_1_residue","condition_hbond_1_atom",
        "condition_index-group_1_token","annotation_expsegmin_1_residue",
        "mask_expsegmin_1_residue","mask_coordinate_1_atom",
        "annotation_src-component_1_atom","atom_array_index","Cartn_x","Cartn_y",
        "Cartn_z","pdbx_PDB_model_num","id"]
L.append("loop_")
for c in COLS:
    L.append("_atom_site." + c + " ")

rows = []
# binder placeholder row (entity 1, chain A) -- mirrors the example exactly
rng = "%d-%d" % (BINDER_MIN, BINDER_MAX)
rows.append(["HETATM","X","UNK",".",rng,"A","1","1",".","1",rng,"A","UNK","nan","0.0","?",
             str(BINDER_MAX),"True","True","False","False","0","1",str(BINDER_MIN),
             "True","False",".","0","nan","nan","nan","1","1"])

idx = 1
for seq_i, r in enumerate(residues, start=1):
    hot = HOTSPOT_ATOMS.get(r.id[1], [])
    for a in r:
        if a.element == "H":
            continue
        x, y, z = a.coord
        rows.append([
            "ATOM", a.element or a.get_id()[0], a.get_id(), ".", r.get_resname(),
            "B", "0", str(seq_i), ".", str(seq_i), r.get_resname(), "B", a.get_id(),
            "%.2f" % a.get_bfactor(), "1.0", "?",
            "-1",                                  # expsegmax
            "True",                                # mask_index
            "False",                               # mask_expsegmax
            "True" if a.get_id() in hot else "False",   # mask_hotspot
            "True",                                # mask_sequence
            "0", "0",
            "-1",                                  # expsegmin
            "False",
            "True",                                # mask_coordinate: target is fixed
            "A1", str(idx),
            "%.3f" % x, "%.3f" % y, "%.3f" % z, "1", str(idx + 1),
        ])
        idx += 1

for r in rows:
    assert len(r) == len(COLS), (len(r), len(COLS))
    L.append(" ".join(r))
L.append("#")

OUT_CIF.write_text("\n".join(L) + "\n")
n_hot = sum(1 for r in rows if r[19] == "True")
print("wrote", OUT_CIF, "| atoms:", len(rows) - 1, "| hotspot atoms:", n_hot,
      "| target residues:", len(residues), "| binder range:", rng)
assert n_hot == 6, "expected 6 hotspot atoms, got %d" % n_hot
