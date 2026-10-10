"""
Identify pH-switch anchor residues on the TNF-alpha trimer interface.

Finds titratable residues (HIS, ASP, GLU) on the target surface at and around
the inter-protomer groove, scores them by burial and human/mouse conservation,
and emits hotspot strings for downstream binder design.

Rationale: titratable groups already present on the target contribute to the
binding-linked proton release without costing any stability in the designed
binder. Engaging them with complementary partners is cheaper than grafting
histidines into the binder.

Usage:
    python target_switch_anchors.py --out anchors.csv
    python target_switch_anchors.py --pdb 1tnf.pdb --crop cropped_target.pdb

Dependencies:
    pip install biopython requests
"""

import argparse
import io
import sys
import warnings
from collections import defaultdict

import requests
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.Align import PairwiseAligner

warnings.filterwarnings("ignore")

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

HUMAN_TNF_ACC = "P01375"
MOUSE_TNF_ACC = "P06804"

TITRATABLE = {"HIS", "ASP", "GLU"}

# For a pH switch what must survive across species is the titratable chemistry,
# not the residue identity. Asp and Glu are interchangeable carboxylates with
# near-identical intrinsic pKa; histidine has no substitute in this window.
TITRATABLE_CLASS = {"D": "carboxylate", "E": "carboxylate", "H": "imidazole"}

# Approximate maximum solvent accessible surface area per residue (A^2),
# used to convert absolute SASA into a relative burial fraction.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


def fetch_text(url):
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return resp.text


def load_structure(pdb_id, pdb_path):
    if pdb_path:
        if pdb_path.lower().endswith((".cif", ".mmcif")):
            return MMCIFParser(QUIET=True).get_structure("target", pdb_path)
        return PDBParser(QUIET=True).get_structure("target", pdb_path)
    text = fetch_text(RCSB_PDB_URL.format(pdb_id=pdb_id.upper()))
    return PDBParser(QUIET=True).get_structure("target", io.StringIO(text))


def read_fasta(text):
    return "".join(line.strip() for line in text.splitlines() if not line.startswith(">"))


def get_sequence(acc, local_path):
    """Read a sequence from a local FASTA if given, otherwise fetch from UniProt."""
    if local_path:
        with open(local_path) as fh:
            return read_fasta(fh.read())
    return read_fasta(fetch_text(UNIPROT_FASTA_URL.format(acc=acc)))


def chain_sequence(chain):
    """Return the one-letter sequence and the parallel list of residue objects."""
    seq = []
    residues = []
    for res in chain:
        if res.id[0] != " ":
            continue
        name = res.get_resname()
        if name not in MAX_ASA:
            continue
        seq.append(protein_letters_3to1.get(name, "X"))
        residues.append(res)
    return "".join(seq), residues


def interface_residues(model, chain_a, chain_b, cutoff):
    """Residues of chain_a with any heavy atom within cutoff of chain_b."""
    atoms_b = [
        atom for res in model[chain_b] if res.id[0] == " "
        for atom in res if atom.element != "H"
    ]
    contacts = set()
    for res in model[chain_a]:
        if res.id[0] != " ":
            continue
        for atom in res:
            if atom.element == "H":
                continue
            for other in atoms_b:
                if (atom - other) <= cutoff:
                    contacts.add(res.id[1])
                    break
            else:
                continue
            break
    return contacts


def compute_relative_sasa(structure):
    """Relative solvent accessibility per residue, in the full assembly context."""
    sr = ShrakeRupley()
    model = structure[0]
    sr.compute(model, level="R")
    rel = {}
    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            name = res.get_resname()
            if name not in MAX_ASA:
                continue
            rel[(chain.id, res.id[1])] = res.sasa / MAX_ASA[name]
    return rel


def align_and_map(struct_seq, human_seq, mouse_seq):
    """
    Map each structure sequence position onto the human UniProt sequence, then
    onto the aligned mouse residue. Returns {struct_index: (human_aa, mouse_aa)}.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = None
    aligner.match_score = 2
    aligner.mismatch_score = -1
    # Sequences differ in length (structure fragment vs soluble domain vs full
    # precursor), so overhanging ends must not be penalised.
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    struct_to_human = {}
    aln = aligner.align(struct_seq, human_seq)[0]
    for (s_start, s_end), (h_start, h_end) in zip(aln.aligned[0], aln.aligned[1]):
        for offset in range(s_end - s_start):
            struct_to_human[s_start + offset] = h_start + offset

    human_to_mouse = {}
    aln2 = aligner.align(human_seq, mouse_seq)[0]
    for (h_start, h_end), (m_start, m_end) in zip(aln2.aligned[0], aln2.aligned[1]):
        for offset in range(h_end - h_start):
            human_to_mouse[h_start + offset] = m_start + offset

    mapping = {}
    for s_idx, h_idx in struct_to_human.items():
        m_idx = human_to_mouse.get(h_idx)
        human_aa = human_seq[h_idx]
        mouse_aa = mouse_seq[m_idx] if m_idx is not None else None
        mapping[s_idx] = (human_aa, mouse_aa)
    return mapping


class CropSelect(Select):
    def __init__(self, keep):
        self.keep = keep

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        return (chain_id, residue.id[1]) in self.keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-id", default="1TNF")
    ap.add_argument("--pdb", default=None, help="local PDB file instead of fetching")
    ap.add_argument("--chain-a", default="A")
    ap.add_argument("--chain-b", default="B")
    ap.add_argument("--contact-cutoff", type=float, default=5.0)
    ap.add_argument("--shell", type=float, default=16.0,
                    help="radius around the exposed groove to search for anchors")
    ap.add_argument("--min-exposure", type=float, default=0.15,
                    help="minimum relative SASA for a residue to count as reachable")
    ap.add_argument("--out", default="anchors.csv")
    ap.add_argument("--crop", default=None, help="write a cropped target PDB here")
    ap.add_argument("--human-fasta", default=None,
                    help="local human TNF-alpha FASTA instead of fetching UniProt")
    ap.add_argument("--mouse-fasta", default=None,
                    help="local mouse TNF-alpha FASTA instead of fetching UniProt")
    args = ap.parse_args()

    structure = load_structure(args.pdb_id, args.pdb)
    model = structure[0]

    human_seq = get_sequence(HUMAN_TNF_ACC, args.human_fasta)
    mouse_seq = get_sequence(MOUSE_TNF_ACC, args.mouse_fasta)

    rel_sasa = compute_relative_sasa(structure)

    iface_a = interface_residues(model, args.chain_a, args.chain_b, args.contact_cutoff)
    iface_b = interface_residues(model, args.chain_b, args.chain_a, args.contact_cutoff)

    # The A/B contact set includes the large buried protomer-protomer packing
    # interface; only its solvent-exposed part is reachable by a binder. That
    # exposed boundary wraps around the protomer junction rather than forming
    # one compact spot, so distances are measured to the nearest exposed
    # boundary residue rather than to a centroid.
    groove_atoms = []
    for chain_id, ids in ((args.chain_a, iface_a), (args.chain_b, iface_b)):
        for rid in ids:
            res = model[chain_id][(" ", rid, " ")]
            if "CA" not in res:
                continue
            if rel_sasa.get((chain_id, rid), 0.0) < args.min_exposure:
                continue
            groove_atoms.append(res["CA"].coord)
    if not groove_atoms:
        sys.exit("No exposed interface contacts found. Lower --min-exposure.")

    def dist_to_groove(coord):
        return min(float(((coord - g) ** 2).sum() ** 0.5) for g in groove_atoms)

    conservation = {}
    for chain_id in (args.chain_a, args.chain_b):
        seq, residues = chain_sequence(model[chain_id])
        mapping = align_and_map(seq, human_seq, mouse_seq)
        for idx, res in enumerate(residues):
            human_aa, mouse_aa = mapping.get(idx, (None, None))
            conservation[(chain_id, res.id[1])] = (human_aa, mouse_aa)

    rows = []
    keep_for_crop = set()
    for chain_id in (args.chain_a, args.chain_b):
        for res in model[chain_id]:
            if res.id[0] != " ":
                continue
            name = res.get_resname()
            if name not in MAX_ASA:
                continue
            if "CA" not in res:
                continue
            dist = dist_to_groove(res["CA"].coord)
            if dist > args.shell:
                continue
            keep_for_crop.add((chain_id, res.id[1]))
            if name not in TITRATABLE:
                continue
            key = (chain_id, res.id[1])
            # A binder can only engage a residue that is solvent-reachable.
            if rel_sasa.get(key, 0.0) < args.min_exposure:
                continue
            human_aa, mouse_aa = conservation.get(key, (None, None))
            rows.append({
                "chain": chain_id,
                "resnum": res.id[1],
                "resname": name,
                "rel_sasa": round(rel_sasa.get(key, float("nan")), 3),
                "dist_to_groove": round(dist, 2),
                "at_interface": key[1] in (iface_a if chain_id == args.chain_a else iface_b),
                "human_aa": human_aa,
                "mouse_aa": mouse_aa,
                "conserved": (human_aa is not None and human_aa == mouse_aa),
                "switch_conserved": (
                    human_aa is not None and mouse_aa is not None
                    and TITRATABLE_CLASS.get(human_aa) is not None
                    and TITRATABLE_CLASS.get(human_aa) == TITRATABLE_CLASS.get(mouse_aa)
                ),
            })

    # Rank: conserved first, then partially buried (shifted pKa more likely),
    # then proximity to the groove.
    def rank_key(r):
        burial_score = abs(r["rel_sasa"] - 0.35) if r["rel_sasa"] == r["rel_sasa"] else 1.0
        return (not r["switch_conserved"], r["dist_to_groove"], burial_score)

    rows.sort(key=rank_key)

    header = ["chain", "resnum", "resname", "rel_sasa", "dist_to_groove",
              "at_interface", "human_aa", "mouse_aa", "conserved", "switch_conserved"]
    with open(args.out, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in header) + "\n")

    conserved_iface = [r for r in rows if r["switch_conserved"] and r["dist_to_groove"] <= 6.0]
    hotspots = ",".join(f"{r['chain']}{r['resnum']}" for r in conserved_iface)

    print(f"Titratable candidates within {args.shell} A of the groove: {len(rows)}")
    print(f"Switch-conserved and within 6 A of the groove: {len(conserved_iface)}")
    print()
    print("Hotspot string (conserved titratable anchors only):")
    print(hotspots if hotspots else "  none found; widen --shell or relax --contact-cutoff")
    print()
    print("Top anchors:")
    for r in rows[:16]:
        if r["conserved"]:
            flag = "identical"
        elif r["switch_conserved"]:
            flag = f"switch kept ({r['human_aa']}/{r['mouse_aa']})"
        else:
            flag = f"SWITCH LOST ({r['human_aa']}/{r['mouse_aa']})"
        print(f"  {r['chain']}{r['resnum']:>4} {r['resname']}  "
              f"relSASA {r['rel_sasa']:<6} groove {r['dist_to_groove']:>6} A  {flag}")

    if args.crop:
        io_writer = PDBIO()
        io_writer.set_structure(structure)
        io_writer.save(args.crop, CropSelect(keep_for_crop))
        print(f"\nCropped target written: {len(keep_for_crop)} residues")


if __name__ == "__main__":
    main()
