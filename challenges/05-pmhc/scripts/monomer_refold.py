"""
Monomer refold self-consistency.

Folds each designed binder ALONE and asks whether it reaches a confident
fold without the target holding it in place. A design that only folds in
complex is an expression and stability risk, and this check was the known
gap in the challenge-1 submission.

Chai-1 single-sequence mode is appropriate here: a de novo binder of 55-80
residues has no MSA by definition, so nothing is lost by omitting one. That
is not true of the MHC, which is why the complex counter-screen needs a
different route.

Self-contained - sequences are embedded, no other files required.

    python monomer_refold.py /kaggle/working/mono
"""
from __future__ import annotations
import csv, shutil, sys
from pathlib import Path

SEQUENCES = {
    "pmhc_g12d_08": "SNLEEKCEEAGKKFVAGAEEVTITLDDGVTYTFKRVDDNTIVLVSADPPVDFDTADLITSLVTYHYLKAL",
    "pmhc_g12d_17": "VRQKLVDAARAGDLEAVEEAIKEGEEKLSKEELIEALGRALGVADAYGHDEIARLLGDRLDALNAEL",
    "pmhc_g12d_24": "ARVVDAVLADDHLGAVLVEGITFAEAERIAEEVLKEVPGVKFAAAGRLTDNRYVISIIAESKEALEAALAKAVELVEA",
    "pmhc_g12d_07": "MLLYASGTLDLKPGDVITVNVYDEDGKLVGKIDLEIKEDGSVVVISEEGDLTEEEIKAAVEWMEQVTAGLR",
    "pmhc_g12d_05": "MGKWEVYLISGETKTFDTLEEAIEAALEDLAAAGELEKVELIKIDDNIYGVVGWPKGSTNPVGALIVKL",
    "pmhc_g12d_06": "SPLLLEVGRTGGADNDWEKMAALLEKIAALPDEELRALAAEDPNTMLVLAEIAEVLGLPEVAAVLRRAA",
    "pmhc_g12d_13": "MLARVVSGSGTSAVLETEGLKPGEAAVLQSLSTGRVVGAVGVDETTASARGASLSVGEVLRVLAQGPAEVSAAEL",
    "pmhc_g12d_10": "SEELEEALRELGKKLLELLEAGRSVEAVDLIEEVAEKYREELSKKYGEKDPAVLKQKVAADAFTLLEEG",
    "pmhc_g12d_01": "DGAAAAALLLTEFAKKYGPDAPLPPELIEQILALATGATREAIRAVIEKNPNATPAEIGAAFAAAPA",
    "pmhc_g12d_22": "SEEEEERAAALVAEAGRKAEELLKKGKEEEAEALLEETLKELEELVKDPELLAAMKAALEDQYEYWKRRYR",
    "pmhc_g12d_09": "GPLTVDEEAALVAELSRYQVEHADDPEAMERLLPIIKTLYEGMNAEDVAERRALLQKAKEQLDALLAE",
    "pmhc_g12d_16": "RTIPATTAWGAEAYVVRAETKEEAEKIAKELSKKTGKAAKVVQMNDTEFLIILADTEEEAEALAREVLGE",
    "pmhc_g12d_00": "MKKFKVTAAGKEFTIEAPSPEEAAREALKKTPGLGELISLTQVDDNTYLAQARMPGPNGGPDRLVGLVVEEIEE",
    "pmhc_g12d_19": "GREELLREVIEELGLTIEDVERLIKFLEDNGGDRGSINFYRDMLEVLKEEEEREK",
    "pmhc_g12d_15": "GVEKLLEAAEAGDADGVAEALGEILKFLDKLEGEERAELERLLGEGGFALGSAATLAEGRAALAEAAKRLKEKLE",
    "pmhc_g12d_14": "AEPQVIQNADLQPGDRLVCGVKGFHPSDIEVELLKNGERIPNVEYSPLSLNSDGSFFILAWTDKTVTSTDVYSCRVRHP",
    "pmhc_g12d_20": "DEERLRELIEEARANKDEIAALLKRFFESEEPSEEELAAVRLLAELGRLANLTRDGPLRNDALTVKLEILEQAAK",
    "pmhc_g12d_04": "GRWEVPEEVIKAVEPAVPLALRAHEDPDDLEAAREAAALRDKVAEEFGVTRDEVGNALGIALEKRLSK",
    "pmhc_g12d_21": "QRTPKIQVYTRVPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHTPLSLNSDGSFYLLYYTEFDTKDTYTCR",
    "pmhc_g12d_18": "IQRTPKIQVYSRVPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHTDLSFNSDGSFYLLYYTEFTPTEKDEYTCR",
    "pmhc_g12d_11": "IQRTPKIQVYSRVPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLAFNSDGSFYLLYYTEFTPTEKDEYTCRV",
    "pmhc_g12d_23": "AETKIKFSEVRQLSCGRYYVVAVSGIGELTIPEVEVSCERALKAKTGETVDTRVLDDTKLLRVYVAETEEAARALI",
    "pmhc_g12d_12": "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHTDLSFNSDGSFYLLYYTEFTPTETDEYTCRVKHPS",
    "pmhc_g12d_03": "QRTPKIQVYTRAPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHTDLSFNADGSFYLLYYTEFDTKDTYTCR",
    "pmhc_g12d_02": "IQRTPKIQVYTRAPAENGKPNFLNCLVSGFHPSDIEVDLLKNGERIEKVEHTDLSFNSDGSFYLLYYTEFDPKDEYACRV",
}

TIMESTEPS = 150
SAMPLES = 2
RECYCLES = 3


def main() -> None:
    from chai_lab.chai1 import run_inference

    out = Path(sys.argv[1] if len(sys.argv) > 1 else "/kaggle/working/mono")
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "monomer_refold.csv"

    done = set()
    if csv_path.exists():
        done = {r["design"] for r in csv.DictReader(csv_path.open())}
        print(f"resuming - {len(done)} already folded")

    fields = ["design", "length", "mono_plddt", "mono_ptm", "agg"]
    for n, (name, seq) in enumerate(SEQUENCES.items(), start=1):
        if name in done:
            continue
        print(f"[{n}/{len(SEQUENCES)}] {name} ({len(seq)} aa)", flush=True)
        fasta = out / f"{name}.fasta"
        fasta.write_text(f">protein|name=binder\n{seq}\n")
        rundir = out / name
        if rundir.exists():
            shutil.rmtree(rundir)
        try:
            cand = run_inference(
                fasta_file=fasta, output_dir=rundir,
                use_esm_embeddings=False, use_msa_server=False,
                num_trunk_recycles=RECYCLES, num_diffn_timesteps=TIMESTEPS,
                num_diffn_samples=SAMPLES, seed=42, device="cuda:0",
            )
        except Exception as error:
            print(f"    FAILED: {error}", flush=True)
            continue

        scores = [float(r.aggregate_score.item()) for r in cand.ranking_data]
        i = max(range(len(scores)), key=scores.__getitem__)
        r = cand.ranking_data[i]
        rec = {
            "design": name,
            "length": len(seq),
            "mono_plddt": round(float(r.plddt_scores.complex_plddt.flatten()[0]), 4),
            "mono_ptm": round(float(r.ptm_scores.complex_ptm.flatten()[0]), 4),
            "agg": round(scores[i], 4),
        }
        exists = csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            if not exists:
                writer.writeheader()
            writer.writerow(rec)
        print(f"    plddt {rec['mono_plddt']:.3f}  ptm {rec['mono_ptm']:.3f}", flush=True)

    print(f"done -> {csv_path}")


if __name__ == "__main__":
    main()
