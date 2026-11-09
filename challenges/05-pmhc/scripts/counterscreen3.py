"""
Three-target counter-screen: negative selection.

Each binder is co-folded against three analytes and the confidence in the
binder-to-target interface is compared across them.

    g12d   binder + alpha1-alpha2 + beta2m + VVVGADGVGK   on-target
    wt     binder + alpha1-alpha2 + beta2m + VVVGAGGVGK   mutation specificity
    apo    binder + alpha1-alpha2 + beta2m                peptide dependence

The scoring rubric is lexicographic: a binder that also grips the empty
groove is worth less than a weaker one that needs the peptide, however good
its on-target number is. The on-target fold is recomputed here rather than
reused from BoltzGen, because a difference between two predictors is not a
difference between two targets.

WHY THE EARLIER ATTEMPT WAS INVALID. It ran with use_msa_server=False. The
MHC heavy chain and beta2-microglobulin are natural proteins and Chai needs
their MSAs; without them the whole complex scored near zero, including
complexes that Boltz-2 placed at ipTM 0.69. Internet must be ON in the Kaggle
notebook for this to work - Settings -> Internet -> on. If it is off, Chai
will fail to reach the MSA server and the run is worthless again.

Self-contained: all sequences are embedded. Resumable: one CSV row is
appended per fold and completed folds are skipped on restart.

    python counterscreen3.py /kaggle/working/cs3
"""
from __future__ import annotations
import csv, shutil, sys, time
from pathlib import Path

A1A2 = "GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQRT"
B2M = "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
PEPTIDE = {"g12d": "VVVGADGVGK", "wt": "VVVGAGGVGK", "apo": None}

BINDERS = {
    "pmhc_g12d_06": "SPLLLEVGRTGGADNDWEKMAALLEKIAALPDEELRALAAEDPNTMLVLAEIAEVLGLPEVAAVLRRAA",
    "pmhc_g12d_09": "GPLTVDEEAALVAELSRYQVEHADDPEAMERLLPIIKTLYEGMNAEDVAERRALLQKAKEQLDALLAE",
}

TIMESTEPS = 150
SAMPLES = 1
RECYCLES = 3
SEEDS = [42, 7]

# One diffusion sample rather than two. This is a COMPARATIVE screen: the
# same seed and the same settings are applied to all three analytes of a
# design, so what is lost is the best-of-two on each absolute number, not the
# validity of the difference between them. Two samples would roughly double
# the wall clock and the session will not survive it.

# Designs are ordered by on-target rank and each design's three analytes run
# consecutively, so an interrupted session always leaves complete triples
# rather than fragments that cannot be compared.

# The top design's triple is then repeated under a second seed before the
# next design starts. A specificity gap that survives a change of seed is a
# property of the structures; one that does not is a property of the sampler.
# Job order:  _06 seed 42 triple, _06 seed 7 triple, _09 seed 42 triple.


def write_fasta(path: Path, binder: str, peptide: str | None) -> None:
    parts = [f">protein|name=binder\n{binder}",
             f">protein|name=mhc\n{A1A2}",
             f">protein|name=b2m\n{B2M}"]
    if peptide:
        parts.append(f">protein|name=peptide\n{peptide}")
    path.write_text("\n".join(parts) + "\n")


def scalar(obj, *names):
    """First attribute that exists and reduces to a float. Chai renames these
    between versions, so never address one by a single hard-coded name."""
    for n in names:
        v = getattr(obj, n, None)
        if v is None:
            continue
        try:
            return round(float(v.flatten()[0]), 4)
        except Exception:
            try:
                return round(float(v), 4)
            except Exception:
                pass
    return None


def pair_iptm(rank, n_chains: int):
    """Binder is chain 0. Return its iptm against each other chain."""
    try:
        m = rank.ptm_scores.per_chain_pair_iptm.flatten().reshape(n_chains, n_chains)
        return [round(float(m[0, j]), 4) for j in range(1, n_chains)]
    except Exception:
        return [None] * (n_chains - 1)


def main() -> None:
    from chai_lab.chai1 import run_inference

    out = Path(sys.argv[1] if len(sys.argv) > 1 else "/kaggle/working/cs3")
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "counterscreen3.csv"

    done = set()
    if csv_path.exists():
        done = {(r["design"], r["target"], int(r["seed"]))
                for r in csv.DictReader(csv_path.open())}
        print(f"resuming - {len(done)} folds already done")

    fields = ["design", "target", "seed", "complex_plddt", "complex_ptm",
              "complex_iptm", "iptm_mhc", "iptm_b2m", "iptm_peptide", "agg"]

    order = list(BINDERS)
    jobs = [(order[0], t, s) for s in SEEDS for t in ("g12d", "apo", "wt")]
    jobs += [(d, t, SEEDS[0]) for d in order[1:] for t in ("g12d", "apo", "wt")]

    printed_header = False
    printed_attrs = False
    for n, (design, target, seed) in enumerate(jobs, start=1):
        if (design, target, seed) in done:
            continue
        t0 = time.time()
        peptide = PEPTIDE[target]
        n_chains = 4 if peptide else 3
        print(f"[{n}/{len(jobs)}] {design} vs {target} seed {seed}", flush=True)

        fasta = out / f"{design}_{target}_s{seed}.fasta"
        write_fasta(fasta, BINDERS[design], peptide)
        rundir = out / f"{design}_{target}_s{seed}"
        if rundir.exists():
            shutil.rmtree(rundir)

        try:
            cand = run_inference(
                fasta_file=fasta, output_dir=rundir,
                use_esm_embeddings=True, use_msa_server=True,
                num_trunk_recycles=RECYCLES, num_diffn_timesteps=TIMESTEPS,
                num_diffn_samples=SAMPLES, seed=seed, device="cuda:0",
            )
        except Exception as error:
            print(f"    FAILED: {type(error).__name__}: {error}", flush=True)
            continue

        try:
            scores = [float(r.aggregate_score.item()) for r in cand.ranking_data]
            best = cand.ranking_data[max(range(len(scores)), key=scores.__getitem__)]
            if not printed_attrs:
                print("    ptm fields:", [a for a in dir(best.ptm_scores)
                                          if not a.startswith("_")], flush=True)
                printed_attrs = True
            pairs = pair_iptm(best, n_chains)
            rec = {
                "design": design, "target": target, "seed": seed,
                "complex_plddt": scalar(best.plddt_scores, "complex_plddt"),
                "complex_ptm": scalar(best.ptm_scores, "complex_ptm"),
                "complex_iptm": scalar(best.ptm_scores, "interface_ptm",
                                       "complex_iptm", "complex_interface_ptm"),
                "iptm_mhc": pairs[0],
                "iptm_b2m": pairs[1],
                "iptm_peptide": pairs[2] if peptide else None,
                "agg": round(max(scores), 4),
            }
        except Exception as error:
            # A scoring failure must never discard a fold that already ran.
            print(f"    SCORING FAILED: {type(error).__name__}: {error}", flush=True)
            continue
        exists = csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=fields)
            if not exists:
                w.writeheader()
            w.writerow(rec)

        # Echo the row to stdout as CSV as well. If the session dies and the
        # filesystem goes with it, the console scrollback still holds every
        # completed measurement and the table can be rebuilt by hand.
        if not printed_header:
            print("CSV," + ",".join(fields), flush=True)
            printed_header = True
        print("CSV," + ",".join("" if rec[f] is None else str(rec[f])
                                for f in fields), flush=True)
        print(f"    {time.time() - t0:.0f}s elapsed", flush=True)

    print(f"done -> {csv_path}")


if __name__ == "__main__":
    main()
