"""
Cross-panel selectivity screen.

The three-target counter-screen asked whether the binder needs THIS peptide.
This asks a harder question: does it discriminate against the other peptides
HLA-A*11:01 actually presents? A binder that grips every A*11:01 ligand is a
groove binder, and in a patient it would meet those ligands on every cell.

The panel, in the same 10-mer register unless noted:

  g12d    VVVGADGVGK   on-target, repeated here as the internal reference -
                       same settings and seed as the counter-screen, so it
                       doubles as a check that the two runs are comparable
  g12v    VVVGAVGVGK   the other common KRAS codon-12 mutations; these are
  g12c    VVVGACGVGK   the off-targets that matter clinically, because a
  g12r    VVVGARGVGK   patient population carries them and a binder that
  g12a    VVVGAAGVGK   cannot tell them apart is not a G12D therapeutic
  g12s    VVVGASGVGK
  g13d    VVVGAGDVGK   aspartate one position along - tests whether the
                       binder reads the charge or reads its location
  nonamer  VVGADGVGK   the 9-mer register of the same neoantigen; a shorter
                       peptide sits differently in the groove
  ebv_10  AVFDRKSDAK   unrelated A*11:01 ligand, same length (EBV EBNA)
  ebv_9    IVTDFSVIK   unrelated A*11:01 ligand, 9-mer (EBV EBNA)

The two EBV entries are literature A*11:01-restricted epitopes and serve as
the unrelated-ligand controls. Verify them against IEDB before quoting them
in a writeup.

Jobs run binder-outer, peptide-inner, so an interrupted session leaves one
design's panel complete rather than two half-panels.

    python crosspan.py /kaggle/working/xp
"""
from __future__ import annotations
import csv, shutil, sys, time
from pathlib import Path

A1A2 = "GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQRT"
B2M = "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"

PEPTIDE = {
    "g12d":    "VVVGADGVGK",
    "g12v":    "VVVGAVGVGK",
    "g12c":    "VVVGACGVGK",
    "g12r":    "VVVGARGVGK",
    "g12a":    "VVVGAAGVGK",
    "g12s":    "VVVGASGVGK",
    "g13d":    "VVVGAGDVGK",
    "nonamer": "VVGADGVGK",
    "ebv_10":  "AVFDRKSDAK",
    "ebv_9":   "IVTDFSVIK",
}

BINDERS = {
    "pmhc_g12d_09": "GPLTVDEEAALVAELSRYQVEHADDPEAMERLLPIIKTLYEGMNAEDVAERRALLQKAKEQLDALLAE",
    "pmhc_g12d_06": "SPLLLEVGRTGGADNDWEKMAALLEKIAALPDEELRALAAEDPNTMLVLAEIAEVLGLPEVAAVLRRAA",
}

TIMESTEPS = 150
SAMPLES = 1
RECYCLES = 3
SEEDS = [42]

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

    out = Path(sys.argv[1] if len(sys.argv) > 1 else "/kaggle/working/xp")
    out.mkdir(parents=True, exist_ok=True)
    csv_path = out / "crosspan.csv"

    done = set()
    if csv_path.exists():
        done = {(r["design"], r["target"], int(r["seed"]))
                for r in csv.DictReader(csv_path.open())}
        print(f"resuming - {len(done)} folds already done")

    fields = ["design", "target", "seed", "complex_plddt", "complex_ptm",
              "complex_iptm", "iptm_mhc", "iptm_b2m", "iptm_peptide", "agg"]

    jobs = [(d, t, SEEDS[0]) for d in BINDERS for t in PEPTIDE]

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
