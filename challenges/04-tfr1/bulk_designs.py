"""Bulk RFD4-Proteina generation + sequence-level quality filtering.

Run in Colab (CPU). Put request.json, request_h1.json, request_h2.json in /content.
Set API_KEY from an env var; do not paste a key into a file you may share.

What it does:
  1. Calls the hosted endpoint several times per epitope (fresh samples each call)
  2. Strips the fixed target sequence off every returned chain
  3. Drops low-complexity artifacts (poly-A runs, short tandem repeats)
  4. Ranks the survivors and writes submission.csv plus a FASTA for folding
"""
import os, json, math, urllib.request, urllib.error, pathlib, csv, collections

API_KEY = os.environ["NVIDIA_API_KEY"]          # export it; never hard-code
URL = "https://health.api.nvidia.com/v1/biology/ipd/rfd4-proteina/generative"

REQUESTS = ["/content/request.json",            # epitope set 0: N206/L209/K374
            "/content/request_h1.json",         # set 1: K240/K241/E244/E266
            "/content/request_h2.json"]         # set 2: K287/P289/V291
ROUNDS = 2                                      # calls per request file
OUT = pathlib.Path("/content/rfd4_out"); OUT.mkdir(exist_ok=True)

# target sequence is stripped off each returned chain
TARGET = ("RLYWDDLKRKLSEKLDSTDFTGTIKLLNENSYVPREAGSQKDENLALYVENQFREFKLSKVWRDQHFVKI"
          "QVKDSAQNSVIIVDKNGRLVYLVENPGGYVAYSKAATVTGKLVHANFGTKKDFEDLYTPVNGSIVIVRAG"
          "KITFAEKVANAESLNAIGVLIYMDQTKFPIVNAELSFFGHAHLGTGDPYTPGFPSFNHTQFPPSRSSGLP"
          "NIPVQTISRAAAEKLFGNMEGDCPSDWKTDSTCRMVTSESKNVKLTVSNVLKE")


def call(path):
    req = urllib.request.Request(
        URL, data=pathlib.Path(path).read_bytes(),
        headers={"Authorization": "Bearer " + API_KEY,
                 "Content-Type": "application/json",
                 "NVCF-POLL-SECONDS": "300"},
        method="POST")
    try:
        with urllib.request.urlopen(req, timeout=600) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        print("  HTTP", e.code, e.read().decode()[:300])
        return None


def strip_target(seq):
    # the model returns binder chain followed by the fixed target chain
    for cut in range(len(TARGET), 150, -1):
        if seq.endswith(TARGET[-cut:]):
            return seq[:-cut]
    return seq


def entropy(s):
    c = collections.Counter(s)
    n = len(s)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def longest_run(s):
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


def kmer_uniqueness(s, k=6):
    if len(s) < k:
        return 1.0
    km = [s[i:i + k] for i in range(len(s) - k + 1)]
    return len(set(km)) / len(km)


def quality(s):
    """Return (ok, reasons). Pure sequence heuristics - no structure involved."""
    bad = []
    if not (40 <= len(s) <= 100):
        bad.append("length %d outside minibinder stratum" % len(s))
    top = collections.Counter(s).most_common(1)[0]
    if top[1] / len(s) > 0.30:
        bad.append("%s is %.0f%% of sequence" % (top[0], 100 * top[1] / len(s)))
    if entropy(s) < 3.2:
        bad.append("low entropy %.2f" % entropy(s))
    if longest_run(s) >= 5:
        bad.append("homopolymer run %d" % longest_run(s))
    if kmer_uniqueness(s) < 0.80:
        bad.append("tandem repeats (6-mer uniqueness %.2f)" % kmer_uniqueness(s))
    return (not bad), bad


seen, kept, rejected = set(), [], []
for path in REQUESTS:
    tag = pathlib.Path(path).stem
    for rnd in range(ROUNDS):
        print("calling %s round %d ..." % (tag, rnd + 1))
        res = call(path)
        if not res:
            continue
        for d in res["designs"]:
            s = strip_target(d["sequence"])
            if s in seen:
                continue
            seen.add(s)
            (OUT / ("%s_%s.cif" % (tag, d["design_id"]))).write_text(
                d["structure"]["content"])
            ok, why = quality(s)
            (kept if ok else rejected).append((tag, d["design_id"], s, why))

print("\nkept %d | rejected %d" % (len(kept), len(rejected)))
for t, i, s, why in rejected:
    print("  REJECT %-14s %-10s %s" % (t, i, "; ".join(why)))

# rank by entropy then k-mer diversity: crude proxies for a real, designable fold
kept.sort(key=lambda r: (entropy(r[2]), kmer_uniqueness(r[2])), reverse=True)
final = kept[:20]

with open("/content/submission.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["name", "sequence", "description"])
    for n, (tag, did, s, _) in enumerate(final):
        w.writerow(["tfr1_%s_%02d" % (tag.replace("request", "ep"), n), s,
                    "De novo minibinder, RFD4-Proteina (NVIDIA/IPD), conditioned on the "
                    "TfR1 apical domain (UniProt P02786 residues 189-383, PDB 11QC); "
                    "mouse-conserved hotspot patch."])

with open("/content/top_designs.fasta", "w") as f:
    for n, (tag, did, s, _) in enumerate(final):
        f.write(">tfr1_%s_%02d\n%s\n" % (tag.replace("request", "ep"), n, s))

print("\nwrote submission.csv (%d designs) and top_designs.fasta" % len(final))
for n, (tag, did, s, _) in enumerate(final[:10]):
    print("  %2d  %-12s len=%2d  H=%.2f  %s" % (n, tag, len(s), entropy(s), s))
