"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""
Build homology-held-out train/val/test splits.

Clusters all sequences with MMseqs2 (easy-cluster, --min-seq-id 0.3 by
default), then assigns ENTIRE clusters to one split, so no test sequence
shares a cluster with any training sequence. This is the split discipline
that makes "works regardless of homology" claims defensible.

Requires: mmseqs (conda install -c conda-forge -c bioconda mmseqs2).

Outputs in <outdir>: clusters.tsv, splits.csv, train.fasta, val.fasta,
test.fasta, train_labels.csv, val_labels.csv, test_labels.csv
"""

import argparse
import csv
import os
import random
import shutil
import subprocess
import sys
from collections import defaultdict

def read_fasta(path):
    seqs, header, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if header:
        seqs[header] = "".join(buf)
    return seqs

def write_fasta(seqs, path, ids):
    with open(path, "w") as fh:
        for acc in ids:
            fh.write(f">{acc}\n")
            for i in range(0, len(seqs[acc]), 60):
                fh.write(seqs[acc][i:i + 60] + "\n")

def cluster_edlib(allseq, ident):
    from collections import Counter
    import edlib
    ids_sorted = sorted(allseq, key=lambda a: -len(allseq[a]))
    reps, rep_len, clusters = [], {}, {}
    index = defaultdict(list)
    for acc in ids_sorted:
        s = allseq[acc]
        L = len(s)
        hits = Counter()
        for i in range(L - 2):
            for ri in index.get(s[i:i + 3], ()):
                hits[ri] += 1
        assigned = False
        for ri, h in hits.most_common(60):
            rl = rep_len[ri]
            if min(L, rl) / max(L, rl) < 0.7:
                continue
            if h < max(2, int(0.10 * (min(L, rl) - 2))):
                break
            r = edlib.align(s, reps[ri], mode="NW", task="distance")
            if 1.0 - r["editDistance"] / max(L, rl) >= ident:
                clusters[ri].append(acc)
                assigned = True
                break
        if not assigned:
            ri = len(reps)
            reps.append(s)
            rep_len[ri] = L
            clusters[ri] = [acc]
            for i in range(L - 2):
                index[s[i:i + 3]].append(ri)
    return {ri: members for ri, members in clusters.items()}

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positives", required=True)
    p.add_argument("--negatives", required=True)
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--min-seq-id", type=float, default=0.3)
    p.add_argument("--val-frac", type=float, default=0.1)
    p.add_argument("--test-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=13)
    p.add_argument("--python-fallback", action="store_true",
                   help="Use built-in CD-HIT-style clustering (edlib) instead of "
                        "mmseqs. Slower; use when mmseqs is unavailable (e.g. Windows). "
                        "Requires: pip install edlib")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    pos, neg = read_fasta(args.positives), read_fasta(args.negatives)
    labels = {**{a: 1 for a in pos}, **{a: 0 for a in neg}}
    allseq = {**pos, **neg}

    if args.python_fallback:
        clusters = cluster_edlib(allseq, args.min_seq_id)
    else:
        if not shutil.which("mmseqs"):
            sys.exit("mmseqs not found. Install: conda install -c bioconda mmseqs2, "
                     "or re-run with --python-fallback (pip install edlib).")
        combined = os.path.join(args.outdir, "_all.fasta")
        write_fasta(allseq, combined, list(allseq))
        clu_dir = os.path.join(args.outdir, "_mmseqs")
        subprocess.run(
            ["mmseqs", "easy-cluster", combined,
             os.path.join(args.outdir, "_clu"), clu_dir,
             "--min-seq-id", str(args.min_seq_id), "-c", "0.5",
             "--cov-mode", "1"],
            check=True)
        clusters = defaultdict(list)
        with open(os.path.join(args.outdir, "_clu_cluster.tsv")) as fh:
            for line in fh:
                rep, member = line.rstrip("\n").split("\t")
                clusters[rep].append(member)
        shutil.rmtree(clu_dir, ignore_errors=True)

    rng = random.Random(args.seed)
    cluster_list = list(clusters.values())
    rng.shuffle(cluster_list)

    cluster_list.sort(key=len, reverse=True)

    def cluster_pos_frac(members):
        return sum(labels[m] for m in members) / len(members)

    splits = {"train": [], "val": [], "test": []}
    counts = {s: [0, 0] for s in splits}
    target = {"test": args.test_frac, "val": args.val_frac}
    total = len(allseq)
    for members in cluster_list:
        n_pos = sum(labels[m] for m in members)
        n_neg = len(members) - n_pos

        best, best_deficit = "train", -1.0
        for s, frac in target.items():
            deficit = frac * total - (counts[s][0] + counts[s][1])
            if deficit > best_deficit:
                best, best_deficit = s, deficit
        splits[best].extend(members)
        counts[best][0] += n_pos
        counts[best][1] += n_neg

    print(f"Clusters: {len(cluster_list)} at min-seq-id {args.min_seq_id}")
    for s in splits:
        n = counts[s]
        print(f"  {s}: {sum(n)} seqs ({n[0]} pos / {n[1]} neg)")

    with open(os.path.join(args.outdir, "splits.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "label", "split"])
        for s, ids in splits.items():
            for a in ids:
                w.writerow([a, labels[a], s])
    for s, ids in splits.items():
        write_fasta(allseq, os.path.join(args.outdir, f"{s}.fasta"), ids)
        with open(os.path.join(args.outdir, f"{s}_labels.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["accession", "label"])
            for a in ids:
                w.writerow([a, labels[a]])

if __name__ == "__main__":
    main()
