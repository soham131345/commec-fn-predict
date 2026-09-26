"""
Copyright (c) 2026 Soham Bhole. All rights reserved.

Drift experiment: the paper's central figure.

Takes held-out test POSITIVES, introduces increasing rates of random amino
acid substitutions (simulating evolutionary/engineered distance from every
known sequence), and measures whether the trained model still flags them at
its calibrated screening threshold. Output: recall vs. sequence identity
curve (PNG + CSV). If the curve stays high where BLAST/HMM homology search
would fail (< ~70% identity), the homology-independence claim is proven.

Usage:
  python drift_experiment.py --model-dir models/kmer \
      --test-fasta splits/test.fasta --test-labels splits/test_labels.csv \
      -o drift/
"""

import argparse
import csv
import json
import os
import random
import sys

import joblib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fn_common

AA = fn_common.ALPHABET


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


def read_labels(path):
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            out[row["accession"]] = int(row["label"])
    return out


def mutate(seq, rate, rng):
    seq = list(seq)
    for i, c in enumerate(seq):
        if c in AA and rng.random() < rate:
            seq[i] = rng.choice([a for a in AA if a != c])
    return "".join(seq)


def score(seqs, model_dir):
    calib = json.load(open(os.path.join(model_dir, "calibration.json")))
    model = joblib.load(os.path.join(model_dir, "model.joblib"))
    if calib.get("model_type") == "esm2_embedding_logreg":
        import torch
        from transformers import AutoModel, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(calib["backbone"])
        lm = AutoModel.from_pretrained(calib["backbone"]).eval()
        feats = []
        with torch.no_grad():
            for i in range(0, len(seqs), 16):
                chunk = [s[:1022] for s in seqs[i:i + 16]]
                enc = tok(chunk, return_tensors="pt", padding=True, truncation=True)
                out = lm(**enc).last_hidden_state
                mask = enc["attention_mask"].unsqueeze(-1)
                feats.append(((out * mask).sum(1) / mask.sum(1).clamp(min=1)).numpy())
        X = np.vstack(feats)
    else:
        X = seqs
    return model.predict_proba(X)[:, 1], float(calib["threshold"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model-dir", required=True)
    p.add_argument("--test-fasta", required=True)
    p.add_argument("--test-labels", required=True)
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--rates", default="0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8")
    p.add_argument("--replicates", type=int, default=5)
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    seqs = read_fasta(args.test_fasta)
    labels = read_labels(args.test_labels)
    positives = [s for a, s in seqs.items() if labels.get(a) == 1]
    if not positives:
        raise SystemExit("No positive sequences in test set.")
    print(f"{len(positives)} test positives")

    rng = random.Random(args.seed)
    rates = [float(x) for x in args.rates.split(",")]
    rows = []
    for rate in rates:
        recalls = []
        for rep in range(args.replicates):
            mutated = [mutate(s, rate, rng) for s in positives]
            scores, thr = score(mutated, args.model_dir)
            recalls.append(float((scores >= thr).mean()))
        rows.append({
            "mutation_rate": rate,
            "mean_identity_pct": round(100 * (1 - rate), 1),
            "recall_mean": float(np.mean(recalls)),
            "recall_min": float(np.min(recalls)),
            "recall_max": float(np.max(recalls)),
        })
        print(f"identity {100*(1-rate):5.1f}%  recall "
              f"{np.mean(recalls):.3f} [{np.min(recalls):.3f}-{np.max(recalls):.3f}]")

    csv_path = os.path.join(args.outdir, "drift_curve.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = [r["mean_identity_pct"] for r in rows]
    y = [r["recall_mean"] for r in rows]
    lo = [r["recall_min"] for r in rows]
    hi = [r["recall_max"] for r in rows]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, y, marker="o")
    ax.fill_between(x, lo, hi, alpha=0.2)
    ax.axvspan(70, 100, color="green", alpha=0.07,
               label="homology search works here")
    ax.axvspan(0, 70, color="red", alpha=0.07,
               label="homology search degrades")
    ax.set_xlabel("identity to nearest known sequence (%)")
    ax.set_ylabel(f"recall at calibrated threshold")
    ax.set_ylim(0, 1.05)
    ax.invert_xaxis()
    ax.legend()
    fig.tight_layout()
    png_path = os.path.join(args.outdir, "drift_curve.png")
    fig.savefig(png_path, dpi=200)
    print(f"Wrote {csv_path} and {png_path}")


if __name__ == "__main__":
    main()
