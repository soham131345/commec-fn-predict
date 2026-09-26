"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""
Evaluate one or more trained models on the homology-held-out test split and
emit a screening-oriented benchmark report.

Metrics per model:
  AUROC, AUPRC, recall at FPR <= 1e-2 and <= 1e-3 (screening operating
  points), Brier score, expected calibration error (ECE, 10 bins), and the
  calibrated threshold's test recall/FPR.

Usage (k-mer model):
  python evaluate.py --test-fasta test.fasta --test-labels test_labels.csv \
      --model kmer_baseline=<dir-with-model.joblib> -o report.md

Usage (ESM model; pass cached test embeddings):
  python evaluate.py --test-fasta test.fasta --test-labels test_labels.csv \
      --model esm=<dir> --esm-test-embeddings <dir>/embeddings_test.npy \
      --esm-test-ids <dir>/ids_test.txt -o report.md

Requires: scikit-learn, joblib, numpy.
"""

import argparse
import json
import os
import sys

import joblib
import numpy as np
from sklearn.metrics import (average_precision_score, brier_score_loss,
                             roc_auc_score)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fn_common

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
    import csv
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            out[row["accession"]] = int(row["label"])
    return out

def recall_at_fpr(y, s, fpr_cap):
    order = np.argsort(-s)
    y = y[order]
    n_neg = max(1, (y == 0).sum())
    n_pos = max(1, (y == 1).sum())
    fp = np.cumsum(y == 0)
    ok = fp / n_neg <= fpr_cap
    if not ok.any():
        return 0.0
    idx = np.where(ok)[0][-1]
    return float((y[:idx + 1] == 1).sum() / n_pos)

def ece(y, p, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    err = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (p >= lo) & (p < hi)
        if m.any():
            err += m.mean() * abs(p[m].mean() - y[m].mean())
    return float(err)

def metrics(y, s, thr=None):
    out = {
        "AUROC": roc_auc_score(y, s),
        "AUPRC": average_precision_score(y, s),
        "recall@FPR1e-2": recall_at_fpr(y, s, 1e-2),
        "recall@FPR1e-3": recall_at_fpr(y, s, 1e-3),
        "Brier": brier_score_loss(y, s),
        "ECE": ece(y, s),
    }
    if thr is not None:
        pred = s >= thr
        out["threshold"] = thr
        out["test_recall@thr"] = float((pred & (y == 1)).sum() / max(1, (y == 1).sum()))
        out["test_FPR@thr"] = float((pred & (y == 0)).sum() / max(1, (y == 0).sum()))
    return out

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--test-fasta", required=True)
    p.add_argument("--test-labels", required=True)
    p.add_argument("--model", action="append", required=True,
                   help="name=<model-dir>; repeatable, one per model to benchmark.")
    p.add_argument("--esm-test-embeddings", default=None)
    p.add_argument("--esm-test-ids", default=None)
    p.add_argument("-o", "--out", default="benchmark_report.md")
    args = p.parse_args()

    seqs, labels = read_fasta(args.test_fasta), read_labels(args.test_labels)
    ids = [a for a in seqs if a in labels]
    y = np.array([labels[a] for a in ids])
    print(f"Test set: {len(ids)} sequences ({y.sum()} pos / {(y == 0).sum()} neg)")

    rows = {}
    for spec in args.model:
        name, mdir = spec.split("=", 1)
        model = joblib.load(os.path.join(mdir, "model.joblib"))
        calib = {}
        cal_path = os.path.join(mdir, "calibration.json")
        if os.path.exists(cal_path):
            calib = json.load(open(cal_path))
        if calib.get("model_type") == "esm2_embedding_logreg":
            if not args.esm_test_embeddings or not args.esm_test_ids:
                raise SystemExit(f"{name}: pass --esm-test-embeddings/--esm-test-ids")
            E = np.load(args.esm_test_embeddings)
            emb_ids = open(args.esm_test_ids).read().split()
            if emb_ids != ids:

                pos = {a: i for i, a in enumerate(emb_ids)}
                E = E[[pos[a] for a in ids]]
            s = model.predict_proba(E)[:, 1]
        else:
            s = model.predict_proba([seqs[a] for a in ids])[:, 1]
        rows[name] = metrics(y, s, calib.get("threshold"))

    cols = ["AUROC", "AUPRC", "recall@FPR1e-2", "recall@FPR1e-3",
            "Brier", "ECE", "test_recall@thr", "test_FPR@thr"]
    lines = ["# Function-prediction benchmark (homology-held-out test set)", "",
             f"Test set: {len(ids)} sequences ({int(y.sum())} positives).", "",
             "| model | " + " | ".join(cols) + " |",
             "|---|" + "---|" * len(cols)]
    for name, m in rows.items():
        lines.append("| " + name + " | " + " | ".join(
            f"{m[c]:.4f}" if isinstance(m.get(c), float) else "-" for c in cols) + " |")
    report = "\n".join(lines) + "\n"
    with open(args.out, "w") as fh:
        fh.write(report)
    print(report)

if __name__ == "__main__":
    main()
