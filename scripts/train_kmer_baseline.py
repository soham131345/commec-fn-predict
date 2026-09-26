"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""
Train the k-mer + logistic-regression baseline for toxin/virulence
function prediction (protein sequences).

Model: amino-acid k-mer counts (default k=4), L2-normalized, multinomial
logistic regression with balanced class weights. Deliberately simple: this
is the benchmark floor that any embedding model must beat on the
homology-held-out test set.

Outputs in <model-dir>:
  model.joblib        - sklearn Pipeline (vectorizer + classifier)
  calibration.json    - threshold chosen at --target-fpr on the validation set
Requires: scikit-learn, joblib.
"""

import argparse
import json
import os
import sys

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fn_common import KmerVectorizer

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

def load_split(fasta, labels_csv):
    seqs, labels = read_fasta(fasta), read_labels(labels_csv)
    ids = [a for a in seqs if a in labels]
    return [seqs[a] for a in ids], np.array([labels[a] for a in ids])

def threshold_at_fpr(y_true, scores, target_fpr):
    """Largest score threshold whose FPR <= target_fpr."""
    order = np.argsort(-scores)
    y = y_true[order]
    s = scores[order]
    n_neg = max(1, (y_true == 0).sum())
    fp = np.cumsum(y == 0)
    ok = fp / n_neg <= target_fpr
    if not ok.any():
        return float(s.max() + 1e-6), 0.0
    idx = np.where(ok)[0][-1]
    recall = (y[:idx + 1] == 1).sum() / max(1, (y_true == 1).sum())
    return float(s[idx]), float(recall)

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--train-fasta", required=True)
    p.add_argument("--train-labels", required=True)
    p.add_argument("--val-fasta", required=True)
    p.add_argument("--val-labels", required=True)
    p.add_argument("-o", "--model-dir", required=True)
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--C", type=float, default=1.0)
    p.add_argument("--target-fpr", type=float, default=0.01,
                   help="Validation FPR at which to fix the screening threshold.")
    args = p.parse_args()

    import os
    os.makedirs(args.model_dir, exist_ok=True)

    Xtr, ytr = load_split(args.train_fasta, args.train_labels)
    Xva, yva = load_split(args.val_fasta, args.val_labels)
    print(f"train: {len(Xtr)} ({ytr.sum()} pos), val: {len(Xva)} ({yva.sum()} pos)")

    model = Pipeline([
        ("kmer", KmerVectorizer(k=args.k)),
        ("clf", LogisticRegression(C=args.C, class_weight="balanced",
                                   max_iter=2000, n_jobs=-1)),
    ])
    model.fit(Xtr, ytr)
    val_scores = model.predict_proba(Xva)[:, 1]
    thr, rec = threshold_at_fpr(yva, val_scores, args.target_fpr)

    joblib.dump(model, os.path.join(args.model_dir, "model.joblib"))
    with open(os.path.join(args.model_dir, "calibration.json"), "w") as fh:
        json.dump({"model_type": "kmer_logreg", "k": args.k,
                   "target_fpr": args.target_fpr, "threshold": thr,
                   "val_recall_at_threshold": rec}, fh, indent=2)
    print(f"Saved. Val threshold {thr:.4f} @ FPR<={args.target_fpr} "
          f"(val recall {rec:.3f}).")

if __name__ == "__main__":
    main()
