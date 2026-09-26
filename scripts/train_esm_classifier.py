"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""
Train an ESM-2 embedding + classifier model for toxin/virulence function
prediction, benchmarked against the k-mer baseline on the same
homology-held-out splits.

Pipeline: embed each protein with a pretrained ESM-2 (mean-pool the last
hidden layer over residues), cache embeddings, fit logistic regression with
balanced class weights on top. Frozen backbone = cheap, reproducible, and
hard to overfit on small labelled sets.

Outputs in <model-dir>:
  embeddings_train.npy / embeddings_val.npy / embeddings_test.npy (cached)
  ids_*.txt            - accessions aligned with embedding rows
  model.joblib         - sklearn LogisticRegression on embeddings
  calibration.json     - threshold at --target-fpr on validation

Requires: torch, transformers, scikit-learn, joblib.
Default model facebook/esm2_t6_8M_UR50D runs on CPU; use t12_35M or t30_150M
with a GPU for stronger baselines.
"""

import argparse
import json
import os

import joblib
import numpy as np
from sklearn.linear_model import LogisticRegression

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
    return ids, [seqs[a] for a in ids], np.array([labels[a] for a in ids])

def embed(ids, seqs, model_name, batch_size, device, cache_prefix):
    """Mean-pooled ESM-2 embeddings, cached to <cache_prefix>.npy + ids."""
    npy, idfile = cache_prefix + ".npy", cache_prefix.replace("embeddings_", "ids_") + ".txt"
    if os.path.exists(npy) and os.path.exists(idfile):
        cached_ids = open(idfile).read().split()
        if cached_ids == ids:
            return np.load(npy)
    import torch
    from transformers import AutoModel, AutoTokenizer
    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(device).eval()
    feats = []
    with torch.no_grad():
        for i in range(0, len(seqs), batch_size):
            chunk = [s[:1022] for s in seqs[i:i + batch_size]]
            enc = tok(chunk, return_tensors="pt", padding=True,
                      truncation=True).to(device)
            out = model(**enc).last_hidden_state
            mask = enc["attention_mask"].unsqueeze(-1)
            pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
            feats.append(pooled.cpu().numpy().astype(np.float32))
            if (i // batch_size) % 20 == 0:
                print(f"  embedded {i}/{len(seqs)}")
    E = np.vstack(feats)
    np.save(npy, E)
    with open(idfile, "w") as fh:
        fh.write("\n".join(ids))
    return E

def threshold_at_fpr(y_true, scores, target_fpr):
    order = np.argsort(-scores)
    y, s = y_true[order], scores[order]
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
    for split in ("train", "val"):
        p.add_argument(f"--{split}-fasta", required=True)
        p.add_argument(f"--{split}-labels", required=True)
    p.add_argument("--test-fasta", default=None,
                   help="Optional: also embed the test split for evaluate.py.")
    p.add_argument("--test-labels", default=None)
    p.add_argument("-o", "--model-dir", required=True)
    p.add_argument("--model-name", default="facebook/esm2_t6_8M_UR50D")
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--device", default="cpu")
    p.add_argument("--C", type=float, default=1.0)
    p.add_argument("--target-fpr", type=float, default=0.01)
    args = p.parse_args()

    os.makedirs(args.model_dir, exist_ok=True)
    splits = {"train": (args.train_fasta, args.train_labels),
              "val": (args.val_fasta, args.val_labels)}
    if args.test_fasta and args.test_labels:
        splits["test"] = (args.test_fasta, args.test_labels)

    data = {}
    for name, (fa, lab) in splits.items():
        ids, seqs, y = load_split(fa, lab)
        print(f"Embedding {name}: {len(ids)} sequences")
        E = embed(ids, seqs, args.model_name, args.batch_size, args.device,
                  os.path.join(args.model_dir, f"embeddings_{name}"))
        data[name] = (E, y)

    Etr, ytr = data["train"]
    Eva, yva = data["val"]
    clf = LogisticRegression(C=args.C, class_weight="balanced",
                             max_iter=2000, n_jobs=-1)
    clf.fit(Etr, ytr)
    val_scores = clf.predict_proba(Eva)[:, 1]
    thr, rec = threshold_at_fpr(yva, val_scores, args.target_fpr)

    joblib.dump(clf, os.path.join(args.model_dir, "model.joblib"))
    with open(os.path.join(args.model_dir, "calibration.json"), "w") as fh:
        json.dump({"model_type": "esm2_embedding_logreg",
                   "backbone": args.model_name, "target_fpr": args.target_fpr,
                   "threshold": thr, "val_recall_at_threshold": rec},
                  fh, indent=2)
    print(f"Saved. Val threshold {thr:.4f} @ FPR<={args.target_fpr} "
          f"(val recall {rec:.3f}).")

if __name__ == "__main__":
    main()
