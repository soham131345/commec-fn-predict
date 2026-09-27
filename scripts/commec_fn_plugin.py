"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""
Drop-in function-prediction plugin for Commec (Common Mechanism) screening.

Commec's pipeline (commec/screen.py) runs fixed steps — biorisk HMM scan,
regulated-pathogen protein/nucleotide search, low-concern clearing — and has
no third-party plugin hook. This tool therefore integrates NON-INVASIVELY:

  1. Takes the same input FASTA that was screened, six-frame-translates it,
     and scores every ORF with a trained function-prediction model
     (k-mer or ESM-embedding classifier from train_*.py).
  2. Reads commec's structured output (<prefix>.output.json) and appends a
     "function_prediction" block per query with calibrated scores, ORF
     coordinates, and a flag/warn recommendation that commec's homology
     steps cannot produce for novel or heavily mutated sequences.
  3. Never alters commec's own hits or statuses — downstream tooling can
     adopt or ignore the extra block. A combined recommendation is written
     as "combined_recommendation" = worst of commec status vs. this layer.

Usage:
  python commec_fn_plugin.py --input-fasta order.fasta \
      --commec-json order.output.json --model-dir models/esm \
      -o order.augmented.json --csv order.function_scores.csv

Model dir must contain model.joblib and calibration.json. For ESM models,
embeddings are computed with the backbone named in calibration.json
(requires torch + transformers); k-mer models need only scikit-learn.
"""

import argparse
import csv
import json
import os
import sys

import joblib
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fn_common

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def read_fasta(path):
    seqs, header, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(buf)))
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line.upper().replace("U", "T"))
    if header is not None:
        seqs.append((header, "".join(buf)))
    return seqs

def six_frame_orfs(nt, min_aa=25):
    """Yield (frame, nt_start, nt_end, protein) for all ORFs >= min_aa."""
    rev = nt.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    for strand, s in (("+", nt), ("-", rev)):
        for frame in range(3):
            prot, start = [], None
            for i in range(frame, len(s) - 2, 3):
                aa = CODON_TABLE.get(s[i:i + 3], "X")
                if aa == "*" or i + 3 > len(s) - frame:
                    if prot and len(prot) >= min_aa and "X" not in "".join(prot):
                        yield (strand, frame, start, i, "".join(prot))
                    prot, start = [], None
                    continue
                if start is None:
                    start = i
                prot.append(aa)
            if prot and len(prot) >= min_aa and "X" not in "".join(prot):
                yield (strand, frame, start, len(s), "".join(prot))

def score_proteins(proteins, model_dir):
    calib = json.load(open(f"{model_dir}/calibration.json"))
    model = joblib.load(f"{model_dir}/model.joblib")
    if calib.get("model_type") == "esm2_embedding_logreg":
        import torch
        from transformers import AutoModel, AutoTokenizer
        tok = AutoTokenizer.from_pretrained(calib["backbone"])
        lm = AutoModel.from_pretrained(calib["backbone"]).eval()
        feats = []
        with torch.no_grad():
            for i in range(0, len(proteins), 16):
                chunk = [s[:1022] for s in proteins[i:i + 16]]
                enc = tok(chunk, return_tensors="pt", padding=True, truncation=True)
                out = lm(**enc).last_hidden_state
                mask = enc["attention_mask"].unsqueeze(-1)
                feats.append(((out * mask).sum(1) / mask.sum(1).clamp(min=1)).numpy())
        X = np.vstack(feats)
    else:
        X = proteins
    return model.predict_proba(X)[:, 1], float(calib["threshold"])

SEVERITY = {"ERROR": 4, "FLAG": 3, "WARN": 2, "PASS": 1, "SKIP": 0,
            "SKIP_SHORT": 0, "SKIP_LONG": 0, "PASS_SKIP_TX": 1, "CLEARED": 0}

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-fasta", required=True)
    p.add_argument("--commec-json", required=True,
                   help="commec screen output (<prefix>.output.json)")
    p.add_argument("--model-dir", required=True)
    p.add_argument("-o", "--out", required=True, help="augmented JSON path")
    p.add_argument("--csv", default=None, help="optional per-ORF score CSV")
    p.add_argument("--min-orf-aa", type=int, default=25)
    args = p.parse_args()

    data = json.load(open(args.commec_json))
    queries = data.get("queries", {})
    fasta = dict(read_fasta(args.input_fasta))

    rows = []
    for qname, qdata in queries.items():

        nt = fasta.get(qname) or fasta.get(qdata.get("query name", ""))
        if not nt:
            for fname, fseq in fasta.items():
                base = fname.split("|")[0].split()[0]
                if qname == base or fname.startswith(qname) or qname.startswith(base):
                    nt = fseq
                    break
        fp_block = {"status": "SKIPPED", "reason": "sequence not found in fasta",
                    "orfs": []}
        if nt:
            orfs = list(six_frame_orfs(nt, args.min_orf_aa))
            if orfs:
                scores, thr = score_proteins([o[4] for o in orfs], args.model_dir)
                orf_rows = []
                for (strand, frame, s, e, prot), sc in zip(orfs, scores):
                    orf_rows.append({
                        "strand": strand, "frame": frame,
                        "nt_start": s, "nt_end": e,
                        "aa_length": len(prot),
                        "toxin_virulence_score": round(float(sc), 4),
                        "above_threshold": bool(sc >= thr),
                    })
                    rows.append([qname, strand, frame, s, e, len(prot),
                                 f"{sc:.4f}", sc >= thr])
                best = max(float(sc) for sc in scores)
                status = "FLAG" if best >= thr else "PASS"
                fp_block = {"status": status,
                            "max_score": round(best, 4),
                            "threshold": thr,
                            "n_orfs_scored": len(orfs),
                            "orfs": orf_rows}
            else:
                fp_block = {"status": "PASS", "reason": "no ORFs >= min length",
                            "orfs": []}
        qdata["function_prediction"] = fp_block

        commec_status = str(qdata.get("status", {}).get("screen_status", "PASS")).upper()
        qdata["combined_recommendation"] = (
            fp_block["status"]
            if SEVERITY.get(fp_block["status"], 0) > SEVERITY.get(commec_status, 0)
            else commec_status)

    data.setdefault("plugin_info", {})["function_prediction"] = {
        "plugin": "commec-fn-predict",
        "model_dir": args.model_dir,
        "note": "Homology-independent ML scores appended per query; commec "
                "statuses unchanged. combined_recommendation merges both layers.",
    }

    with open(args.out, "w") as fh:
        json.dump(data, fh, indent=2)
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["query", "strand", "frame", "nt_start", "nt_end",
                        "aa_length", "score", "above_threshold"])
            w.writerows(rows)
    n_flag = sum(1 for q in queries.values()
                 if q.get("function_prediction", {}).get("status") == "FLAG")
    print(f"Wrote {args.out}: {n_flag}/{len(queries)} queries flagged by "
          f"function prediction.")

if __name__ == "__main__":
    main()
