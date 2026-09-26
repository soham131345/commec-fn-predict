# Function-prediction benchmark (homology-held-out test set)

Test set: 4003 sequences (3267 positives).

| model | AUROC | AUPRC | recall@FPR1e-2 | recall@FPR1e-3 | Brier | ECE | test_recall@thr | test_FPR@thr |
|---|---|---|---|---|---|---|---|---|
| kmer | 0.9121 | 0.9785 | 0.4778 | 0.0260 | 0.1369 | 0.2033 | 0.3294 | 0.0041 |
| esm | 0.9328 | 0.9823 | 0.4463 | 0.0860 | 0.1059 | 0.1323 | 0.3642 | 0.0068 |
