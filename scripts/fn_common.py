"""
Copyright (c) 2026 Soham Bhole. All rights reserved.
"""

"""Shared, importable components for the function-prediction pipeline.

KmerVectorizer lives here (not in a training script) so that joblib-pickled
models unpickle correctly in evaluate.py and commec_fn_plugin.py — classes
defined in __main__ are not portable across processes.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

class KmerVectorizer(BaseEstimator, TransformerMixin):
    """k-mer count vectors over the 20 canonical amino acids, L2-normalized."""

    def __init__(self, k=4):
        self.k = k
        self.vocab = None

    def fit(self, X, y=None):
        vocab, seen = {}, set()
        for seq in X:
            for i in range(len(seq) - self.k + 1):
                kmer = seq[i:i + self.k]
                if all(c in ALPHABET for c in kmer) and kmer not in seen:
                    seen.add(kmer)
                    vocab[kmer] = len(vocab)
        self.vocab = vocab
        return self

    def transform(self, X):
        from scipy import sparse
        rows, cols, vals = [], [], []
        for r, seq in enumerate(X):
            for i in range(len(seq) - self.k + 1):
                j = self.vocab.get(seq[i:i + self.k])
                if j is not None:
                    rows.append(r)
                    cols.append(j)
                    vals.append(1.0)
        M = sparse.csr_matrix((vals, (rows, cols)),
                              shape=(len(X), len(self.vocab)),
                              dtype=np.float32)
        norms = np.asarray(M.multiply(M).sum(axis=1)).ravel() ** 0.5
        norms[norms == 0] = 1.0
        return M.multiply(1.0 / norms[:, None]).tocsr()
