[![DOI](https://zenodo.org/badge/1388954909.svg)](https://doi.org/10.5281/zenodo.22999841)

# commec-fn-predict

Homology-independent toxin/virulence function prediction for DNA synthesis screening, packaged as a non-invasive plugin for the Common Mechanism (commec).

**Author:** Soham Bhole, independent researcher
**Preprint:** In works
**License:** MIT (code and documentation)

---

## 1. What this project does

Current DNA synthesis screening asks: "does this sequence resemble a known hazard?"
Screening tools such as the Common Mechanism (IBBIS) and SecureDNA compare incoming
orders against databases of known dangerous sequences using BLAST and profile HMMs.
This works only when the sequence looks like something already catalogued.

This project asks a different question: "does this sequence encode a protein that
does something hazardous?" It scores DNA orders directly by machine learning models
trained on public toxin and virulence proteins, so it can flag:

- toxins whose sequence has drifted beyond database recognition,
- synonymously recoded DNA (identical protein, completely different DNA letters),
- novel proteins with toxin-like function and no database neighbor.

It ships as a plugin that augments commec's output JSON without modifying commec
itself, so screening providers can adopt it without forking anything.

## 2. Headline results

Benchmark: cluster-disjoint splits at 30 percent identity (no test sequence shares
a cluster with any training sequence), test set of 4,003 sequences
(3,267 positives, 736 negatives), thresholds frozen on validation at FPR <= 1 percent.

| Metric | k-mer (4-mer) + logistic regression | ESM-2 (8M) embedding + logistic regression |
|---|---|---|
| AUROC | 0.9121 | 0.9328 |
| AUPRC | 0.9785 | 0.9823 |
| Recall at FPR <= 1e-2 | 0.4778 | 0.4463 |
| Recall at FPR <= 1e-3 | 0.0260 | 0.0860 |
| Brier score | 0.1369 | 0.1059 |
| ECE (10 bins) | 0.2033 | 0.1323 |
| Test recall at frozen threshold | 0.3294 | 0.3642 |
| Test FPR at frozen threshold | 0.0041 | 0.0068 |

Drift experiment: each test positive was randomly mutated and scored at the frozen
threshold. The k-mer model collapses to noise by 40 percent identity. The ESM-2
model still recalls 7 percent of hazards at 20 percent identity, deep inside the
regime where homology search fails.

In-situ probe: 40 held-out toxins in fresh DNA — commec flagged 35/40 (87.5%);
the ESM-2 layer ranked all 4 of commec's misses in the 0.69–0.94 band and the
combined stack flagged 39/39 evaluable queries. See manuscript Section 3.3.

| Mean identity to original | k-mer recall | ESM-2 recall |
|---|---|---|
| 100% | 0.329 | 0.364 |
| 80% | 0.086 | 0.256 |
| 60% | 0.013 | 0.160 |
| 40% | 0.002 | 0.103 |
| 20% | 0.001 | 0.071 |

See `benchmark_report.md`, `drift_kmer/drift_curve.png`, and
`drift_esm/drift_curve.png` for the full tables and figure.

## 3. Repository layout

```
scripts/                  all pipeline code
  fetch_training_data.py      downloads labeled proteins from the UniProt REST API
  make_homology_splits.py     cluster-disjoint train/val/test splits at 30% identity
  fn_common.py                shared k-mer vectorizer (needed to load saved models)
  train_kmer_baseline.py      trains the 4-mer logistic regression baseline
  train_esm_classifier.py     trains the frozen ESM-2 embedding classifier
  evaluate.py                 AUROC, AUPRC, recall at fixed FPR, Brier, ECE
  drift_experiment.py         random-substitution robustness curves
  make_demo.py                builds the original / recoded / drifted demo sequences
  commec_fn_plugin.py         the commec plugin (six-frame ORF scoring, JSON output)
  subsample_fasta.py          utility to downsample a FASTA file
run_pipeline.bat          one-click Windows runner for the full pipeline
run_pipeline_part2.bat    resume script (skips download if data exists)
benchmark_report.md       the full benchmark table
drift_kmer/               k-mer drift curve (PNG + CSV)
drift_esm/                ESM-2 drift curve (PNG + CSV)
demo/                     demonstration sequences and plugin outputs
manuscript.md, manuscript.pdf   the preprint
```

Large files are excluded on purpose: raw UniProt FASTAs (re-downloaded by the
pipeline) and trained model weights (retrained by the pipeline). `data/README.txt`
explains how to regenerate them.

## 4. Requirements

- Windows 10/11 with Miniconda (the .bat files), or any Linux/macOS with Python 3.10+
- 8 GB RAM minimum, 16+ GB recommended for the ESM-2 step
- About 5 GB free disk space (databases, embeddings cache)
- No GPU required. ESM-2 embedding runs on CPU (slow but fine, a few hours once)

Python packages (installed automatically by the .bat files):
`requests scikit-learn scipy joblib numpy matplotlib edlib torch transformers`

## 5. Quickstart: reproduce everything

Windows:

1. Install Miniconda from https://docs.conda.io/en/latest/miniconda.html
   (default options are fine).
2. Download or clone this repository.
3. Double-click `run_pipeline.bat`. It will:
   - find conda and create the environment,
   - download positives (UniProt Tox-Prot and Virulence keywords) and negatives,
   - build cluster-disjoint splits,
   - train both models,
   - write `benchmark_report.md`,
   - run both drift experiments,
   - build the demo sequences in `demo/`.
4. Total time on a typical laptop: several hours, mostly the ESM-2 embedding step.
   You can close the lid problem-free only after it finishes; keep it plugged in.

Linux/macOS: run the scripts in `scripts/` in the order listed in section 3
(fetch, splits, train k-mer, train ESM, evaluate, drift, demo). Each script has
`--help`.

## 6. Using the plugin on your own DNA order

The plugin mirrors a commec run: it takes the same input FASTA and commec's
output JSON, and writes an augmented copy with ML verdicts added. Commec's own
hits and statuses are never modified.

```
python scripts/commec_fn_plugin.py \
    --input-fasta order.fasta \
    --commec-json order.output.json \
    --model-dir models/esm \
    -o order.augmented.json \
    --csv order.function_scores.csv
```

What it does:

1. Translates each query in all six reading frames and extracts ORFs of 25 aa or more.
2. Scores every ORF with the trained model (`models/kmer` also works, no torch needed).
3. Appends a `function_prediction` block per query: per-ORF scores and coordinates,
   the calibrated threshold, and a FLAG or PASS status.
4. Adds `combined_recommendation`: the more severe of commec's verdict and this
   layer's verdict.

Output columns in the CSV: query, strand, frame, nucleotide coordinates, ORF length,
score (0 to 1), and whether it exceeds the frozen threshold.

Known limitation (reported in the preprint): very short ORFs (25 to 70 aa) lie far
outside the training distribution and can produce extreme scores. Query-level
max-score aggregation is conservative (it flags more, not less), but
ORF-length-aware aggregation is recommended before operational use.

## 7. The demo

`demo/` contains a real reviewed venom toxin as DNA in three forms:

- `original.fasta`: the toxin reverse-translated with random synonymous codons,
- `recoded.fasta`: the identical protein with a different random codon choice
  (different DNA, same protein),
- `drifted.fasta`: the protein mutated to 58.4 percent identity, then
  reverse-translated,
- `mock_commec.output.json`: a commec-format report marking all three PASS
  (no database resemblance, by construction),
- `*.augmented.json` and `*.scores.csv`: this layer's output on each construct.

Result: the k-mer model flags original and recoded with a bit-identical score
(synonymous recoding is a non-event for a protein-function model), while commec
format reports PASS for all three. The drifted construct escapes on its true
reading frame, an honest, measured frontier that the drift curves quantify.

To run the real commec screen on these files (Linux or WSL):

```
conda create -y -n commec -c conda-forge -c bioconda commec
conda activate commec
commec setup -d commec_dbs
commec screen demo/combined.fasta -d commec_dbs -o demo_out
python scripts/commec_fn_plugin.py --input-fasta demo/combined.fasta \
    --commec-json demo_out/combined.output.json \
    --model-dir models/kmer -o demo/combined.augmented.json
```

## 8. Data sources and ethics

- All training data is public and already annotated: reviewed UniProtKB/Swiss-Prot
  entries with the keywords Toxin (KW-0800, the Tox-Prot program) or Virulence
  (KW-0843) as positives, and keyword-excluded reviewed entries as negatives.
- This is defensive work. It releases classifiers, not generators. Drifted
  sequences are used only to measure screening robustness. No tooling for
  optimizing sequences against screens is included.
- The intended adopters are DNA synthesis providers and screening maintainers
  (IBBIS, SecureDNA), alongside whom findings of this kind should be shared.

## 9. Citing this work

```
Bhole, S. (2026). Homology-independent function prediction as a drop-in layer
for DNA synthesis screening: calibrated k-mer and protein language model
baselines with a non-invasive plugin for the Common Mechanism.
bioRxiv. DOI to be added after posting.
```

## 10. Contact

Soham Bhole, independent researcher.
Issues and discussion: please use the GitHub Issues tab on this repository.
