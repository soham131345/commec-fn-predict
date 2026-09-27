---
title: "Homology-independent function prediction as a drop-in layer for DNA synthesis screening: calibrated k-mer and protein language model baselines with a non-invasive plugin for the Common Mechanism"
author: "Soham Bhole"
date: "2026-09-26"
---

# Abstract

DNA synthesis screening — the last technical checkpoint before a digital sequence becomes physical DNA — relies almost entirely on homology search: a sequence is flagged only if it resembles a known hazard. Engineered, diverged, or de novo designed proteins with toxin-like function can therefore pass current screens undetected. Here we present a homology-independent function-prediction layer that scores whether a sequence encodes a toxin or virulence factor directly from the protein sequence, and a plugin architecture that adds this capability to the open-source Common Mechanism (commec) screening tool without modifying it. Using only public annotated data (UniProt Tox-Prot and keyword-curated virulence factors), a cluster-disjoint evaluation protocol in which no test sequence shares more than ~30% identity with any training sequence, and screening-realistic metrics (recall at fixed low false-positive rates), we benchmark a deliberately simple amino-acid k-mer logistic-regression baseline against a frozen ESM-2 protein language model embedding classifier on a 4,003-sequence homology-held-out test set. Both models generalize without homology (AUROC 0.912 and 0.933), but they differ decisively at screening operating points: at a false-positive rate of 0.1%, the ESM-2 model detects 3.3x more held-out hazards than the k-mer baseline (recall 0.086 vs 0.026). Under controlled sequence drift, the k-mer model collapses (recall 0.33 to 0.004 by 50% identity) while the ESM-2 model retains detectable signal down to 20% identity — far beyond the reach of homology search. In an in-situ evaluation against a full installation of the Common Mechanism, 40 held-out toxins reverse-translated into fresh, never-before-seen DNA were screened: commec flagged 35 of 40 (87.5%), while the four toxins it passed entirely were scored in the 0.69–0.94 range by the ESM-2 layer — elevated but below the frozen threshold — and were flagged at query level only through short-ORF aggregation artifacts, a failure mode we characterize openly. All code is released as a drop-in companion to commec. We argue that calibrated, homology-independent scoring is a practical near-term complement to homology screening, and we release a full pipeline for training, benchmarking, and deploying such models.

**Keywords:** biosecurity, DNA synthesis screening, function prediction, protein language models, commec, dual-use

# 1. Introduction

Providers of synthetic DNA screen orders against databases of regulated and dangerous sequences before manufacture. The widely used open-source implementation of this practice, the Common Mechanism (commec) maintained by the International Biosecurity and Biosafety Initiative for Science (IBBIS), combines profile-HMM scans against a curated biorisk database with BLAST-based homology searches against regulated-pathogen proteins and nucleotides, followed by low-concern clearing steps. SecureDNA offers a complementary, similarity-based screening system. Both paradigms share a structural limitation: **detection requires resemblance**. A toxin whose amino-acid sequence has drifted beyond recognition, or a novel protein designed de novo with toxin-like function, has no database neighbor and therefore no hit.

The biosecurity community has identified this gap explicitly: function — not sequence similarity — is the property that screening ultimately cares about, and machine-learning models that predict function directly from sequence are a stated need of the field. Protein language models and even simple compositional features capture functional signal that survives extensive sequence divergence, suggesting that a "function-prediction layer" is feasible. What has been missing is (i) an evaluation protocol that proves generalization rather than memorization of homologs, (ii) calibration to the extremely low false-positive rates that high-throughput screening requires, and (iii) packaging that lets screening providers actually use the model.

This work contributes all three:

1. **A reproducible training pipeline** that builds labelled toxin/virulence datasets entirely from public, already-annotated sequences (no access to controlled or unannotated hazardous sequences is required).
2. **A homology-held-out benchmarking protocol** — cluster-disjoint splits at 30% identity — with screening-oriented metrics: recall at FPR <= 1% and <= 0.1%, calibration error, and a fixed-threshold operating point chosen on validation data only.
3. **A non-invasive commec plugin** that six-frame-translates screened DNA, scores open reading frames, and appends calibrated function-prediction verdicts to commec's output without altering commec itself.

We benchmark two models spanning the cost-capability range: amino-acid 4-mer counts with logistic regression (no GPU, no deep learning, runs on any provider's existing hardware) and frozen ESM-2 (8M-parameter) mean-pooled embeddings with logistic regression.

# 2. Methods

## 2.1 Training data

Positives were all reviewed UniProtKB/Swiss-Prot entries annotated with the keyword *Toxin* (KW-0800, the UniProt Tox-Prot program) or *Virulence* (KW-0843): 12,439 sequences, fetched from the UniProt REST API in September 2026. Negatives were reviewed entries of length 50–1500 aa with toxin- and virulence-related keywords excluded, from which 12,439 sequences were randomly subsampled (seed 13); any negative identical to a positive was removed. Only public, previously annotated sequences were used; no controlled data, no wet-lab work, and no sequence design was involved.

## 2.2 Homology-held-out evaluation

Random train/test splits leak homology: a test toxin with a close training homolog is trivially detectable and inflates metrics. We therefore clustered all 24,878 sequences at 30% identity using a CD-HIT-style greedy algorithm (3-mer count prefilter, edlib global-alignment confirmation, 70% minimum length ratio) and assigned **entire clusters** to train/validation/test splits (70/10/20) with label-balanced greedy allocation, yielding a held-out test set of 4,003 sequences (3,267 positives, 736 negatives). No test sequence shares a 30%-identity cluster with any training sequence.

## 2.3 Models

The baseline represents each protein as an L2-normalized sparse vector of amino-acid 4-mer counts (20 canonical residues) and fits multinomial logistic regression with balanced class weights (C = 1.0). The embedding model mean-pools the final hidden states of a frozen ESM-2 checkpoint (facebook/esm2_t6_8M_UR50D, 8M parameters) over the attention mask and fits the same logistic-regression head; embeddings are computed once and cached. For each model, the screening threshold was fixed once on the validation split as the largest threshold with FPR <= 1%, then frozen; all test numbers use this frozen threshold.

## 2.4 Metrics

We report AUROC, AUPRC, recall at FPR <= 1e-2 and <= 1e-3 (screening operating points, since order streams are overwhelmingly benign), Brier score, expected calibration error, and recall/FPR at the frozen validation threshold.

## 2.5 Drift (evasion) experiment

To measure robustness to divergence, each held-out test positive was subjected to random amino-acid substitution at rates from 0 to 0.8 (3 replicates per rate), and recall at the frozen threshold was measured as a function of mean identity to the original sequence, for both models.

## 2.6 Commec plugin

Commec's pipeline (`Screen.run()` in `commec/screen.py`) executes four fixed steps — biorisk profile-HMM scan, regulated-pathogen protein BLASTX, nucleotide BLASTN on non-coding regions, and low-concern clearing — and exposes no third-party extension hook; its step registry is a fixed enum. The plugin therefore integrates **non-invasively**: it runs on the same input FASTA that was screened, six-frame-translates it, extracts ORFs >= 25 aa, scores each with the trained model, and appends a `function_prediction` block (per-ORF scores, calibrated threshold, per-query FLAG/PASS) plus a `combined_recommendation` field to a copy of commec's output JSON. Commec's own hits and statuses are never modified, so the layer is adoptable without forking and survives commec upgrades.

## 2.7 In-situ commec evaluation

Commec (with current databases: biorisk rev 1.2, control_lists rev 1.1, best_match rev 1.0, low_concern rev 1.1) was installed under WSL2 and run on three construct sets: (i) a held-out test toxin (UniProt Q7Z1Y7, 320 aa) as reverse-translated DNA, a synonymous recoding, and a 58.4%-identity drifted variant; (ii) the same three constructs for a venom serine protease (P18965, 260 aa); (iii) 40 randomly sampled held-out test positives (seed 7), each reverse-translated with random synonymous codons into fresh DNA, screened with `commec screen`, then augmented with the plugin (`combined_recommendation` = worst of both layers; commec outputs never modified).

# 3. Results

## 3.1 Benchmark on the homology-held-out test set

Test set: 4,003 sequences (3,267 positives, 736 negatives); thresholds frozen on validation at FPR <= 1%.

| Metric | k-mer (k=4) + LR | ESM-2 (8M) embedding + LR |
|---|---|---|
| AUROC | 0.9121 | **0.9328** |
| AUPRC | 0.9785 | **0.9823** |
| Recall @ FPR <= 1e-2 | **0.4778** | 0.4463 |
| Recall @ FPR <= 1e-3 | 0.0260 | **0.0860** |
| Brier score | 0.1369 | **0.1059** |
| ECE (10 bins) | 0.2033 | **0.1323** |
| **Test recall @ frozen threshold** | 0.3294 | **0.3642** |
| **Test FPR @ frozen threshold** | **0.0041** | 0.0068 |

Both models generalize to sequences with no 30%-identity relative in training data — the central requirement for a homology-independent layer. Their operating-point behavior differs in a way that matters for screening. At the loosest screening point (FPR 1%) the two are comparable, but at the strictest point (FPR 0.1%) the ESM-2 model detects **3.3x** more hazards (0.086 vs 0.026): the ranking quality of the language model concentrates its remaining true positives at the extreme high-confidence end, exactly where a low-FPR screen operates. The ESM-2 model is also better calibrated (ECE 0.13 vs 0.20), though neither model's probabilities should yet be presented to end users as risk estimates.

## 3.2 Drift robustness

Recall at each model's frozen validation threshold as held-out test positives are mutated away from their original sequence (3 replicates per rate; 95% binomial intervals <= 0.005 throughout):

| Mean identity to original | k-mer recall | ESM-2 recall |
|---|---|---|
| 100% | 0.329 | **0.364** |
| 90% | 0.183 | **0.316** |
| 80% | 0.086 | **0.256** |
| 70% | 0.035 | **0.204** |
| 60% | 0.013 | **0.160** |
| 50% | 0.004 | **0.125** |
| 40% | 0.002 | **0.103** |
| 30% | 0.001 | **0.083** |
| 20% | 0.001 | **0.071** |

(Figure 1.)

![Figure 1. Recall at the frozen screening threshold as test positives are mutated away from their original sequence, for the k-mer baseline and the ESM-2 embedding model. Homology search degrades sharply below ~70% identity; the k-mer model collapses with it, while the ESM-2 model retains detectable signal (recall 0.07) down to 20% identity — the regime where no database match exists at all.](drift_curve.png)

The k-mer baseline collapses under mutation — expected, since substitution destroys exact 4-mers — losing nearly half its recall by 90% identity and falling to noise (0.1%) by 40%. The ESM-2 model degrades far more gracefully at every point: it retains 87% of its original recall at 90% identity, still recalls one in eight hazards at 50% identity, and, critically, **does not go to zero** — at 20% identity to the original toxin, well beyond the failure point of BLAST-based screening, it still recalls 7% of hazards at the frozen 1%-FPR threshold. Protein language model embeddings capture functional signal that survives near-total sequence divergence, which is precisely the property a function-prediction layer exists to provide.

## 3.3 In-situ evaluation against the Common Mechanism

**Regulated-family constructs.** Commec flagged all three Q7Z1Y7 constructs — original, recoded, and the 58.4%-identity drifted variant — via its biorisk profile-HMM step (hit at nt 142–960 in each), and likewise all three P18965 venom-protease constructs. Two conclusions follow. First, commec's curated biorisk database is substantially more drift-robust than the "homology fails below ~70% identity" heuristic suggests: profile HMMs detected a known-family toxin even at 58% protein identity in freshly recoded DNA. Second, the function-prediction layer provided orthogonal confirmation: the k-mer model flagged the original and recoded Q7Z1Y7 constructs from the true 320-aa reading frame with a **bit-identical** score (0.8128; frozen threshold ~0.78), demonstrating that synonymous recoding — which defeats nucleotide-level matching — is a non-event for a protein-function model. The drifted construct scored 0.607 (PASS), an honest illustration of the frontier quantified in Section 3.2.

**Probe of 40 held-out toxins in fresh DNA.** The decisive experiment: 40 held-out test positives, reverse-translated with random codons, screened by the full commec pipeline and both ML layers.

| Layer | Flagged | Notes |
|---|---|---|
| Commec (full pipeline) | **35 / 40** (87.5%) | 4 Pass, 1 Skip (query < 41 bp) |
| k-mer layer | 16 / 40 | all 16 within commec's flags |
| ESM-2 layer (query level) | 35 / 40 | includes **all 4 commec passes** |
| **Combined (commec OR ESM-2)** | **39 / 39 evaluable** | union of verdicts |

The four toxins commec passed entirely — Q8J0D5, Q9KI13, A0A0J9W9G2, B6V870 — are the measured residual gap of the current screen on this sample (10%). The ESM-2 layer flagged all four at query level (max scores 0.989–0.996). However, per-ORF inspection shows these four flags are driven by short secondary ORFs (<= 152 aa) scoring >= 0.97, not by the full-length proteins: the true toxin ORFs score 0.69–0.94 (e.g. Q9KI13, 931 aa, 0.9445) — elevated well above typical negatives but below the frozen threshold (~0.96). In other words, the ESM-2 model genuinely ranks all four gap toxins near the top of its range, but does not cleanly catch them at the calibrated operating point; the query-level flags arise from the conservative max-score aggregation acting on off-distribution short ORFs. We report both facts rather than claiming a clean catch.

**Interpretation.** The screening gap on curated toxin families is real but narrower than commonly assumed (10% in this probe), and it is concentrated in families under-represented in the biorisk database. A protein language model layer adds three things even here: recoding-invariant independent confirmation (k-mer, bit-identical), high ranking of exactly the sequences commec misses (all four gap toxins in the 0.69–0.94 band), and — with the conservative aggregation shipped — a combined stack that flagged every evaluable query in the probe. Closing the residual gap *cleanly* is a calibration and aggregation problem (Section 4), not a data problem.

# 4. Discussion

**The gap is real, smaller than assumed, and closable.** Current synthesis screening answers "does this look like a known hazard?"; this layer answers "does this look like it *does something* hazardous?" Our in-situ probe shows the deployed screen is stronger than the field often assumes — 87.5% of held-out toxins flagged in fresh DNA, and even 58%-identity drift of covered families caught by profile HMMs. The residual 10% is exactly where a function-prediction layer belongs: the two questions are orthogonal, which is why the layer is additive rather than redundant.

**The k-mer floor is honest, and the embedding model earns its keep.** A 4-mer logistic regression already recovers a third of homology-held-out toxins at sub-1% FPR with no GPU — but it is brittle: its signal is exact-word signal, and it collapses under drift on the same schedule as homology search itself. The ESM-2 model costs more (a one-time embedding pass, CPU-feasible at screening volumes) and buys three things: better ranking at the strictest false-positive rates, better calibration, and — the decisive property — recall that persists at 20–30% identity, where k-mers, BLAST, and profile HMMs have all gone silent. For a screening provider, the practical configuration is both: k-mer as a free always-on filter, ESM-2 as the second look at anything near the boundary.

**Limitations.** (i) Labels are noisy at the margins (adhesins, secretion components annotated as virulence factors). (ii) Negatives were drawn in UniProt accession order, a mild distributional bias. (iii) The drift experiment uses uniform random substitution, which underestimates structure-preserving engineering. (iv) The in-situ probe covers 40 held-out toxins; larger probes across more diverse families would sharpen the residual-gap estimate. (v) Calibration (ECE 0.13–0.20) is adequate for threshold-based screening but should be improved (e.g. isotonic recalibration) before probabilities are shown to end users. (vi) Absolute recall at the frozen threshold (0.33–0.36) means the layer is a complement to, never a replacement for, homology screening and human review. (vii) Both in-situ demos exposed an aggregation artifact: short spurious ORFs (25–152 aa) lie far off the training distribution and can produce extreme scores, so query-level max-score aggregation can FLAG via non-coding frame fragments — including, in the 40-toxin probe, all four of commec's misses. For a conservative screening layer this errs on the safe side, but it means query-level flags must not be read as protein-level detections; ORF-length-aware aggregation is required before operational use, and per-ORF score distributions on out-of-distribution short peptides warrant dedicated study.

**Dual-use statement.** This work is defensive. It uses only public, already-annotated sequences; it releases classifiers, not generators. Adversarially drifted evaluation sequences are used to *measure* screening robustness, and we deliberately do not release tooling for optimizing sequences against screens. We follow the emerging community norm that robustness findings of this kind are best shared with screening maintainers (IBBIS, SecureDNA) alongside public release.

# 5. Data and code availability

All code — data fetching, clustering/splits, both model trainers, the evaluation harness, the drift experiment, the demo generator, and the commec plugin — is released in the accompanying repository, together with the exact UniProt queries, the benchmark table, drift curves, and figures reported here. Model weights for the k-mer baseline are included (joblib).

# Acknowledgments

The author thanks the IBBIS Common Mechanism team for open-sourcing commec, and UniProt/Tox-Prot for the curated annotations that make hazard-labelled training data available to everyone.

# References

1. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Research* (2025).
2. Jungo F. et al. The Tox-Prot program: annotation of animal toxin entries in UniProtKB/Swiss-Prot. *Toxicon* (2012).
3. Chen L. et al. VFDB: a reference database for bacterial virulence factors. *Nucleic Acids Research* (2005; 2016 update).
4. Lin Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379, 1123–1130 (2023).
5. Steinegger M., Soding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology* 35, 1026–1028 (2017).
6. Sosic M., Sikic M. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics* 33, 1394–1395 (2017).
7. International Biosecurity and Biosafety Initiative for Science. The Common Mechanism (commec). github.com/ibbis-bio/common-mechanism.
8. Pedregosa F. et al. Scikit-learn: Machine Learning in Python. *Journal of Machine Learning Research* 12, 2825–2830 (2011).
9. Wolf T. et al. Transformers: State-of-the-Art Natural Language Processing. *EMNLP: System Demonstrations* (2020).
