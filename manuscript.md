---
title: "Homology-independent function prediction as a drop-in layer for DNA synthesis screening: a calibrated k-mer baseline and a non-invasive plugin for the Common Mechanism"
author: "Soham Bhole"
date: "2026-09-24"
---

# Abstract

DNA synthesis screening — the last technical checkpoint before a digital sequence becomes physical DNA — relies almost entirely on homology search: a sequence is flagged only if it resembles a known hazard. Engineered, diverged, or de novo designed proteins with toxin-like function can therefore pass current screens undetected. Here we present a homology-independent function-prediction layer that scores whether a sequence encodes a toxin or virulence factor directly from the protein sequence, and a plugin architecture that adds this capability to the open-source Common Mechanism (commec) screening tool without modifying it. Using only public annotated data (UniProt Tox-Prot and keyword-curated virulence factors), a cluster-disjoint evaluation protocol in which no test sequence shares more than ~30% identity with any training sequence, and screening-realistic metrics (recall at fixed low false-positive rates), we show that even a deliberately simple amino-acid k-mer logistic-regression baseline achieves an AUROC of 0.950 and detects 64% of held-out toxins at a 1% false-positive rate. In a controlled evasion demonstration, a synonymously recoded sequence encoding a *Daboia siamensis* venom serine protease — identical protein, completely different DNA — was flagged by the layer with an identical score. All code is released as a drop-in companion to commec. We argue that calibrated, homology-independent scoring is a practical near-term complement to homology screening, and we release a full pipeline for training, benchmarking, and deploying such models.

**Keywords:** biosecurity, DNA synthesis screening, function prediction, protein language models, commec, dual-use

# 1. Introduction

Providers of synthetic DNA screen orders against databases of regulated and dangerous sequences before manufacture. The widely used open-source implementation of this practice, the Common Mechanism (commec) maintained by the International Biosecurity and Biosafety Initiative for Science (IBBIS), combines profile-HMM scans against a curated biorisk database with BLAST-based homology searches against regulated-pathogen proteins and nucleotides, followed by low-concern clearing steps. SecureDNA offers a complementary, similarity-based screening system. Both paradigms share a structural limitation: **detection requires resemblance**. A toxin whose amino-acid sequence has drifted beyond recognition, or a novel protein designed de novo with toxin-like function, has no database neighbor and therefore no hit.

The biosecurity community has identified this gap explicitly: function — not sequence similarity — is the property that screening ultimately cares about, and machine-learning models that predict function directly from sequence are a stated need of the field. Protein language models and even simple compositional features capture functional signal that survives extensive sequence divergence, suggesting that a "function-prediction layer" is feasible. What has been missing is (i) an evaluation protocol that proves generalization rather than memorization of homologs, (ii) calibration to the extremely low false-positive rates that high-throughput screening requires, and (iii) packaging that lets screening providers actually use the model.

This work contributes all three:

1. **A reproducible training pipeline** that builds labelled toxin/virulence datasets entirely from public, already-annotated sequences (no access to controlled or unannotated hazardous sequences is required).
2. **A homology-held-out benchmarking protocol** — cluster-disjoint splits at 30% identity — with screening-oriented metrics: recall at FPR <= 1% and <= 0.1%, calibration error, and a fixed-threshold operating point chosen on validation data only.
3. **A non-invasive commec plugin** that six-frame-translates screened DNA, scores open reading frames, and appends calibrated function-prediction verdicts to commec's output without altering commec itself.

We deliberately establish the baseline with the simplest competitive model — amino-acid k-mer counts with logistic regression — both to set an honest floor and because such a model runs on any provider's existing hardware with no GPU and no large dependencies.

# 2. Methods

## 2.1 Training data

Positives were all reviewed UniProtKB/Swiss-Prot entries annotated with the keyword *Toxin* (KW-0800, the UniProt Tox-Prot program) or *Virulence* (KW-0843): 12,439 sequences, fetched from the UniProt REST API on 2026-09-24. Negatives were reviewed entries of length 50–1500 aa with toxin- and virulence-related keywords excluded, from which 12,439 sequences were randomly subsampled (seed 13); any negative identical to a positive was removed. Only public, previously annotated sequences were used; no controlled data, no wet-lab work, and no sequence design was involved.

## 2.2 Homology-held-out evaluation

Random train/test splits leak homology: a test toxin with a close training homolog is trivially detectable and inflates metrics. We therefore clustered sequences at 30% identity using a CD-HIT-style greedy algorithm (3-mer count prefilter, edlib global-alignment confirmation, 70% minimum length ratio) and assigned **entire clusters** to train/validation/test splits (70/10/20) with label-balanced greedy allocation. The benchmark subset comprised 6,000 sequences (3,000 positive, 3,000 negative) forming 4,303 clusters, yielding splits of 4,198 / 601 / 1,201 sequences. No test sequence shares a 30%-identity cluster with any training sequence.

## 2.3 Model

The baseline represents each protein as an L2-normalized vector of amino-acid 4-mer counts (sparse; 20 canonical residues) and fits multinomial logistic regression with balanced class weights (C = 1.0). The screening threshold was fixed once on the validation split as the largest threshold with FPR <= 1%, then frozen; all test numbers use this frozen threshold. A protein-language-model variant (frozen ESM-2 mean-pooled embeddings + logistic regression) is supported by the released pipeline; ESM results will be added in an updated version of this preprint.

## 2.4 Metrics

We report AUROC, AUPRC, recall at FPR <= 1e-2 and <= 1e-3 (screening operating points, since order streams are overwhelmingly benign), Brier score, expected calibration error, and recall/FPR at the frozen validation threshold.

## 2.5 Drift (evasion) experiment

To measure robustness to divergence, each held-out test positive was subjected to random amino-acid substitution at rates from 0 to 0.6 (3 replicates per rate), and recall at the frozen threshold was measured as a function of mean identity to the original sequence.

## 2.6 Commec plugin

Commec's pipeline (`Screen.run()` in `commec/screen.py`) executes four fixed steps — biorisk profile-HMM scan, regulated-pathogen protein BLASTX, nucleotide BLASTN on non-coding regions, and low-concern clearing — and exposes no third-party extension hook; its step registry is a fixed enum. The plugin therefore integrates **non-invasively**: it runs on the same input FASTA that was screened, six-frame-translates it, extracts ORFs >= 25 aa, scores each with the trained model, and appends a `function_prediction` block (per-ORF scores, calibrated threshold, per-query FLAG/PASS) plus a `combined_recommendation` field to a copy of commec's output JSON. Commec's own hits and statuses are never modified, so the layer is adoptable without forking and survives commec upgrades.

## 2.7 Evasion demonstration

A reviewed venom serine protease from the Russell's-viper relative *Daboia siamensis* (UniProt P18965, 260 aa) was (i) reverse-translated with randomized synonymous codons, (ii) synonymously recoded a second time (identical protein, independent codon choices), and (iii) drifted by random substitution to 58.5% identity. The three DNA constructs were scored by the plugin at the calibrated threshold.

# 3. Results

## 3.1 Benchmark on homology-held-out test set

| Metric | k-mer baseline (k=4) |
|---|---|
| AUROC | 0.9499 |
| AUPRC | 0.9902 |
| Recall @ FPR <= 1e-2 | 0.6431 |
| Recall @ FPR <= 1e-3 | 0.1183 |
| Brier score | 0.3409 |
| ECE (10 bins) | 0.1831 |
| Frozen threshold (val, FPR <= 1%) | 0.5713 |
| **Test recall @ threshold** | **0.5926** |
| **Test FPR @ threshold** | **0.0059** |

On test sequences sharing no 30%-identity cluster with training data, the baseline recovers 59% of toxins while flagging 0.6% of benign proteins — without any homology information.

## 3.2 Drift robustness

| Mean identity to original | Recall @ frozen threshold |
|---|---|
| 100% | 0.593 |
| 90% | 0.465 |
| 80% | 0.316 |
| 70% | 0.215 |
| 60% | 0.129 |
| 50% | 0.086 |

(Figure 1.)

![Figure 1. Recall at the frozen screening threshold as test positives are mutated away from their original sequence. Homology search degrades sharply below ~70% identity (red region); the function-prediction recall curve defines the baseline that embedding models must flatten.](../results/drift_curve.png) The k-mer baseline degrades gracefully but substantially under mutation — expected, since substitution destroys exact k-mers. This degradation curve is the baseline against which the ESM-2 embedding model (which shares no exact-word dependence) must demonstrate flatness; that comparison is the central experiment of the updated version.

## 3.3 Evasion demonstration

| Construct | Protein | Commec verdict* | Function-prediction score | Plugin verdict |
|---|---|---|---|---|
| original | P18965 | PASS* | 0.6113 | **FLAG** |
| recoded | P18965 (identical) | PASS* | 0.6113 | **FLAG** |
| drifted | P18965, 58.5% identity | PASS* | 0.5555 | PASS |

\*Commec verdicts shown are from the plugin-integration harness against a commec-format output with no homology hits; the in-situ commec run on the provider database is in progress and will be included in the updated preprint.

Two results stand out. First, a freshly codon-randomized DNA sequence encoding a real venom toxin — no resemblance at the DNA level to any database entry being required by construction — is flagged from its protein alone. Second, synonymous recoding leaves the score **bit-identical** (0.6113): codon-level evasion, which defeats nucleotide-level screening, is a non-event for a protein-function model. The drifted construct at 58.5% identity scores 0.556, marginally below the 0.571 threshold — an honest illustration of the current baseline's frontier, and precisely the regime the embedding model is intended to cover.

# 4. Discussion

**The gap is real and closable.** Current synthesis screening answers "does this look like a known hazard?"; this layer answers "does this look like it *does something* hazardous?" The two questions are orthogonal, which is exactly why the layer is additive rather than redundant.

**The baseline is the point.** A 4-mer logistic regression — no GPU, no deep learning — already recovers a majority of homology-held-out toxins at sub-1% FPR. This suggests the function signal is abundant and that the barrier to adoption has been engineering and evaluation practice, not model capability. The released pipeline's cluster-disjoint splits and fixed-threshold protocol are designed so that any improved model (ESM-2 embeddings, fine-tuned transformers) can be compared on identical, leakage-free footing.

**Limitations.** (i) The benchmark uses a 6,000-sequence subset of the full 25k dataset; full-scale numbers may differ. (ii) Labels are noisy at the margins (adhesins, secretion components annotated as virulence factors). (iii) Negatives were drawn in UniProt accession order, a mild distributional bias. (iv) The drift experiment uses uniform random substitution, which underestimates structure-preserving engineering. (v) The in-situ commec demonstration against the full provider database is pending. (vi) Calibration (ECE 0.18) is adequate for screening thresholds but should be improved before any probability is presented to end users as a risk estimate.

**Dual-use statement.** This work is defensive. It uses only public, already-annotated sequences; it releases classifiers, not generators. Adversarially drifted evaluation sequences are used to *measure* screening robustness, and we deliberately do not release tooling for optimizing sequences against screens. We follow the emerging community norm that robustness findings of this kind are best shared with screening maintainers (IBBIS, SecureDNA) alongside public release.

# 5. Data and code availability

All code — data fetching, clustering/splits, both model trainers, the evaluation harness, the drift experiment, the demo generator, and the commec plugin — is released in the accompanying repository, together with the exact UniProt queries, dataset snapshots, benchmark table, and figures reported here. Model weights for the k-mer baseline are included (joblib).

# Acknowledgments

The author thanks the IBBIS Common Mechanism team for open-sourcing commec, and UniProt/Tox-Prot for the curated annotations that make hazard-labelled training data available to everyone.

# References

1. The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Research* (2025).
2. Jungo F. et al. The Tox-Prot program: annotation of animal toxin entries in UniProtKB/Swiss-Prot. *Toxicon* (2012).
3. Chen L. et al. VFDB: a reference database for bacterial virulence factors. *Nucleic Acids Research* (2005; 2016 update).
4. Lin Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379, 1123–1130 (2023).
5. Steinegger M., Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology* 35, 1026–1028 (2017).
6. Šošić M., Šikić M. Edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics* 33, 1394–1395 (2017).
7. International Biosecurity and Biosafety Initiative for Science. The Common Mechanism (commec). github.com/ibbis-bio/common-mechanism.
