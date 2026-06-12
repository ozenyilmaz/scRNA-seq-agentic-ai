# Multi-Criteria Therapeutic Target Prioritization from single-cell RNA-Seq
Using Agentic AI Knowledge-Harmonization and Shannon Entropy Weighting

[![Open In Colab](https://img.shields.io/badge/Colab-Run%20Pipeline-orange?style=for-the-badge&logo=google-colab)](https://colab.research.google.com/github/KULLANICI_ADIN/scRNA-seq-target-prioritization/blob/main/)

## Executive Summary
[cite_start]This repository contains a modular, production-grade single-cell RNA-sequencing (scRNA-seq) framework designed for drug target prioritization[cite: 301]. [cite_start]Engineered as an **Agentic AI Knowledge-Harmonization Layer**, the pipeline translates unsupervised cluster markers into clinically validated therapeutic candidates. 

[cite_start]By bridging single-cell transcriptomics with real-time semantic queries across **NCBI PubMed** and the **ChEMBL REST API** [cite: 334, 335][cite_start], the architecture calculates a composite multi-criteria priority score[cite: 318]. [cite_start]It incorporates robust mathematical shields against un-druggable target populism, algorithmic sifting of technical noise, and information-theoretic dynamic criteria weighting[cite: 287, 315].

---

## 🛠️ Workflow Architecture & Data Transformation Mechanics

[cite_start]The core engine implements a rigid, seven-stage sequential data transformation network. [cite_start]Every phase transition is guarded by explicit state validation checkpoints to enforce mathematical integrity over the cell-to-target mapping.

```text
[Raw Count Matrix] ──> Stage 1: Quality Control & Filtering
                             │ (Mitochondrial & Ribosomal Gating)
                             ▼
                      Stage 2: Library-Size Equalization & log1p
                             │ (Total count normalization to 1e4)
                             ▼
                      Stage 3: Highly Variable Gene (HVG) Selection
                             │ (Mean-Variance dispersion pruning)
                             ▼
                      Stage 4: Linear Regression & ARPACK PCA
                             │ (Covariate removal & Dimensionality Reduction)
                             ▼
                      Stage 5: Community Detection & Manifold Embedding
                             │ (kNN Graph -> Leiden Partitioning -> UMAP)
                             ▼
                      Stage 6: Non-Parametric DEG Testing
                             │ (Wilcoxon Rank-Sum Cluster Profiling)
                             ▼
                      Stage 7: Agentic AI Prioritization Matrix
                               (Multi-API Harvesting & Shannon Entropy Scoring)
