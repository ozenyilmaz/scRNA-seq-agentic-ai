# Multi-Criteria Therapeutic Target Prioritization from single-cell RNA-Seq
Using Agentic AI Knowledge-Harmonization and Shannon Entropy Weighting

[![Open In Colab](https://img.shields.io/badge/Colab-Run%20Pipeline-orange?style=for-the-badge&logo=google-colab)](https://colab.research.google.com/github/KULLANICI_ADIN/scRNA-seq-target-prioritization/blob/main/)

## Executive Summary
his repository contains a modular, production-grade single-cell RNA-sequencing (scRNA-seq) framework designed for drug target prioritization Engineered as an **Agentic AI Knowledge-Harmonization Layer**, the pipeline translates unsupervised cluster markers into clinically validated therapeutic candidates. 

By bridging single-cell transcriptomics with real-time semantic queries across **NCBI PubMed** and the **ChEMBL REST API**  the architecture calculates a composite multi-criteria priority scoreIt incorporates robust mathematical shields against un-druggable target populism, algorithmic sifting of technical noise, and information-theoretic dynamic criteria weighting.

---

## 🛠️ Workflow Architecture & Data Transformation Mechanics

The core engine implements a rigid, seven-stage sequential data transformation network. Every phase transition is guarded by explicit state validation checkpoints to enforce mathematical integrity over the cell-to-target mapping.

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


┌─────────────────────────────────────────────────────────┐
       │             KNOWLEDGE-HARMONIZATION AGENT               │
       └────────────────────────────┬────────────────────────────┘
                                    │
    [SENSE]  ──────> Reads Dynamic DEG Cluster Matrices
                                    │
    [THINK]  ──────> Invokes Shannon Entropy to calculate Diversity
             ──────> Applies Logistic Sigmoid Transforms on Clinical Data
                                    │
    [ACT]    ──────> Hits PubMed & ChEMBL REST Endpoints
             ──────> Auto-corrects Type Mismatches ('4.0' -> 4)
             ──────> Generates Verified Drug Target Priorities

The four-phase roadmap (§6) provides a technically rigorous path toward a high-sensitivity cardiac drug-discovery instrument: replacing the PBMC reference with a cardiac-native single-cell atlas (§6.1) will simultaneously resolve the biological domain gap and unlock the full discriminative potential of the ChEMBL and PubMed enrichment layers, while per-cluster normalization (§6.2), curve selection (§6.3), and multi-database extension (§6.4) will furnish the rank-permutation sensitivity required for a publication-grade cardiovascular genomics analysis.
