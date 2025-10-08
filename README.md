# 🧬 Biomarker Discovery Pipeline

An end-to-end, reproducible computational framework for transcriptomic biomarker discovery and validation.  
This repository demonstrates how raw RNA-seq or microarray data can be processed, analyzed, and modeled using a combination of statistical and machine learning methods.

---

## 🧠 Experience Gained & Responsibilities

- Analyzed raw RNA-seq data from **FASTQ** to quantified **counts** using splice-aware aligners (*STAR*, *HISAT2*) and quasi-mappers (*Salmon*, *Kallisto*) for robust gene expression profiling.
- Gained foundational understanding of transcriptomic technologies (**RNA-seq**, **microarray**) and primary data repositories including **GEO**, **SRA**, **NCBI**, and **ArrayExpress**.  
- Explored expression data formats (**FASTQ**) and applied preprocessing tools (*FastQC*, *Trimmomatic*, *HISAT2*, *STAR*, *featureCounts*) for RNA-seq alignment and quantification.  
- Retrieved and processed GEO datasets using **GEO2R**, followed by manual metadata curation for comparative transcriptomic analysis.  
- Performed **differential gene expression analysis** using *edgeR*, *limma*, and *DESeq2*, including normalization, outlier detection, and batch effect correction.  
- Applied downstream visualization and statistical methods: **heatmaps**, **boxplots**, **PCA**, **UMAP**, **t-SNE**, **hierarchical clustering**, **volcano plots**, and **phylogenetic analyses**.  
- Conducted **functional enrichment analysis** using *DAVID*, *STRING*, *Enrichr*, *Reactome*, and *GSEA* to identify enriched pathways and biological processes associated with **Mpox virus infection**.  
- Predicted **candidate therapeutic compounds** using *DSigDB* based on DEG signatures.  
- Validated **key DEGs** as potential biomarkers via **ROC curve analysis (AUC-based)**.  
- Developed and benchmarked **machine learning models** (*Random Forest*, *Extra Trees*, *SVM*) to identify predictive biomarkers, including model tuning and feature importance ranking.  
- Utilized a range of bioinformatics tools and libraries:  
  **R** (*DESeq2*, *limma*, *edgeR*, *GOplot*, *ggplot2*, *pROC*),  
  **Cytoscape** (*ClueGO*, *CluePedia*, *cytoHubba*, *MCODE*),  
  **Python** (*PyCaret*, *scikit-learn*, *pandas*), and **TBtools**.  
- **First-author publication:**  
  *Debnath, J.P., et al.* (2025). *Identification of potential biomarkers for 2022 Mpox virus infection: a transcriptomic network analysis and machine learning approach.* **Scientific Reports, 15**, 2922.

---

## 🧩 Project Workflow

```
1️⃣ Preprocessing & Quality Control    →  FastQC, Trimmomatic  
2️⃣ Alignment & Quantification         →  STAR / Salmon / featureCounts  
3️⃣ Differential Expression Analysis   →  DESeq2, edgeR, limma  
4️⃣ Functional Enrichment Analysis     →  GSEA, Enrichr, GOstats  
5️⃣ Network Modeling                   →  Identify hub genes, functional modules  
6️⃣ Machine Learning Benchmarking      →  PyCaret (AUPRC, AUROC, F1, MCC)  
```

---

## 📥 Data Access

Raw RNA-seq and microarray datasets used in this pipeline can be downloaded from public repositories:

- **SRA (Sequence Read Archive)**
- **GEO (Gene Expression Omnibus)**

---

## 📜 License

This project is distributed under the MIT License.  
See [LICENSE](LICENSE) for details.
