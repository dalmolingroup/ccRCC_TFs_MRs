# Transcriptional master regulators associated with metastatic ccRCC

Code for the manuscript **"Cross-Cohort Network Analysis Identifies Transcription Factor Master Regulators of Aggressiveness and Survival in Clear Cell Renal Cell Carcinoma"** (Farias *et al.*, submitted to *Computational and Structural Biotechnology Journal*, CSBJ-D-26-01039, revised version).

The workflow infers a gene regulatory network (GRN) for clear cell renal cell carcinoma (ccRCC) from TCGA-KIRC, identifies transcription factors (TFs) acting as master regulators (MRs) of a tumour–normal signature derived from paired ICGC-RECA samples, selects with Boruta/Random Forest a signature of 14 MRs associated with metastasis at diagnosis, groups them into two regulon clusters (proliferative and differentiation-associated) and characterises these clusters by regulon activity, expression, functional enrichment, external validation of regulon targets, clinicopathological associations, survival, VHL status and tumour composition.

An archived snapshot of this repository is available on Zenodo: 10.5281/zenodo.24316449 .

---

## Repository structure

```
data/                         GDC manifest and instructions to obtain the input data
R/01_data/                    data download and preparation
R/02_network_MRA/             regulatory network (RTN), differential expression, master regulator analysis
python/03_machine_learning/   signature selection and evaluation (Boruta, Random Forest, nested CV, baselines)
R/04_activity_clusters/       regulon activity, regulon clusters, ICGC activity, two-tailed GSEA
R/05_expression/              expression of the signature TFs (tumour vs normal; M1 vs M0)
R/06_enrichment_validation/   functional enrichment and validation of regulon targets (ChIP-Atlas, DoRothEA)
R/07_survival_clinical/       clinical associations, Cox models, VHL subgroups, purity/deconvolution, ICGC survival
R/08_tables_flowchart/        Tables 1 and 2, patient flowchart
```

## Running order

Scripts are numbered in the order in which they should be run. All R scripts assume that the project root is `~/Kidney/`, with input and intermediate data in `~/Kidney/0_data/`, and that the working directory is `~/Kidney/` (`setwd("~/Kidney")`). Paths and parameters are defined at the top of each script (`CONFIG` block).

| # | Script | Main outputs | Manuscript |
|---|---|---|---|
| 01 | `R/01_data/01_datasets_construction.R` | TCGA-KIRC counts and clinical data; ICGC-RECA counts and donor data (`0_data/*.RData`) | Methods 2.1 |
| 02 | `R/01_data/02_transcription_factor_list.R` | TF list (AnimalTFDB 4.0) | Methods 2.2 |
| 03 | `R/02_network_MRA/03_rtn_network.R` | RTN network (`rtni_TF.RData`) | Methods 2.2; Results 3.1 |
| 04a–c | `R/02_network_MRA/04a_de_deseq2.R`, `04b_de_edger.R`, `04c_de_limma.R` | paired tumour–normal DE in ICGC-RECA | Methods 2.2 |
| 04d | `R/02_network_MRA/04d_upset_de_FigS2.R` | intersection of DE genes | Fig. S2 |
| 05 | `R/02_network_MRA/05_master_regulator_analysis.R` | MRA (`MRA.RData`, 356 MRs); input matrices for the ML step (`expr_MR_KIRC_full.csv`, `metadata_KIRC_full.csv`, `expr_all_regulators_voom_KIRC.csv`) | Results 3.1; Table S2 |
| 06 | `python/03_machine_learning/ml_pipeline_revision.py` | split, Boruta signature, test-set evaluation, nested CV, stability, ICGC evaluation (`results_ml_revision/`). Run `python ml_pipeline_revision.py`, then `python ml_pipeline_revision.py baselines` | Methods 2.3; Results 3.2 |
| 07 | `python/03_machine_learning/clinical_baseline.ipynb` | clinical model, combined model, multivariable logistic regression | Results 3.2 |
| 08 | `python/03_machine_learning/plots_ml_revision.py` | Figure 3; `python plots_ml_revision.py compare` for Fig. S4 | Fig. 3, Fig. S4 |
| 09 | `python/03_machine_learning/make_supp_tables_ml.py` | split/fold and baseline tables | Tables S3, S4 |
| 10 | `R/04_activity_clusters/10_regulon_activity_graphs_Fig2_Fig4_FigS3.R` | network graphs (`tni.graph`), association map of the 14 regulons (clusters), regulon activity heatmap | Fig. 2, Fig. 4A–B, Fig. S3 |
| 11 | `R/04_activity_clusters/11_activity_vs_clinical_Fig4C.R` | cluster activity vs metastasis, grade, pT and stage | Fig. 4C; Table S16 |
| 12 | `R/04_activity_clusters/12_FigS5_ICGC_activity.R` | regulon activity in ICGC-RECA, adjusted Rand index | Fig. S5 |
| 13 | `R/04_activity_clusters/13_FigS6_GSEA2.R` | two-tailed GSEA of the 14 regulons | Fig. S6 |
| 14 | `R/05_expression/14_expression_Fig5AB_TabS6.R` | expression tumour vs normal (paired) and M1 vs M0 | Fig. 5A–B; Table S6 |
| 15 | `R/05_expression/15_Fig5C_concordance.R` | concordance of M1 vs M0 effects, TCGA vs ICGC | Fig. 5C |
| 16 | `R/06_enrichment_validation/16_enrichment_and_regulon_validation.R` | ORA per regulon and per cluster/sign (network background); CollecTRI/DoRothEA and ChIP-Atlas overlap | Methods 2.6–2.7 |
| 17 | `R/06_enrichment_validation/17_Fig6_TabS7-S12.R` | Figure 6 and enrichment/validation tables | Fig. 6; Tables S7–S12 |
| 18 | `R/07_survival_clinical/18_survival_cox_TCGA.R` | multivariable Cox models (OS, PFI, DSS; all and M0 patients) | Methods 2.8 |
| 19 | `R/07_survival_clinical/19_VHL_subgroups.R` | VHL-mutant vs wild-type analyses | Methods 2.8 |
| 20 | `R/07_survival_clinical/20_purity_deconvolution.R` | correlations with purity and immune fractions; adjusted models | Methods 2.9 |
| 21 | `R/07_survival_clinical/21_survival_ICGC.R` | exploratory OS analysis in ICGC-RECA | Methods 2.8 |
| 22 | `R/07_survival_clinical/22_Fig7_FigS7-S9_TabS13-S25.R` | Figure 7, Kaplan–Meier, VHL and purity figures, survival/clinical tables | Fig. 7; Figs. S7–S9; Tables S13–S25 |
| 23 | `R/08_tables_flowchart/23_Table1.R` | clinical characteristics | Table 1 |
| 24 | `R/08_tables_flowchart/24_Table2.R` | summary of the 14 signature TFs | Table 2 |
| 25 | `R/08_tables_flowchart/25_FigS1_flowchart.R` | patient flowchart | Fig. S1 |

Scripts 18–21 can be run in any order; script 22 must be run after them. Before running script 22, the per-patient data files must exist (`dados_por_paciente.csv`, `dados_sobrevivencia.csv`, `dados_VHL_por_paciente.csv`), which are written by scripts 11, 18 and 19.

**Figures not generated by code.** Figure 1 was drawn with BioRender. The network layouts of Figure 2 and Fig. S3 were rendered interactively in RedeR from the graph objects created by script 10 (`rtn_graph_Boruta_sign.RData`). The panels of Figure 5 were combined with `cowplot`.

## Input data

Input data are not redistributed in this repository.

- **TCGA-KIRC**: RNA-seq (STAR counts), clinical and somatic mutation data from the NCI Genomic Data Commons (data release 40.0), downloaded with TCGAbiolinks. The file list is given in `data/gdc_manifest_TCGA-KIRC.txt`.
- **TCGA Pan-Cancer Clinical Data Resource** (Liu *et al.*, *Cell* 2018; Table S1, sheet "TCGA-CDR") for survival endpoints: `0_data/TCGA-CDR-SupplementalTableS1.xlsx`.
- **ICGC RECA-EU** (release 28): `exp_seq.RECA-EU.tsv.gz` and `donor.RECA-EU.tsv.gz`.
- **Tumour purity and immune deconvolution**: ESTIMATE, quanTIseq and CIBERSORT (LM22, 1,000 permutations) estimates computed with IOBR (v2.2.4) from TCGA-KIRC STAR counts (UCSC Xena) converted to TPM, saved as `0_data/10_tme_combine.rds`.
- **External regulatory evidence**: ChIP-Atlas target genes (hg38, TSS ± 5 kb; downloaded by script 16) and DoRothEA (confidence levels A–C; R package `dorothea`).

## Software

- **R 4.5.1**: RTN 2.26.0, RTNsurvival, RTNduals, TCGAbiolinks 2.30.4, GDCRNATools, DESeq2, edgeR, limma, clusterProfiler 4.18.4, org.Hs.eg.db 3.22.0, GO.db 3.22.0, dorothea, decoupleR, OmnipathR, survival 3.8.12, survminer 0.5.2, IOBR 2.2.4, gtsummary 2.6.1, tidyverse, ggplot2, cowplot, patchwork, ggh4x, ggrepel, pheatmap, ComplexHeatmap, paletteer, readxl, openxlsx, kableExtra, flextable, finalfit, vroom.
- **Python 3.11.15**: see `requirements.txt`.

The exact package versions of the R session are listed in `sessionInfo.txt`.

## Reproducibility

- Random seeds are fixed in all stochastic steps (`SEED = 2026` in the machine-learning workflow; `set.seed(2026)` in R scripts with permutations or resampling).
- The train/test split, the selected features, the selection frequencies across resamplings and nested-CV folds, and all predictions are written to `results_ml_revision/`.
- The regulon membership of the 14 signature TFs, with the mode of regulation of each target, is provided as Supplementary Table S12.

## License

MIT (see `LICENSE`).

## Citation

If you use this code, please cite the article [reference/DOI] and the Zenodo archive 10.5281/zenodo.24316449.
