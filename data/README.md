# Input data

Place the input files in `~/Kidney/0_data/` (see the main README, section "Input data").

- `gdc_manifest_TCGA-KIRC.txt`: GDC manifest of the TCGA-KIRC RNA-seq files (STAR counts) used in the study.
  It can be used with the GDC Data Transfer Tool (`gdc-client download -m gdc_manifest_TCGA-KIRC.txt`);
  script `R/01_data/01_datasets_construction.R` downloads the same data with TCGAbiolinks.
