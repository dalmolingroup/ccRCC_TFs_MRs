#### TF Signature Construction ####
## Date: 09/14/2023
## Author: Epitácio Farias

#### Loading packages and data ####
library("RTN")

load("~/Kidney/0_data/DEG_DESeq2.RData")
load("~/Kidney/0_data/DEG_intersec.RData")
load("~/Kidney/0_data/tf_list.RData")
load("~/Kidney/0_data/rtn_TF_data.RData")

#### Correcting some issues ####
fold_paired <- DEG_genes_final$log2FoldChange
names(fold_paired) <- DEG_genes_final$gene_name

symbol_paired <- DEG_genes_final$gene_name

#### Creating RTN Analysis Object #### 
rtna_kirc_prim_norm <- tni2tna.preprocess(object = rtni_TF,
                                          phenotype = fold_paired,
                                          hits = symbol_paired)
#### Master Regulator Analysis ####
rtna_kirc_prim_norm <- tna.mra(rtna_kirc_prim_norm, pValueCutoff = 0.01, pAdjustMethod = "BH",
                              minRegulonSize = 15, tnet = "dpi")

rtna_kirc_prim_norm <- tna.gsea1(rtna_kirc_prim_norm, nPermutations=1000)
gsea1 <- tna.get(rtna_kirc_prim_norm, what="gsea1", ntop = -1)

rtna_kirc_prim_norm <- tna.gsea2(rtna_kirc_prim_norm, nPermutations=1000)
gsea2 <- tna.get(rtna_kirc_prim_norm, what="gsea2", ntop = -1)
gsea2_diff <- gsea2$differential

mra_kirc_prim_norm <- tna.get(rtna_kirc_prim_norm, what = "mra")

gsea1_MR <- gsea1[gsea1$Regulon %in% mra_kirc_prim_norm$Regulon,]
gsea2_diff_MR <- gsea2_diff[gsea2_diff$Regulon %in% mra_kirc_prim_norm$Regulon,]

mra_kirc_prim_norm[!(mra_kirc_prim_norm$Regulon %in% gsea1_MR$Regulon),]
mra_kirc_prim_norm[!(mra_kirc_prim_norm$Regulon %in% gsea2_diff_MR$Regulon),]

#### Saving Files ####
save(rtna_kirc_prim_norm, mra_kirc_prim_norm, file = "~/Kidney/0_data/MRA.RData")
write.csv(mra_kirc_prim_norm, file = "~/Kidney/0_data/MRA_TFs.csv")

#### Input matrix and metadata for the machine-learning step (06_ml_pipeline.py) ####
library("tidyverse")
mra_tf <- vroom::vroom("~/Kidney/0_data/MRA_TFs.csv")
metadata = vroom::vroom("~/Kidney/0_data/metadata_KIRC.csv")

load("~/Kidney/0_data/rtni_constructor_TF.RData")
load(file = "~/Kidney/0_data/KIRC_muta_dataset.RData")
kirc_clinic2 <- read_csv("~/Kidney/0_data/kirc.clinic2.csv")

kirc_clinic2 <- kirc_clinic2[,c(1,6)]
metadata$...1 <- stringr::str_sub(metadata$...1, 1, 15)
clin <- merge(metadata, kirc_clinic2, by.x = "...1", by.y = "sampleID", all.x = T)
clin <- clin[!duplicated(clin$...1),]

count_data <- t(rtni_constructor_TF@gexp) |> as.data.frame()
count_data <- count_data[, colnames(count_data) %in% mra_tf$...1]

meta <- clin |> filter(ajcc_pathologic_m %in% c("0", "1")) |> column_to_rownames("...1")
meta$ajcc_pathologic_m <- ifelse(meta$ajcc_pathologic_m == 0, "M0", "M1")

vhl_muta <- muta[muta$Hugo_Symbol == "VHL",c(1,10,16)]

vhl_muta$Tumor_Sample_Barcode <- vhl_muta$Tumor_Sample_Barcode %>% stringr::str_sub(1,16)
vhl_muta$VHL_Status <- 1

library(dplyr)

# Uma linha por amostra (junta os tipos se houver mais de uma mutação)
vhl_resumo <- vhl_muta %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(
    Variant_Type = paste(unique(Variant_Type), collapse = ";"),
    VHL_Status   = 1,
    .groups = "drop"
  )

# Junta ao meta; quem não está no vhl_muta recebe 0 / "None"
meta <- meta %>%
  left_join(vhl_resumo, by = c("sample" = "Tumor_Sample_Barcode")) %>%
  mutate(
    VHL_Status   = ifelse(is.na(VHL_Status), 0, VHL_Status),
    Variant_Type = ifelse(is.na(Variant_Type), "None", Variant_Type)
  )

table(meta$VHL_Status)

rownames(meta) <- stringr::str_sub(meta$sample,1,15)

count_data <- count_data %>% rownames_to_column("Sample")
count_data$Sample <- stringr::str_sub(count_data$Sample, 1,15)
count_data <- count_data[!duplicated(count_data$Sample),]
rownames(count_data) <- NULL
count_data <- count_data %>% column_to_rownames("Sample")
count_data <- count_data[intersect(rownames(meta), rownames(count_data)), ]

write.csv(meta, "~/Kidney/0_data/metadata_KIRC_full.csv")
write.csv(count_data, "0_data/expr_MR_KIRC_full.csv")

#### Expression of all regulators (baseline without network prioritisation) ####
# Usado em `python ml_pipeline_revision.py baselines` (Boruta + RF sobre todos os reguladores).
# CONFIRA: deve reproduzir o arquivo usado na revisão (mesmas amostras e normalização do expr_MR_KIRC_full.csv).
tfs_all <- unname(tni.get(rtni_constructor_TF, "regulatoryElements"))
expr_all <- t(rtni_constructor_TF@gexp[rownames(rtni_constructor_TF@gexp) %in% tfs_all, ]) |> as.data.frame()
rownames(expr_all) <- stringr::str_sub(rownames(expr_all), 1, 15)
expr_all <- expr_all[!duplicated(rownames(expr_all)), ]
expr_all <- expr_all[rownames(count_data), ]
write.csv(expr_all, "0_data/expr_all_regulators_voom_KIRC.csv")
