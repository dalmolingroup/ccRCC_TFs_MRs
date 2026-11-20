#### Regulons Activity Analysis ####
## Date: 04/08/2024
## Author: Epitácio Farias 

#### Loading Packages ####
if(!require("RTN")){BiocManager::install("RTN")}
if(!require("RTNduals")){BiocManager::install("RTNduals")}
if(!require("RTNsurvival")){BiocManager::install("RTNsurvival")}
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("finalfit")){install.packages("finalfit")}
if(!require("survminer")){install.packages("survminer")}
if(!require("clusterProfiler")){BiocManager::install("clusterProfiler")}
if(!require("org.Hs.eg.db")){BiocManager::install("org.Hs.eg.db")}

#### Loading the data ####
# 14 TFs da assinatura final (saída do 06_ml_pipeline.py; 1ª coluna = Ensembl ID)
selected_features <- vroom::vroom("~/Kidney/results_ml_revision/importancia_features_modelo_final_class_weight.csv")
colnames(selected_features)[1] <- "ENSEMBL_ID"
load("~/Kidney/0_data/rtn_TF_data.RData")
load("~/Kidney/0_data/MRA.RData")

selected_features <- bitr(selected_features$ENSEMBL_ID, fromType= "ENSEMBL", toType = c("SYMBOL","ENTREZID","ENSEMBL"),
                          OrgDb = org.Hs.eg.db, drop = T)
selected_features <- selected_features[!duplicated(selected_features$ENSEMBL), ]

#### Retrieving the Regulons ####
targets_regs <- regulons_TF[names(regulons_TF) %in% selected_features$ENSEMBL]

extract_names <- function(my_list, vector_name) {
  if (vector_name %in% names(my_list)) {
    names(my_list[[vector_name]])
  } else {
    NULL
  }
}

for (i in 1:length(selected_features$ENSEMBL)){
  targets_regs[[i]] <- paste(extract_names(targets_regs, names(targets_regs[i])))
  names(targets_regs)[i] <- names(targets_regs)[i]
}

rm(regulons_TF, extract_names, i)

#### Regulon Interactions ####
g <- tni.graph(rtni_TF, regulatoryElements = names(targets_regs))

g_all <- tni.graph(rtni_TF, regulatoryElements = regs$GENE_ID,tnet = "dpi", 
                   minRegulonSize = 15, gtype = "amapDend")

g_MR <- tni.graph(rtni_TF, tnet = "dpi", minRegulonSize = 15, gtype = "amapDend",
                  regulatoryElements = mra_kirc_prim_norm$Regulon)

g_r <- tni.graph(rtni_TF, tnet = "dpi", minRegulonSize = 15,
                 gtype = "amap",
                 regulatoryElements = names(targets_regs), amapFilter = "phyper",
                 amapCutoff = 1)

g_tal <- tni.graph(rtni_TF, tnet = "dpi", minRegulonSize = 15, gtype = "amapDend")

# g_MR: rede dos 356 MRs (Figura 2, renderizada no RedeR); g_all/g_tal: rede completa (Fig S3, RedeR);
# g_r: mapa de associação (Jaccard) dos 14 regulons -> clusters (Figura 4A)
save(g, g_r, g_tal, g_MR, g_all, file = "~/Kidney/0_data/rtn_graph_Boruta_sign.RData")
rm(regs)

#### Co-regulation analysis and inference of dual regulons ####
rmbr_TF <- tni2mbrPreprocess(rtni_TF)

rmbr_TF <- mbrAssociation(rmbr_TF)

mbrGet(rmbr_TF, what="summary")

## Filtering by the TF from the best model constructed
## Overlap Situation
overlap <- mbrGet(rmbr_TF, what="dualsOverlap")

overlap2 <- overlap[overlap$Regulon1 %in% selected_features$SYMBOL | overlap$Regulon2 %in% selected_features$SYMBOL, ]

## Correlation
correlation <- mbrGet(rmbr_TF, what="dualsCorrelation")
correlation2 <- correlation[correlation$Regulon1 %in% selected_features$SYMBOL | correlation$Regulon2 %in% selected_features$SYMBOL, ]

## Duals Relation
all(order(rownames(correlation2)) %in% order(rownames(overlap2)))
all(order(rownames(overlap2)) %in% order(rownames(correlation2)))

duals <- mbrGet(rmbr_TF, what="dualRegulons")

duals2 <- duals[duals %in% rownames(overlap2)]

save(rmbr_TF,correlation,correlation2, duals, duals2, file = "~/Kidney/0_data/coregulation_data_Boruta_sign.RData")


#### Computing Regulons Activity (Figura 4B) ####
library("tidyverse")
if(!require("GDCRNATools")){BiocManager::install("GDCRNATools")}
library("RTN")
library("pheatmap")
library("patchwork")
library("paletteer")
load("~/Kidney/0_data/DEG_DESeq2.RData")
load("~/Kidney/0_data/rtni_TF.RData")
load("~/Kidney/0_data/ICGC_DEA_Data.RData")

selected_features <- vroom::vroom("~/Kidney/results_ml_revision/importancia_features_modelo_final_class_weight.csv")

colnames(selected_features)[1] <- "ENSEMBL_ID"

###### HEATMAP TCGA-KIRC DATASET #####
# Compute regulon activity for individual samples
rtni_TF <- tni.gsea2(rtni_TF, regulatoryElements = selected_features$ENSEMBL_ID)
tcga_regact <- tni.get(rtni_TF, what = "regulonActivity")

c2 <- c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4")
annotation_row = data.frame(
  MRs = rownames(t(tcga_regact$dif))
)

annotation_row$cluster <- ifelse(annotation_row$MRs %in% c2, "C2","C1")
annotation_row <- annotation_row %>%  column_to_rownames("MRs")
my_colors <- paletteer_d("MexBrewer::Alacena")

# Get sample attributes from the 'rtni1st' dataset
kirc_annot <- tni.get(rtni_TF, "colAnnotation")
kirc_annot$ajcc_pathologic_m <- ifelse(kirc_annot$ajcc_pathologic_m == "1", "M1","M0")

# Get ER+/- and PAM50 attributes for pheatmap
attribs <- c("ajcc_pathologic_m","vital_status")
kirc_annot <- kirc_annot[,attribs]

reorder_callback <- function(hc, mat) {
  
  if (!any(hc$labels %in% rownames(annotation_row))) {
    return(hc)
  }
  
  dend <- as.dendrogram(hc)
  clusters <- cutree(hc, k = 2)
  
  c1_members <- clusters[names(clusters) %in% rownames(annotation_row)[annotation_row$cluster == "C1"]]
  c2_members <- clusters[names(clusters) %in% rownames(annotation_row)[annotation_row$cluster == "C2"]]
  
  c1_cluster_id <- as.integer(names(sort(table(c1_members), decreasing = TRUE))[1])
  c2_cluster_id <- as.integer(names(sort(table(c2_members), decreasing = TRUE))[1])
  
  # Rotaciona quando C1 está embaixo (id menor que C2)
  if (c1_cluster_id < c2_cluster_id) {
    dend <- rotate(dend, order = rev(labels(dend)))
  }
  
  as.hclust(dend)
}


# Figura 4B - 1 tumor primário (01) por paciente (501), como em todas as outras análises
# (o objeto da rede inclui alíquotas extras do mesmo paciente; revisão R2, minor 1).
# Agrupamento dos TFs como na Seção 2.4 (1 - Spearman, ligação média, k = 2).
dif <- tcga_regact$dif
dif <- dif[substr(rownames(dif), 14, 15) == "01", ]
dif <- dif[!duplicated(substr(rownames(dif), 1, 12)), ]
stopifnot(nrow(dif) == 501)
d_rows <- as.dist(1 - cor(dif, method = "spearman"))

tiff(filename = "~/Kidney/Figures/Boruta_Sign_Figures/Heatmap_Act_TCGA.tiff",
     res = 300, height = 1000, width = 3000)
a <- pheatmap(t(dif),
         main = sprintf("TCGA-KIRC (n = %d tumours)", nrow(dif)),
         show_colnames = FALSE, annotation_legend = TRUE,
         cluster_rows = hclust(d_rows, method = "average"), cutree_rows = 2,
         clustering_distance_cols = "manhattan", clustering_method = "ward.D2",   # colunas (tumores)
         annotation_row = annotation_row,
         annotation_colors = list(cluster = c("C1" = "#542788", "C2" = "#1a9850")),  # C1 roxo, C2 verde
         color = my_colors, border_color = NA, annotation_names_row = FALSE)
dev.off()
