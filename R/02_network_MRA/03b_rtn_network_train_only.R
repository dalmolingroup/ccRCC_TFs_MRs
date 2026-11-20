#### Análise de sensibilidade: rede e MRA inferidas SÓ no treino ####
## Revisão CSBJ-D-26-01039 - Revisor 2, Comentário 2
## "perform network inference and master regulator prioritization within the corresponding
##  training partitions, at least as a sensitivity analysis"
##
## Repete 03_rtn_network.R + 05_master_regulator_analysis.R com os mesmos parâmetros, mas usando
## apenas os tumores do TREINO do split principal (70%; IDs de export_split_ids.py).
## A normalização (voom) também é feita só no treino. A assinatura tumor x normal do MRA
## (pares do ICGC-RECA) não usa dados do TCGA e fica inalterada.
##
## Saídas (~/Kidney/0_data/trainnet/) - consolidadas na Supplementary Table S26 pelo 03c_TableS26_trainnet.R:
##   rtni_TF_train.RData, MRA_train.RData, MRA_TFs_train.csv
##   overlap_MR_train_vs_full.txt   sobreposição com os 356 candidatos e com os 14 TFs da assinatura
##   expr_MR_KIRC_trainnet.csv      entrada do Python: mesmas amostras do expr_MR_KIRC_full.csv,
##                                  colunas = MRs da rede só-treino
##
## Rodar no HPC em segundo plano, ex.:  nohup Rscript 03b_rtn_network_train_only.R > trainnet.log 2>&1 &
## Cada etapa salva um checkpoint; se cair, o script retoma da última etapa concluída.

# sem tidyverse/GDCRNATools: evitam o conflito de versões do ggplot2 no container (3.5.2 x >= 4.0.0)
suppressPackageStartupMessages({ library(RTN); library(tibble); library(parallel); library(edgeR); library(limma) })

# idêntico ao GDCRNATools::gdcVoomNormalization(counts, filter = FALSE): TMM + voom, sem design
voom_norm <- function(counts) {
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts))
  limma::voom(dge, design = NULL, plot = FALSE)$E
}

## ---------------------------------------------------------------- CONFIG
OUT      <- path.expand("~/Kidney/0_data/trainnet")
TRAIN_IDS <- path.expand("~/Kidney/0_data/split_train_ids.csv")
N_PERM   <- 1000     # igual ao original. 100 reduz bastante o tempo (o RTN usa nulo agrupado
                     # entre todos os pares TF-gene, pooledNullDistribution = TRUE); se usar 100,
                     # declare na resposta ao revisor
N_CORES  <- 32       # aumente se o nó tiver mais núcleos
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ck <- function(f) file.path(OUT, f)

## ---------------------------------------------------------------- 1. dados (como no 03_rtn_network.R)
load("~/Kidney/0_data/KIRC_data_updated.RData")            # clin_data3, gene_library, count_data_kirc2
load("~/Kidney/0_data/tf_list.RData")
tf_ensembl <- tf_list$Ensembl[!is.na(tf_list$Ensembl)]

ColAnnot <- clin_data3 |> rownames_to_column("patient_ID")
RowAnnot <- gene_library[gene_library$gene_id %in% rownames(count_data_kirc2), ]
names(RowAnnot)[3] <- "SYMBOL"

# amostras do treino (IDs de 15 caracteres no Python)
train_ids <- read.csv(TRAIN_IDS)$sample
cols15    <- substr(colnames(count_data_kirc2), 1, 15)
keep      <- colnames(count_data_kirc2) %in% ColAnnot$patient_ID & cols15 %in% train_ids
keep      <- keep & !duplicated(cols15)
cat(sprintf("Treino: %d IDs no arquivo | %d amostras encontradas na matriz de contagens\n",
            length(train_ids), sum(keep)))
stopifnot(sum(keep) == length(train_ids))

# voom SÓ no treino
exp_train <- voom_norm(count_data_kirc2[, keep])
ColAnnot_train <- ColAnnot[match(colnames(exp_train), ColAnnot$patient_ID), ]

## ---------------------------------------------------------------- 2. rede (permutação, bootstrap, DPI)
if (!file.exists(ck("rtni_TF_train.RData"))) {
  if (file.exists(ck("rtni_train_boot.RData"))) {
    load(ck("rtni_train_boot.RData"))
  } else {
    if (file.exists(ck("rtni_train_perm.RData"))) {
      load(ck("rtni_train_perm.RData"))
    } else {
      rtni_train <- tni.constructor(expData = exp_train, regulatoryElements = tf_ensembl,
                                    colAnnotation = ColAnnot_train, rowAnnotation = RowAnnot)
      options(cluster = snow::makeCluster(spec = N_CORES, "SOCK"))
      message(Sys.time(), " permutação (", N_PERM, ")...")
      rtni_train <- tni.permutation(rtni_train, nPermutations = N_PERM)
      save(rtni_train, file = ck("rtni_train_perm.RData"))
    }
    if (is.null(getOption("cluster"))) options(cluster = snow::makeCluster(spec = N_CORES, "SOCK"))
    message(Sys.time(), " bootstrap...")
    rtni_train <- tni.bootstrap(rtni_train)
    save(rtni_train, file = ck("rtni_train_boot.RData"))
  }
  message(Sys.time(), " DPI...")
  rtni_train <- tni.dpi.filter(rtni_train)
  try(stopCluster(getOption("cluster")), silent = TRUE)
  save(rtni_train, file = ck("rtni_TF_train.RData"))
} else {
  load(ck("rtni_TF_train.RData"))
}
print(tni.regulon.summary(rtni_train))

## ---------------------------------------------------------------- 3. MRA (como no 05_master_regulator_analysis.R)
load("~/Kidney/0_data/DEG_DESeq2.RData")                 # DEG_genes_final (pares do ICGC-RECA)
fold_paired   <- setNames(DEG_genes_final$log2FoldChange, DEG_genes_final$gene_name)
symbol_paired <- DEG_genes_final$gene_name

rtna_train <- tni2tna.preprocess(object = rtni_train, phenotype = fold_paired, hits = symbol_paired)
rtna_train <- tna.mra(rtna_train, pValueCutoff = 0.01, pAdjustMethod = "BH",
                      minRegulonSize = 15, tnet = "dpi")
mra_train <- tna.get(rtna_train, what = "mra")
save(rtna_train, mra_train, file = ck("MRA_train.RData"))
write.csv(mra_train, ck("MRA_TFs_train.csv"))

## ---------------------------------------------------------------- 4. sobreposição com a análise original
strip_v <- function(x) sub("\\..*", "", x)
mra_full <- read.csv("~/Kidney/0_data/MRA_TFs.csv", row.names = 1)
sig      <- read.csv("~/Kidney/results_ml_revision/importancia_features_modelo_final_class_weight.csv")[[1]]

id_train <- strip_v(rownames(mra_train)); id_full <- strip_v(rownames(mra_full)); id_sig <- strip_v(sig)
inter    <- intersect(id_train, id_full)
sym      <- setNames(mra_full$Regulon, strip_v(rownames(mra_full)))

txt <- c(
  sprintf("MRs - rede completa (501): %d | rede só-treino (%d): %d", length(id_full), ncol(exp_train), length(id_train)),
  sprintf("Em comum: %d (%.1f%% dos originais; %.1f%% dos novos); Jaccard = %.2f",
          length(inter), 100 * length(inter) / length(id_full), 100 * length(inter) / length(id_train),
          length(inter) / length(union(id_train, id_full))),
  sprintf("Assinatura (14 TFs) entre os MRs da rede só-treino: %d/14", sum(id_sig %in% id_train)),
  paste("Ausentes:", paste(sym[id_sig[!id_sig %in% id_train]], collapse = ", ")))
writeLines(txt, ck("overlap_MR_train_vs_full.txt")); cat(txt, sep = "\n")

## ---------------------------------------------------------------- 5. entrada do Python
# mesma matriz e mesmas amostras do expr_MR_KIRC_full.csv; só muda o conjunto de candidatos
load("~/Kidney/0_data/rtni_constructor_TF.RData")
full_rows <- rownames(read.csv("~/Kidney/0_data/expr_MR_KIRC_full.csv", row.names = 1, check.names = FALSE))
X <- t(rtni_constructor_TF@gexp) |> as.data.frame()
X <- X[, strip_v(colnames(X)) %in% id_train, drop = FALSE]
rownames(X) <- substr(rownames(X), 1, 15)
X <- X[!duplicated(rownames(X)), , drop = FALSE][full_rows, , drop = FALSE]
stopifnot(!anyNA(X))
write.csv(X, ck("expr_MR_KIRC_trainnet.csv"))
cat(sprintf("\nexpr_MR_KIRC_trainnet.csv: %d amostras x %d MRs\n", nrow(X), ncol(X)))
