#### Supplementary Table S26: rede e MRA inferidas só no treino (análise de sensibilidade) ####
## Revisão CSBJ-D-26-01039 - Revisor 2, Comentário 2
##
## Junta as saídas da análise de sensibilidade em uma tabela (xlsx, 4 abas):
##   Summary           sobreposição dos MRs (rede completa x rede só-treino)
##   Signature_MRA     FDR do MRA dos 14 TFs da assinatura nas duas redes
##   Selected_features reguladores confirmados pelo Boruta com os candidatos da rede só-treino
##   Holdout           métricas no MESMO teste intocado: análise original x rede só-treino
##
## Ordem: export_split_ids.py -> 03b_rtn_network_train_only.R -> ml_pipeline_trainnet.py -> este script
## Saída: ~/Kidney/results_ml_trainnet/TableS26_trainnet_sensitivity.xlsx

suppressPackageStartupMessages({
  library(RTN); library(jsonlite); library(openxlsx); library(AnnotationDbi); library(org.Hs.eg.db)
})

K        <- path.expand("~/Kidney")
OUT      <- file.path(K, "results_ml_trainnet")
BAL      <- "class_weight"
strip_v  <- function(x) sub("\\..*", "", x)
sym      <- function(ens) unname(mapIds(org.Hs.eg.db, strip_v(ens), "SYMBOL", "ENSEMBL", multiVals = "first"))

## ---------------------------------------------------------------- MRA: rede completa x só-treino
mra_full  <- read.csv(file.path(K, "0_data/MRA_TFs.csv"), row.names = 1, check.names = FALSE)
mra_train <- read.csv(file.path(K, "0_data/trainnet/MRA_TFs_train.csv"), row.names = 1, check.names = FALSE)
load(file.path(K, "0_data/trainnet/MRA_train.RData"))                    # rtna_train
mra_train_all <- tna.get(rtna_train, what = "mra", ntop = -1)            # todos os regulons testados

id_full <- strip_v(rownames(mra_full)); id_train <- strip_v(rownames(mra_train))
inter   <- intersect(id_full, id_train)
tab_sum <- data.frame(
  Measure = c("Tumours used for network inference", "Master regulators (adjusted p < 0.01)",
              "Shared master regulators", "Shared, % of original", "Shared, % of training-only",
              "Jaccard index", "Signature regulators retained as candidates"),
  `Full cohort (original)` = c("501", length(id_full), length(inter), "", "", "", "14/14"),
  `Training set only` = c("350", length(id_train), length(inter),
                          sprintf("%.1f", 100 * length(inter) / length(id_full)),
                          sprintf("%.1f", 100 * length(inter) / length(id_train)),
                          sprintf("%.2f", length(inter) / length(union(id_full, id_train))), ""),
  check.names = FALSE)

## ---------------------------------------------------------------- 14 TFs da assinatura
sig <- strip_v(read.csv(file.path(K, "results_ml_revision", paste0("importancia_features_modelo_final_", BAL, ".csv")))[[1]])
res_tr <- fromJSON(file.path(OUT, paste0("resultados_", BAL, ".json")))
res_or <- fromJSON(file.path(K, "results_ml_revision", paste0("resultados_", BAL, ".json")))
conf_tr <- strip_v(res_tr$holdout$boruta_confirmed)

adjp <- function(tab, ids) {
  r <- match(ids, strip_v(rownames(tab)))
  ifelse(is.na(r), NA, tab$Adjusted.Pvalue[r])
}
tab_sig <- data.frame(TF = sym(sig), Ensembl = sig,
                      `MRA adj. p, full cohort` = signif(adjp(mra_full, sig), 2),
                      `MRA adj. p, training only` = signif(adjp(mra_train_all, sig), 2),
                      `Candidate in training-only network` = ifelse(sig %in% id_train, "Yes", "No"),
                      `Selected by Boruta (training-only candidates)` = ifelse(sig %in% conf_tr, "Yes", "No"),
                      check.names = FALSE)
tab_sig <- tab_sig[order(tab_sig$`Candidate in training-only network`, tab_sig$`MRA adj. p, training only`), ]
tab_sum$`Training set only`[tab_sum$Measure == "Signature regulators retained as candidates"] <-
  sprintf("%d/14", sum(sig %in% id_train))

## ---------------------------------------------------------------- features confirmadas (rede só-treino)
tab_feat <- data.frame(Ensembl = conf_tr, TF = sym(conf_tr),
                       `In original signature` = ifelse(conf_tr %in% sig, "Yes", "No"),
                       `MRA adj. p, training only` = signif(adjp(mra_train_all, conf_tr), 2),
                       check.names = FALSE)
tab_feat <- tab_feat[order(tab_feat$`In original signature`, decreasing = TRUE), ]

## ---------------------------------------------------------------- holdout (mesmo teste)
met <- function(m, label, analysis, n_feat) {
  data.frame(Analysis = analysis, Model = label, Features = n_feat,
             AUC = sprintf("%.2f (%.2f–%.2f)", m$AUC, m$AUC_ci95[1], m$AUC_ci95[2]),
             `PR-AUC` = sprintf("%.2f", m$PR_AUC), Sensitivity = sprintf("%.2f", m$sensitivity),
             Specificity = sprintf("%.2f", m$specificity), check.names = FALSE)
}
tab_hold <- rbind(
  met(res_or$holdout$modelo_selecionado, "Signature (Boruta)", "Original (full-cohort network)",
      length(res_or$holdout$boruta_confirmed)),
  met(res_or$holdout$modelo_356_TFs, "All master regulators", "Original (full-cohort network)", length(id_full)),
  met(res_tr$holdout$modelo_selecionado, "Signature (Boruta)", "Training-only network",
      length(res_tr$holdout$boruta_confirmed)),
  met(res_tr$holdout$modelo_356_TFs, "All master regulators", "Training-only network", length(id_train)))
tab_hold$`Test set` <- sprintf("n = %d (%d M1)", res_or$holdout$n_teste, res_or$holdout$modelo_selecionado$n_M1)

## ---------------------------------------------------------------- salvar
wb <- createWorkbook()
for (sh in list(list("Summary", tab_sum), list("Signature_MRA", tab_sig),
                list("Selected_features", tab_feat), list("Holdout", tab_hold))) {
  addWorksheet(wb, sh[[1]]); writeData(wb, sh[[1]], sh[[2]])
  setColWidths(wb, sh[[1]], cols = seq_len(ncol(sh[[2]])), widths = "auto")
}
saveWorkbook(wb, file.path(OUT, "TableS26_trainnet_sensitivity.xlsx"), overwrite = TRUE)
print(tab_sum, row.names = FALSE); print(tab_hold, row.names = FALSE)
cat("\nSalvo:", file.path(OUT, "TableS26_trainnet_sensitivity.xlsx"), "\n")
