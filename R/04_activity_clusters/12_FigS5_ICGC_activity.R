###############################################################################
# Supplementary Figure S5 e Supplementary Table S27
# Coerência de atividade dos clusters regulatórios em TCGA-KIRC e ICGC-RECA
# Revisão CSBJ-D-26-01039 - Revisor 2, Comentários 1 e 3
#
# Network clusters (C1 = proliferative; C2 = differentiation-associated):
#   definidos no TCGA pela similaridade de Jaccard dos alvos (tni.graph, amap; script 10).
# Activity groups: agrupamento hierárquico dos 14 TFs pela atividade (dES) entre tumores,
#   feito SEPARADAMENTE em cada coorte.
#   Especificação principal: distância 1 - Spearman, ligação média, k = 2.
#   Sensibilidade: 1 - Spearman/complete, Manhattan/ward.D2, Euclidiana/ward.D2.
#
# ENTRADAS
#   TCGA: rtni_TF.RData (rede; atividade recalculada com tni.gsea2), 1 tumor primário por paciente (501)
#   ICGC: dES_ICGC_RT91.csv (91 tumores RT; gerado pelo R/07_survival_clinical/21_survival_ICGC.R)
#         + donor_data_final (ICGC_DEA_Data.RData) para tecido, doador e estágio M
#
# SAÍDAS (OUT_DIR)
#   FigS5_cluster_activity.(tiff|pdf|png)
#     A heatmap da atividade no ICGC-RECA | B, C correlação de Spearman entre os 14 TFs (TCGA, ICGC)
#     D escores dos clusters em M0 x M1 no ICGC-RECA (Cliff's delta; ORs na Table S28, script 11b)
#   TableS27_cluster_assignments.xlsx (abas Assignments, Sensitivity, Summary) + CSVs
#
# Ordem: rodar o 21_survival_ICGC.R antes (gera dES_ICGC_RT91.csv).
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(cowplot)
  library(openxlsx)
  library(RTN)
  library(paletteer)
})

set.seed(2026)
OUT_DIR <- path.expand("~/Kidney/Figures/FigS5_cluster_activity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ============================================================ 0. Definições
TF14 <- data.frame(
  symbol  = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A",
              "IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  cluster = rep(c("C1", "C2"), c(5, 9)))
CL_NAME <- c(C1 = "Proliferative", C2 = "Differentiation-associated")
CL_COL  <- c(C1 = "#542788", C2 = "#1a9850")
ref     <- setNames(TF14$cluster, TF14$symbol)       # network clusters (TCGA, Jaccard)

# Adjusted Rand index
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  s <- sum(choose(tab, 2))
  a <- sum(choose(rowSums(tab), 2)); b <- sum(choose(colSums(tab), 2))
  e <- a * b / choose(n, 2)
  (s - e) / ((a + b) / 2 - e)
}

# Métodos de agrupamento dos TFs (dif = amostras x TFs)
METODOS <- list(
  "1-Spearman / average (primary)" = function(d) hclust(as.dist(1 - cor(d, method = "spearman")), "average"),
  "1-Spearman / complete"          = function(d) hclust(as.dist(1 - cor(d, method = "spearman")), "complete"),
  "Manhattan / ward.D2"            = function(d) hclust(dist(t(d), "manhattan"), "ward.D2"),
  "Euclidean / ward.D2"            = function(d) hclust(dist(t(d)), "ward.D2"))

# Corta em k = 2 e rotula os grupos (C1/C2) pelo alinhamento de maior concordância
activity_groups <- function(hc) {
  g <- cutree(hc, k = 2)[names(ref)]
  a <- ifelse(g == 1, "C1", "C2"); b <- ifelse(g == 1, "C2", "C1")
  setNames(if (sum(a == ref) >= sum(b == ref)) a else b, names(ref))
}

## ============================================================ 1. Atividade e checagens
ENS <- c(FOXM1 = "ENSG00000111206", MYBL2 = "ENSG00000101057", NFE2L3 = "ENSG00000050344",
         PARP12 = "ENSG00000059378", TFAP2A = "ENSG00000137203", IRF6 = "ENSG00000117595",
         MEF2D = "ENSG00000116604", NR3C2 = "ENSG00000151623", PBX1 = "ENSG00000185630",
         RARB = "ENSG00000077092", RORC = "ENSG00000143365", SATB1 = "ENSG00000182568",
         TAL1 = "ENSG00000162367", TCF4 = "ENSG00000196628")
to_symbol <- function(m) {
  key <- sub("\\..*", "", colnames(m))
  colnames(m) <- ifelse(key %in% ENS, names(ENS)[match(key, ENS)], colnames(m)); m
}

# TCGA: atividade dos 14 regulons, 1 tumor primário (01) por paciente -> 501
load(path.expand("~/Kidney/0_data/rtni_TF.RData"))                       # rtni_TF
rtni_TF  <- tni.gsea2(rtni_TF, regulatoryElements = unname(ENS))
tcga_dif <- to_symbol(tni.get(rtni_TF, what = "regulonActivity")$differential)
tcga_dif <- tcga_dif[substr(rownames(tcga_dif), 14, 15) == "01", ]
tcga_dif <- tcga_dif[!duplicated(substr(rownames(tcga_dif), 1, 12)), ]
stopifnot(nrow(tcga_dif) == 501)

# ICGC: atividade única dos 91 tumores RT (21_survival_ICGC.R)
icgc_dif <- to_symbol(as.matrix(read.csv(path.expand("~/Kidney/0_data/dES_ICGC_RT91.csv"),
                                         row.names = 1, check.names = FALSE)))
load(path.expand("~/Kidney/0_data/ICGC_DEA_Data.RData"))                  # donor_data_final
icgc_annot <- as.data.frame(donor_data_final)
rownames(icgc_annot) <- paste0(icgc_annot$submitted_donor_id, icgc_annot$tissue_type)

act <- list(TCGA = tcga_dif[, TF14$symbol], ICGC = icgc_dif[, TF14$symbol])

# ICGC: somente tumores primários (RT), 1 por doador
icgc_s <- rownames(act$ICGC)
stopifnot(all(icgc_s %in% rownames(icgc_annot)),
          all(icgc_annot[icgc_s, "tissue_type"] == "RT"),
          !anyDuplicated(icgc_annot[icgc_s, "submitted_donor_id"]))
cat(sprintf("TCGA-KIRC: %d tumours | ICGC-RECA: %d tumours\n", nrow(act$TCGA), nrow(act$ICGC)))

## ============================================================ 2. Agrupamento principal + estatísticas
prim <- lapply(act, function(d) {
  hc  <- METODOS[[1]](d)
  grp <- activity_groups(hc)
  rho <- cor(d, method = "spearman")
  ct  <- table(network = ref, activity = grp)
  perm <- replicate(10000, ari(sample(ref), grp))
  ut   <- upper.tri(rho); same <- outer(ref, ref, "==")[ut]
  list(hc = hc, grp = grp, rho = rho, ct = ct,
       ARI      = ari(ref, grp),
       fisher_p = fisher.test(ct)$p.value,
       perm_p   = (sum(perm >= ari(ref, grp)) + 1) / (length(perm) + 1),
       rho_within  = mean(rho[ut][same]),
       rho_between = mean(rho[ut][!same]))
})

for (co in names(prim)) {
  cat(sprintf("\n%s: ARI = %.2f | Fisher p = %.2g | permutation p = %.2g | mean rho within = %.2f, between = %.2f\n",
              co, prim[[co]]$ARI, prim[[co]]$fisher_p, prim[[co]]$perm_p,
              prim[[co]]$rho_within, prim[[co]]$rho_between))
  print(prim[[co]]$ct)
}

## ============================================================ 3. Supplementary Table S27
# Aba 1 - atribuição de cada TF + coesão (rho médio com o próprio cluster vs o oposto)
tf_profile <- function(rho, tf) {
  own   <- setdiff(names(ref)[ref == ref[tf]], tf)
  other <- names(ref)[ref != ref[tf]]
  c(own = mean(rho[tf, own]), opposite = mean(rho[tf, other]))
}

tab_tf <- TF14 |>
  transmute(TF = symbol,
            `Network cluster (TCGA, target overlap)` = paste0(cluster, " - ", CL_NAME[cluster]))
for (co in names(prim)) {
  p  <- prim[[co]]
  pr <- t(sapply(TF14$symbol, function(tf) tf_profile(p$rho, tf)))
  tab_tf[[paste0(co, ": activity group")]]      <- p$grp[TF14$symbol]
  tab_tf[[paste0(co, ": concordant")]]          <- ifelse(p$grp[TF14$symbol] == ref, "Yes", "No")
  tab_tf[[paste0(co, ": mean rho, own cluster")]] <- round(pr[, "own"], 2)
  tab_tf[[paste0(co, ": mean rho, opposite")]]  <- round(pr[, "opposite"], 2)
}

# Aba 2 - sensibilidade (método x coorte)
tab_sens <- do.call(rbind, lapply(names(act), function(co) {
  do.call(rbind, lapply(names(METODOS), function(m) {
    grp  <- activity_groups(METODOS[[m]](act[[co]]))
    disc <- names(ref)[grp != ref]
    data.frame(Cohort = co, Method = m, ARI = round(ari(ref, grp), 2),
               `Discordant TFs` = ifelse(length(disc), paste(disc, collapse = ", "), "none"),
               check.names = FALSE)
  }))
}))
print(tab_sens, row.names = FALSE)

# Aba 3 - resumo (especificação principal)
tab_sum <- do.call(rbind, lapply(names(prim), function(co) {
  p <- prim[[co]]
  data.frame(Cohort = co, `Tumours (n)` = nrow(act[[co]]),
             ARI = round(p$ARI, 2), `Fisher exact p` = signif(p$fisher_p, 2),
             `Permutation p (10,000)` = signif(p$perm_p, 2),
             `Mean rho within clusters` = round(p$rho_within, 2),
             `Mean rho between clusters` = round(p$rho_between, 2), check.names = FALSE)
}))

wb <- createWorkbook()
for (sh in list(list("Assignments", tab_tf), list("Sensitivity", tab_sens), list("Summary", tab_sum))) {
  addWorksheet(wb, sh[[1]]); writeData(wb, sh[[1]], sh[[2]])
  setColWidths(wb, sh[[1]], cols = seq_len(ncol(sh[[2]])), widths = "auto")
}
saveWorkbook(wb, file.path(OUT_DIR, "TableS27_cluster_assignments.xlsx"), overwrite = TRUE)
write.csv(tab_tf,   file.path(OUT_DIR, "TableS27_assignments.csv"), row.names = FALSE)
write.csv(tab_sens, file.path(OUT_DIR, "TableS27_sensitivity.csv"), row.names = FALSE)
write.csv(tab_sum,  file.path(OUT_DIR, "TableS27_summary.csv"),     row.names = FALSE)

## ============================================================ 4. Supplementary Figure S5
ann_row <- data.frame(Cluster = ref, row.names = names(ref))
m_lab <- function(x) { x <- as.character(x); ifelse(x %in% c("M0", "M1"), x, "MX/NA") }
icgc_M <- m_lab(icgc_annot[rownames(act$ICGC), "pathologic_m"])
ann_col <- data.frame(M = icgc_M, row.names = rownames(act$ICGC))
ann_colors <- list(Cluster = CL_COL, M = c(M0 = "#2b6cb0", M1 = "#d73027", "MX/NA" = "grey85"))
my_colors <- as.character(paletteer_d("MexBrewer::Alacena"))

# A - heatmap do ICGC (o do TCGA é a Figura 4B)
hmA <- pheatmap(t(act$ICGC),
                main = sprintf("ICGC-RECA (n = %d tumours) - ARI = %.2f", nrow(act$ICGC), prim$ICGC$ARI),
                cluster_rows = prim$ICGC$hc, cutree_rows = 2,
                clustering_distance_cols = "euclidean", clustering_method = "average",
                annotation_row = ann_row, annotation_col = ann_col,
                annotation_colors = ann_colors, annotation_names_row = FALSE,
                color = my_colors, border_color = NA,
                show_colnames = FALSE, fontsize = 7, fontsize_row = 8, silent = TRUE)$gtable

# B, C - correlação de Spearman entre os TFs (sem números: valores na Table S27 e nos CSVs)
cor_plot <- function(co, n) {
  rho <- prim[[co]]$rho
  ord <- prim[[co]]$hc$labels[prim[[co]]$hc$order]
  write.csv(round(rho, 3), file.path(OUT_DIR, paste0(co, "_TF_correlation.csv")))
  df  <- as.data.frame(as.table(rho)) |> setNames(c("a", "b", "rho")) |>
    mutate(a = factor(a, ord), b = factor(b, rev(ord)))
  ggplot(df, aes(a, b, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-1, 1), name = expression(rho)) +
    labs(x = NULL, y = NULL, title = sprintf("%s (n = %d)", c(TCGA = "TCGA-KIRC", ICGC = "ICGC-RECA")[co], n)) +
    coord_fixed() + theme_minimal(base_size = 7) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = suppressWarnings(element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                      color = CL_COL[ref[ord]])),
          axis.text.y = suppressWarnings(element_text(color = CL_COL[ref[rev(ord)]])),
          legend.key.width = unit(2.5, "mm"))
}

# D - escores dos clusters (média dos dES) em M0 x M1 no ICGC-RECA
cliff <- function(x, y) (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
dD <- data.frame(M = icgc_M,
                 `Cluster 1` = rowMeans(act$ICGC[, names(ref)[ref == "C1"]]),
                 `Cluster 2` = rowMeans(act$ICGC[, names(ref)[ref == "C2"]]), check.names = FALSE) |>
  filter(M %in% c("M0", "M1"))
longD <- bind_rows(lapply(c("Cluster 1", "Cluster 2"), function(cl)
  data.frame(cluster = cl, M = dD$M, score = dD[[cl]])))
annD <- bind_rows(lapply(c("Cluster 1", "Cluster 2"), function(cl) {
  x <- dD[[cl]][dD$M == "M1"]; y <- dD[[cl]][dD$M == "M0"]
  b <- replicate(2000, cliff(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  data.frame(cluster = cl, txt = sprintf("δ = %.2f (%.2f to %.2f)", cliff(x, y),
                                         quantile(b, 0.025), quantile(b, 0.975)))
}))
pD <- ggplot(longD, aes(M, score, fill = cluster)) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.45, color = "#8a8983") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, color = "#0b0b0b", linewidth = 0.3) +
  geom_text(data = annD, aes(x = 1.5, y = Inf, label = txt), vjust = 1.3, size = 2.2,
            color = "#52514e", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Cluster 1" = unname(CL_COL["C1"]), "Cluster 2" = unname(CL_COL["C2"])),
                    guide = "none") +
  scale_x_discrete(labels = c(M0 = sprintf("M0\n(n = %d)", sum(dD$M == "M0")),
                              M1 = sprintf("M1\n(n = %d)", sum(dD$M == "M1")))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  facet_wrap(~cluster, scales = "free_y") +
  labs(x = NULL, y = "Cluster activity (mean dES)", title = "ICGC-RECA") +
  theme_minimal(base_size = 7) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))

fig <- plot_grid(
  hmA,
  plot_grid(cor_plot("TCGA", nrow(act$TCGA)), cor_plot("ICGC", nrow(act$ICGC)), pD,
            ncol = 3, rel_widths = c(1, 1, 1.1), labels = c("B", "C", "D"), label_size = 10),
  ncol = 1, rel_heights = c(1, 1), labels = c("A", ""), label_size = 10)

W <- 190; H <- 190
ggsave(file.path(OUT_DIR, "FigS5_cluster_activity.tiff"), fig, width = W, height = H, units = "mm",
       dpi = 500, compression = "lzw", bg = "white")
ggsave(file.path(OUT_DIR, "FigS5_cluster_activity.pdf"), fig, width = W, height = H, units = "mm",
       device = cairo_pdf)
ggsave(file.path(OUT_DIR, "FigS5_cluster_activity.png"), fig, width = W, height = H, units = "mm",
       dpi = 300, bg = "white")

cat("\nArquivos salvos em:", OUT_DIR, "\n")
