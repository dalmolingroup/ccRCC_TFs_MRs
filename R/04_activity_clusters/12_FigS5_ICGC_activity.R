### Figura suplementar: atividade dos 14 regulons no ICGC-RECA e concordância com os clusters da rede ###
# Revisão CSBJ-D-26-01039 (R2.M5, R2.M2; Resultados 3.3)
#
# Os clusters foram definidos no TCGA pela sobreposição dos alvos (Jaccard, amap). Aqui se verifica se,
# no ICGC-RECA (não usado para defini-los), os 14 TFs se agrupam da mesma forma pela ATIVIDADE (dES).
#   1. agrupamento hierárquico dos TFs pela correlação de Spearman da atividade entre os tumores
#      (distância 1 - rho, ligação média), cortado em k = 2
#   2. concordância com os clusters da rede: índice de Rand ajustado (ARI)
#
# Painéis:  A  heatmap da atividade (TFs x tumores; barra = cluster da rede; tumores ordenados por M)
#           B  matriz de correlação de Spearman entre os 14 TFs
#
# Saídas (OUT_DIR): FigS_ICGC_activity.(tiff|png|pdf), ICGC_TF_grouping.csv, ICGC_TF_correlation.csv
# Pacotes: ggplot2, dplyr, cowplot, patchwork (não precisa de pheatmap/ComplexHeatmap).

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(cowplot); library(patchwork) })

## ---------------------------------------------------------------- CONFIG
TNS_RDATA   <- "~/Kidney/0_data/rtni_TF_survival_ICGC_Boruta_Newsignature.RData"
TNS_OBJECT  <- "rtns_TF_ICGC"
META_FILE   <- "~/Kidney/0_data/metadata_ICGC_full.csv"
TUMOR_COL   <- "tissue_type"; TUMOR_VALUE <- "RT"
DONOR_COL   <- "submitted_donor_id"
M_COL       <- "pathologic_m"
OUT_DIR     <- "~/Kidney/Figures/FigS_ICGC_activity"
ACTIVITY_CSV <- NULL       # alternativa ao TNS: CSV amostras x regulons (dES), 1ª coluna = amostra

LINKAGE <- "average"       # o ARI = 0,12 da carta corresponde a este agrupamento; confira na saída
K       <- 2

# Cluster 1 = proliferativo; Cluster 2 = diferenciação (nomes só a partir da seção 3.4)
TF14 <- data.frame(
  symbol  = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A",
              "IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  ensembl = c("ENSG00000111206", "ENSG00000101057", "ENSG00000050344", "ENSG00000059378", "ENSG00000137203",
              "ENSG00000117595", "ENSG00000116604", "ENSG00000151623", "ENSG00000185630", "ENSG00000077092",
              "ENSG00000143365", "ENSG00000182568", "ENSG00000162367", "ENSG00000196628"),
  cluster = rep(c("Cluster 1", "Cluster 2"), c(5, 9)))
CL_COL <- c("Cluster 1" = "#542788", "Cluster 2" = "#1a9850")

pe <- path.expand
dir.create(pe(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
strip_v <- function(x) sub("\\.\\d+$", "", x)

## ---------------------------------------------------------------- 1. atividade
if (is.null(ACTIVITY_CSV)) {
  suppressPackageStartupMessages(library(RTNsurvival))
  env <- new.env(); load(pe(TNS_RDATA), envir = env)
  tns <- get(TNS_OBJECT, envir = env)
  act <- tryCatch(tnsGet(tns, "regulonActivity"), error = function(e) NULL)
  reg <- tryCatch(tnsGet(tns, "regulatoryElements"),
                  error = function(e) tryCatch(tni.get(tnsGet(tns, "TNI"), "regulatoryElements"),
                                               error = function(e) NULL))
} else {
  act <- list(dif = as.matrix(read.csv(pe(ACTIVITY_CSV), row.names = 1, check.names = FALSE)))
  reg <- NULL
}

# nome do regulon no objeto -> símbolo (aceita símbolo, ENSG com/sem versão ou o vetor regulatoryElements)
resolve <- function(cols) {
  key <- strip_v(cols)
  if (!is.null(reg)) {                               # regulatoryElements: nomes e valores = símbolo/ID
    alt <- c(setNames(strip_v(unname(reg)), names(reg)), setNames(names(reg), strip_v(unname(reg))))
    key2 <- ifelse(is.na(alt[cols]), key, alt[cols])
  } else key2 <- key
  sym <- ifelse(key %in% TF14$symbol, key,
         ifelse(key %in% TF14$ensembl, TF14$symbol[match(key, TF14$ensembl)],
         ifelse(key2 %in% TF14$symbol, key2,
         ifelse(strip_v(key2) %in% TF14$ensembl, TF14$symbol[match(strip_v(key2), TF14$ensembl)], NA))))
  setNames(sym, cols)
}
need_gsea <- is.null(act$dif) || sum(!is.na(resolve(colnames(act$dif)))) < nrow(TF14)
if (need_gsea) {
  if (!is.null(ACTIVITY_CSV)) stop("O CSV não tem os 14 TFs")
  if (is.null(reg)) stop("Não consegui ler os regulatoryElements do TNS")
  cand <- unique(c(names(reg), unname(reg)))
  ids <- cand[!is.na(resolve(cand))]
  message("Calculando a atividade (tnsGSEA2) para ", length(ids), " TFs...")
  tns <- tnsGSEA2(tns, regulatoryElements = ids, verbose = FALSE)
  act <- tnsGet(tns, "regulonActivity")
}
dES <- act$dif
sym <- resolve(colnames(dES))
dES <- dES[, !is.na(sym), drop = FALSE]; colnames(dES) <- sym[!is.na(sym)]
miss <- setdiff(TF14$symbol, colnames(dES))
if (length(miss)) stop("TFs sem atividade no ICGC: ", paste(miss, collapse = ", "))
dES <- dES[, TF14$symbol]

## ---------------------------------------------------------------- 2. tumores (1 por doador) e status M
meta <- if (file.exists(pe(META_FILE))) read.csv(pe(META_FILE), check.names = FALSE, row.names = 1) else NULL
samp <- data.frame(sample = rownames(dES))
if (!is.null(meta)) {
  meta$sample <- rownames(meta)
  samp <- left_join(samp, meta, by = "sample")
  if (TUMOR_COL %in% names(samp)) samp <- samp[is.na(samp[[TUMOR_COL]]) | samp[[TUMOR_COL]] == TUMOR_VALUE, ]
  if (DONOR_COL %in% names(samp)) samp <- samp[is.na(samp[[DONOR_COL]]) | !duplicated(samp[[DONOR_COL]]), ]
}
samp$M <- if (M_COL %in% names(samp)) ifelse(samp[[M_COL]] %in% c("M0", "M1"), samp[[M_COL]], "NA") else "NA"
dES <- dES[samp$sample, ]
cat(sprintf("ICGC-RECA: %d tumores (M0 = %d, M1 = %d, M desconhecido = %d)\n", nrow(dES),
            sum(samp$M == "M0"), sum(samp$M == "M1"), sum(samp$M == "NA")))

## ---------------------------------------------------------------- 3. agrupamento dos TFs e ARI
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  s <- sum(choose(tab, 2)); a <- sum(choose(rowSums(tab), 2)); b <- sum(choose(colSums(tab), 2))
  e <- a * b / choose(n, 2); (s - e) / ((a + b) / 2 - e)
}
rho <- cor(dES, method = "spearman")
hc  <- hclust(as.dist(1 - rho), method = LINKAGE)
grp <- cutree(hc, k = K)
# nomeia os grupos de atividade pelo cluster da rede mais frequente em cada um (só para a tabela)
maj <- tapply(TF14$cluster, grp[TF14$symbol], function(z) names(which.max(table(z))))
tf_tab <- TF14 |> mutate(activity_group = paste0("Group ", grp[symbol]),
                         activity_group_majority = maj[as.character(grp[symbol])],
                         same_as_network = activity_group_majority == cluster)
ARI <- ari(tf_tab$cluster, tf_tab$activity_group)
write.csv(tf_tab, file.path(pe(OUT_DIR), "ICGC_TF_grouping.csv"), row.names = FALSE)
write.csv(round(rho, 3), file.path(pe(OUT_DIR), "ICGC_TF_correlation.csv"))

# correlação média dentro e entre clusters da rede (frase opcional para o texto)
cl <- setNames(TF14$cluster, TF14$symbol); ut <- upper.tri(rho)
same <- outer(cl, cl, "==")[ut]
cat(sprintf("rho médio dentro dos clusters = %.2f | entre clusters = %.2f\n", mean(rho[ut][same]), mean(rho[ut][!same])))

## ---------------------------------------------------------------- 4. figura
ord_tf <- hc$labels[hc$order]
z <- scale(dES)                                     # só para cor (dES já é relativo à média da coorte)
z[z > 2.5] <- 2.5; z[z < -2.5] <- -2.5
# tumores: M0 e depois M1; dentro de cada, pela diferença média C1 - C2
score <- rowMeans(z[, TF14$symbol[TF14$cluster == "Cluster 1"]]) - rowMeans(z[, TF14$symbol[TF14$cluster == "Cluster 2"]])
ord_s <- samp$sample[order(factor(samp$M, levels = c("M0", "M1", "NA")), score)]

hm <- as.data.frame(as.table(z)) |> setNames(c("sample", "tf", "z")) |>
  mutate(sample = factor(sample, ord_s), tf = factor(tf, rev(ord_tf)))
bar_tf <- data.frame(tf = factor(ord_tf, rev(ord_tf)), cluster = cl[ord_tf])
bar_m  <- data.frame(sample = factor(ord_s, ord_s), M = samp$M[match(ord_s, samp$sample)])

th <- theme_minimal(base_size = 7) + theme(panel.grid = element_blank())
pA_main <- ggplot(hm, aes(sample, tf, fill = z)) + geom_tile() +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-2.5, 2.5),
                       name = "Activity\n(z)") +
  labs(x = sprintf("ICGC-RECA tumours (n = %d)", nrow(dES)), y = NULL) +
  th + theme(axis.text.x = element_blank(), legend.key.width = unit(2.5, "mm"))
pA_cl <- ggplot(bar_tf, aes(1, tf, fill = cluster)) + geom_tile() +
  scale_fill_manual(values = CL_COL, name = "Network cluster") + theme_void() +
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 7))
pA_m <- ggplot(bar_m, aes(sample, 1, fill = M)) + geom_tile() +
  scale_fill_manual(values = c(M0 = "#2b6cb0", M1 = "#d73027", "NA" = "grey85"), name = "M stage") + theme_void() +
  theme(legend.text = element_text(size = 6), legend.title = element_text(size = 7))
pA <- (pA_m + theme(plot.margin = margin(0, 0, 1, 0))) + plot_spacer() + pA_main + pA_cl +
  plot_layout(ncol = 2, widths = c(1, 0.035), heights = c(0.045, 1), guides = "collect") &
  theme(legend.position = "right", legend.key.size = unit(3, "mm"))
pA <- patchworkGrob(pA)            # vira grob para entrar no plot_grid do cowplot

cm <- as.data.frame(as.table(rho)) |> setNames(c("a", "b", "rho")) |>
  mutate(a = factor(a, ord_tf), b = factor(b, rev(ord_tf)))
pB <- ggplot(cm, aes(a, b, fill = rho)) + geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", rho)), size = 1.7) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = expression(rho)) +
  labs(x = NULL, y = NULL, subtitle = sprintf("Adjusted Rand index vs network clusters = %.2f", ARI)) +
  th + theme(axis.text.x = suppressWarnings(element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                         color = CL_COL[cl[ord_tf]])),
             axis.text.y = suppressWarnings(element_text(color = CL_COL[cl[rev(ord_tf)]])),
             legend.key.width = unit(2.5, "mm")) + coord_fixed()

fig <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1.35, 1), labels = c("A", "B"), label_size = 10)
W <- 190; H <- 95
ggsave(file.path(pe(OUT_DIR), "FigS_ICGC_activity.tiff"), fig, width = W, height = H, units = "mm",
       dpi = 500, compression = "lzw", bg = "white")
ggsave(file.path(pe(OUT_DIR), "FigS_ICGC_activity.png"), fig, width = W, height = H, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(pe(OUT_DIR), "FigS_ICGC_activity.pdf"), fig, width = W, height = H, units = "mm", device = cairo_pdf)

## ---------------------------------------------------------------- 5. resumo
cat(sprintf("\nAgrupamento pela atividade (1 - Spearman, ligação %s, k = %d): ARI = %.2f\n", LINKAGE, K, ARI))
print(tf_tab |> select(symbol, cluster, activity_group, same_as_network), row.names = FALSE)
cat("\nTFs no grupo oposto ao da rede: ",
    paste(tf_tab$symbol[!tf_tab$same_as_network], collapse = ", "), "\n")
cat("(Esperado pela carta: TFAP2A, MEF2D, TAL1, TCF4; ARI = 0,12)\n")
