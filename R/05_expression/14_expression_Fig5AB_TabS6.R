###########################################################################
##### Figura 5 — Expressão dos MRs (tumor vs normal e M1 vs M0)       #####
##### Autor: Farias, Epitácio — revisão CSBJ-D-26-01039 (out/2026)    #####
###########################################################################
# Junta os scripts 9_Expression_Evaluation_MR_Boruta_Sign.R e
# 9_Expression_Evaluation_MR_Boruta_Sign_in_Metastasis_Samples.R e corrige:
#   1. M1 vs M0 só com amostras TUMORAIS (antes as normais entravam com o M do paciente)
#   2. Um tumor primário por paciente; no TCGA, a mesma coorte do ML (422 M0 / 79 M1)
#   3. Tumor vs normal com teste PAREADO (Wilcoxon signed-rank) nos pacientes com as duas
#      amostras (R2, Fig 5); o teste não pareado com todas as amostras vai para a tabela
#   4. Tamanho de efeito, n por grupo e p ajustado (BH entre os TFs) numa tabela suplementar
#      (R2 Fig 5, R3.15)
#   5. Clusters definidos explicitamente (sem inversão manual) e cores consistentes
#   6. Rascunho de legenda gerado com os n reais
# Também salva 0_data/ensembl_to_symbol.csv (para os eixos das outras figuras).
#
# Entrada: os .RData que os scripts antigos já salvavam (o download do GDC não é refeito).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG
TCGA_RDATA   <- "~/Kidney/0_data/KIRC_MR_expression_eval_data.RData"   # count_data_kirc2, gene_library2, clin_data
ICGC_RDATA   <- "~/Kidney/0_data/ICGC_DEA_Data.RData"                  # count_final_ICGC, donor_data_final
ML_COHORT    <- "~/Kidney/0_data/metadata_KIRC_full.csv"               # coorte do ML (opcional; define M e pacientes)
ML_ID_COL    <- "sample"
SIGNATURE    <- "~/Kidney/results_ml_revision/importancia_features_modelo_final_class_weight.csv"           # um ENSG por linha
OUT_DIR      <- "~/Kidney/Figures/Figure5"
SYMBOL_OUT   <- "~/Kidney/0_data/ensembl_to_symbol.csv"

# Clusters por símbolo (a ordem aqui é a ordem dos painéis). Ajuste se necessário.
CLUSTERS <- list(
  "cluster_2"        = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  "cluster_1" = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A")
)

# Use as MESMAS cores em todos os painéis da Fig 5 (inclusive nos KMs feitos no Inkscape)
CLUSTER_COLORS <- c("cluster_1" = "#542788", "cluster_2" = "#1a9850")
GROUP_COLORS   <- c(Normal = "#bdbdbd", Tumor = "#c51b7d", M0 = "#2166ac", M1 = "#d73027")

ICGC_TUMOR <- "RT"; ICGC_NORMAL <- "RA"       # códigos de tissue_type no ICGC-RECA
ICGC_HISTOLOGY_FILTER <- NULL                 # ex.: c("Clear cell renal cell carcinoma") se houver coluna de histologia

dir.create(path.expand(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
pe <- path.expand

## ---------------------------------------------------------------- utilitários
log_norm <- function(counts) {
  # mesma normalização dos scripts antigos: log2(contagem normalizada por DESeq2 + 0.5)
  counts <- as.matrix(counts); storage.mode(counts) <- "numeric"
  sf <- DESeq2::estimateSizeFactorsForMatrix(round(counts))
  log2(sweep(counts, 2, sf, "/") + 0.5)
}

cliffs_delta <- function(x, y) {             # P(x > y) - P(x < y), a partir do U de Mann-Whitney
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))$statistic
  unname(2 * w / (length(x) * length(y)) - 1)
}
paired_rb <- function(d) {                   # correlação rank-biserial pareada
  d <- d[d != 0]; n <- length(d)
  if (n == 0) return(NA_real_)
  v <- suppressWarnings(wilcox.test(d, exact = FALSE))$statistic
  unname(2 * v / (n * (n + 1) / 2) - 1)
}
stars <- function(p) cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), c("****", "***", "**", "*", "ns"))

## ---------------------------------------------------------------- assinatura e símbolos
env_t <- new.env(); load(pe(TCGA_RDATA), envir = env_t)
gene_lib <- env_t$gene_library2
gene_lib$gene_id <- sub("\\..*", "", gene_lib$gene_id)
write.csv(data.frame(ensembl = gene_lib$gene_id, symbol = gene_lib$gene_name),
          pe(SYMBOL_OUT), row.names = FALSE)

sig <- vroom::vroom(SIGNATURE)
sig <- sig$...1
sym <- setNames(gene_lib$gene_name, gene_lib$gene_id)[sig]
if (any(is.na(sym))) stop("Sem símbolo para: ", paste(sig[is.na(sym)], collapse = ", "))
cluster_of <- setNames(rep(names(CLUSTERS), lengths(CLUSTERS)), unlist(CLUSTERS))
if (length(setdiff(sym, names(cluster_of))))
  stop("TFs da assinatura sem cluster definido em CLUSTERS: ", paste(setdiff(sym, names(cluster_of)), collapse = ", "))
gene_order <- unlist(CLUSTERS)[unlist(CLUSTERS) %in% sym]
message("Assinatura: ", length(sig), " TFs — ", paste(gene_order, collapse = ", "))

## ---------------------------------------------------------------- TCGA
cnt <- env_t$count_data_kirc2
cnt <- cnt[, vapply(cnt, is.numeric, logical(1)), drop = FALSE]     # tira colunas de anotação, se houver
clin <- env_t$clin_data
smp <- intersect(colnames(cnt), rownames(clin))
clin <- clin[smp, ]
expr_t <- log_norm(cnt[, smp])[sig, , drop = FALSE]

ann_t <- data.frame(sample = smp, patient = substr(smp, 1, 12),
                    type = ifelse(clin$sample_type == "Primary Tumor", "Tumor",
                                  ifelse(clin$sample_type == "Solid Tissue Normal", "Normal", NA)),
                    M = as.character(clin$ajcc_pathologic_m)) |>
  filter(!is.na(type))

# coorte do ML: define pacientes e M dos tumores (mesmos 422/79 do classificador)
if (file.exists(pe(ML_COHORT))) {
  ml <- read.csv(pe(ML_COHORT), check.names = FALSE)
  ml <- data.frame(ml_sample = ml[[ML_ID_COL]], patient = substr(ml[[ML_ID_COL]], 1, 12),
                   M_ml = ml$ajcc_pathologic_m)
  ann_t <- left_join(ann_t, distinct(ml, patient, .keep_all = TRUE), by = "patient")
  ann_t$M <- ifelse(!is.na(ann_t$M_ml), ann_t$M_ml, NA)
  ann_t$pref <- ann_t$sample == ann_t$ml_sample
} else {
  message("[aviso] coorte do ML não encontrada; usando ajcc_pathologic_m do colData")
  ann_t$pref <- FALSE
}
pick_one <- function(a) a |> arrange(patient, desc(pref), sample) |> distinct(patient, .keep_all = TRUE)
tum_t <- pick_one(filter(ann_t, type == "Tumor"))
nor_t <- pick_one(filter(ann_t, type == "Normal"))

## ---------------------------------------------------------------- ICGC
env_i <- new.env(); load(pe(ICGC_RDATA), envir = env_i)
ci <- env_i$count_final_ICGC
if ("gene_id" %in% names(ci)) { rownames(ci) <- sub("\\..*", "", ci$gene_id); ci$gene_id <- NULL }
ci <- ci[, vapply(ci, is.numeric, logical(1)), drop = FALSE]
don <- env_i$donor_data_final
hist_cols <- grep("histolog", names(don), ignore.case = TRUE, value = TRUE)
if (length(hist_cols)) { message("ICGC — histologia (R1.m1):"); print(table(don[[hist_cols[1]]], useNA = "ifany")) }
if (!is.null(ICGC_HISTOLOGY_FILTER) && length(hist_cols)) don <- don[don[[hist_cols[1]]] %in% ICGC_HISTOLOGY_FILTER, ]
don$sample <- paste0(don$submitted_donor_id, don$tissue_type)
don <- don[don$sample %in% colnames(ci), ]
sig_i <- intersect(sig, rownames(ci))
if (length(setdiff(sig, sig_i))) message("[ICGC] ausentes: ", paste(sym[setdiff(sig, sig_i)], collapse = ", "))
expr_i <- log_norm(ci[, don$sample])[sig_i, , drop = FALSE]
ann_i <- data.frame(sample = don$sample, patient = as.character(don$submitted_donor_id),
                    type = ifelse(don$tissue_type == ICGC_TUMOR, "Tumor",
                                  ifelse(don$tissue_type == ICGC_NORMAL, "Normal", NA)),
                    M = as.character(don$pathologic_m), pref = FALSE) |> filter(!is.na(type))
tum_i <- pick_one(filter(ann_i, type == "Tumor"))
nor_i <- pick_one(filter(ann_i, type == "Normal"))

## ---------------------------------------------------------------- dados longos
long <- function(expr, ann, cohort) {
  as.data.frame(t(expr[, ann$sample, drop = FALSE])) |>
    mutate(sample = ann$sample) |>
    pivot_longer(-sample, names_to = "ensembl", values_to = "expression") |>
    left_join(ann, by = "sample") |>
    mutate(cohort = cohort, gene = unname(sym[ensembl]))
}
tn <- list(
  `TCGA-KIRC` = list(t = tum_t, n = nor_t, e = expr_t),
  `ICGC-RECA` = list(t = tum_i, n = nor_i, e = expr_i)
)

## ---------------------------------------------------------------- testes
stats_tn <- list(); stats_m <- list(); plot_tn <- list(); plot_m <- list(); counts <- list()
for (co in names(tn)) {
  d <- tn[[co]]
  pairs <- intersect(d$t$patient, d$n$patient)
  tp <- d$t[match(pairs, d$t$patient), ]; np <- d$n[match(pairs, d$n$patient), ]
  tum_m <- filter(d$t, M %in% c("M0", "M1"))
  counts[[co]] <- c(pairs = length(pairs), tumor = nrow(d$t), normal = nrow(d$n),
                    M0 = sum(tum_m$M == "M0"), M1 = sum(tum_m$M == "M1"))
  message(co, ": ", paste(names(counts[[co]]), counts[[co]], sep = " = ", collapse = ", "))
  
  for (g in rownames(d$e)) {
    xt <- d$e[g, tp$sample]; xn <- d$e[g, np$sample]; dif <- xt - xn
    pw <- suppressWarnings(wilcox.test(xt, xn, paired = TRUE, exact = FALSE))$p.value
    at <- d$e[g, d$t$sample]; an <- d$e[g, d$n$sample]
    pu <- suppressWarnings(wilcox.test(at, an, exact = FALSE))$p.value
    stats_tn[[length(stats_tn) + 1]] <- data.frame(
      cohort = co, gene = sym[[g]], ensembl = g, cluster = cluster_of[[sym[[g]]]],
      n_pairs = length(pairs), median_normal = median(xn), median_tumor = median(xt),
      median_log2FC_paired = median(dif), rank_biserial_paired = paired_rb(dif), p_paired = pw,
      n_tumor_all = length(at), n_normal_all = length(an), cliffs_delta_all = cliffs_delta(at, an),
      p_unpaired_all = pu)
    
    x1 <- d$e[g, tum_m$sample[tum_m$M == "M1"]]; x0 <- d$e[g, tum_m$sample[tum_m$M == "M0"]]
    stats_m[[length(stats_m) + 1]] <- data.frame(
      cohort = co, gene = sym[[g]], ensembl = g, cluster = cluster_of[[sym[[g]]]],
      n_M0 = length(x0), n_M1 = length(x1), median_M0 = median(x0), median_M1 = median(x1),
      log2FC_median = median(x1) - median(x0), cliffs_delta = cliffs_delta(x1, x0),
      p = suppressWarnings(wilcox.test(x1, x0, exact = FALSE))$p.value)
  }
  plot_tn[[co]] <- long(d$e, bind_rows(tp, np), co)          # figura: só os pares (mesmos do teste)
  plot_m[[co]]  <- long(d$e, tum_m, co) |> mutate(type = M)
}
stats_tn <- bind_rows(stats_tn) |> group_by(cohort) |>
  mutate(p_adj_paired = p.adjust(p_paired, "BH"), p_adj_unpaired_all = p.adjust(p_unpaired_all, "BH"),
         signif = as.character(stars(p_adj_paired))) |> ungroup()
stats_m <- bind_rows(stats_m) |> group_by(cohort) |>
  mutate(p_adj = p.adjust(p, "BH"), signif = as.character(stars(p_adj))) |> ungroup()

write.csv(stats_tn, file.path(pe(OUT_DIR), "SuppTable_expression_tumor_vs_normal.csv"), row.names = FALSE)
write.csv(stats_m,  file.path(pe(OUT_DIR), "SuppTable_expression_M1_vs_M0.csv"), row.names = FALSE)

# dados da figura para a montagem final no formato da revista (15b_Figure5_layout.R)
saveRDS(list(plot_tn = bind_rows(plot_tn), plot_m = bind_rows(plot_m), stats_tn = stats_tn, stats_m = stats_m,
             counts = counts, gene_order = gene_order, cluster_of = cluster_of,
             CLUSTER_COLORS = CLUSTER_COLORS, GROUP_COLORS = GROUP_COLORS),
        file.path(pe(OUT_DIR), "fig5_data.rds"))

## ---------------------------------------------------------------- figuras
cohort_lab <- function(co, kind) {
  k <- counts[[co]]
  if (kind == "tn") sprintf("%s\n(%d pairs)", co, k["pairs"]) else sprintf("%s\n(%d M0 / %d M1)", co, k["M0"], k["M1"])
}

make_plot <- function(df, st, groups, title, kind) {
  df <- df |> mutate(cohort_lab = vapply(cohort, cohort_lab, "", kind = kind),
                     gene = factor(gene, levels = gene_order), type = factor(type, levels = groups))
  labs_order <- vapply(names(tn), cohort_lab, "", kind = kind)
  df$cohort_lab <- factor(df$cohort_lab, levels = labs_order)
  br <- df |> group_by(gene, cohort_lab) |>
    summarise(top = max(expression), rng = diff(range(expression)), .groups = "drop") |>
    mutate(y = top + 0.08 * rng, x = as.integer(cohort_lab)) |>
    left_join(st |> mutate(cohort_lab = factor(vapply(cohort, cohort_lab, "", kind = kind), levels = labs_order),
                           gene = factor(gene, levels = gene_order)) |>
                select(gene, cohort_lab, signif), by = c("gene", "cohort_lab"))
  p <- ggplot(df, aes(cohort_lab, expression, fill = type)) +
    geom_boxplot(outlier.size = 0.6, width = 0.6, position = position_dodge(0.8), linewidth = 0.35) +
    geom_segment(data = br, aes(x = x - 0.3, xend = x + 0.3, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = br, aes(x = x, y = y, label = signif), inherit.aes = FALSE, vjust = -0.3, size = 3) +
    scale_fill_manual(values = GROUP_COLORS[groups], name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(title = title, x = NULL, y = "Expression (log2 normalized counts)") +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(face = "bold", color = "white"),
          axis.text.x = element_text(size = 8), legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  strip_fill <- unname(CLUSTER_COLORS[cluster_of[gene_order]])
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    p + ggh4x::facet_wrap2(~gene, ncol = 3, scales = "free_y",
                           strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_fill)))
  } else {
    message("[aviso] ggh4x indisponível: faixas sem cor por cluster")
    p + facet_wrap(~gene, ncol = 3, scales = "free_y") +
      theme(strip.text = element_text(face = "bold", color = "black"))
  }
}

p_tn <- make_plot(bind_rows(plot_tn), stats_tn, c("Normal", "Tumor"),
                  "Master regulator expression: tumour vs matched normal tissue", "tn")
p_m  <- make_plot(bind_rows(plot_m), stats_m, c("M0", "M1"),
                  "Master regulator expression by metastatic status (tumours only)", "m")
h <- 1.1 + 1.9 * ceiling(length(gene_order) / 3)
for (nm in c("Figure5B_tumor_vs_normal", "Figure5C_M1_vs_M0")) {
  p <- if (grepl("normal", nm)) p_tn else p_m
  ggsave(file.path(pe(OUT_DIR), paste0(nm, ".tiff")), p, width = 9, height = h, dpi = 300, compression = "lzw")
  ggsave(file.path(pe(OUT_DIR), paste0(nm, ".pdf")), p, width = 9, height = h)
}

## ---------------------------------------------------------------- rascunho de legenda
color_name <- function(hex) {
  known <- c("#542788" = "purple", "#1a9850" = "green", "#2166ac" = "blue", "#d73027" = "red")
  ifelse(tolower(hex) %in% names(known), known[tolower(hex)], hex)
}
k1 <- counts[["TCGA-KIRC"]]; k2 <- counts[["ICGC-RECA"]]
legend <- sprintf(paste0(
  "Figure 5. Expression of the master regulators in TCGA-KIRC and ICGC-RECA. ",
  "(B) Tumour versus matched adjacent normal tissue (TCGA-KIRC, n = %d pairs; ICGC-RECA, n = %d pairs); ",
  "two-sided Wilcoxon signed-rank test. (C) Primary tumours from patients with (M1) and without (M0) distant ",
  "metastasis at diagnosis (TCGA-KIRC, %d M0 / %d M1; ICGC-RECA, %d M0 / %d M1); two-sided Mann-Whitney U test. ",
  "Expression is shown as log2 DESeq2-normalized counts. P values were adjusted across the regulators available in each ",
  "cohort (%d in TCGA-KIRC), separately for each comparison (Benjamini-Hochberg): ns, not significant; *p < 0.05; **p < 0.01; ***p < 0.001; ****p < 0.0001. ",
  "Strip colours indicate the regulatory cluster (%s). Effect sizes and exact values are given in Supplementary Tables SX and SY."),
  k1["pairs"], k2["pairs"], k1["M0"], k1["M1"], k2["M0"], k2["M1"], length(gene_order),
  paste(sprintf("%s, %s", color_name(CLUSTER_COLORS[names(CLUSTERS)]), tolower(names(CLUSTERS))), collapse = "; "))
writeLines(legend, file.path(pe(OUT_DIR), "Figure5_legend_draft.txt"))

cat("\n=== Resumo ===\n")
print(stats_tn |> transmute(cohort, gene, n_pairs, log2FC = round(median_log2FC_paired, 2),
                            r = round(rank_biserial_paired, 2), p_adj = signif(p_adj_paired, 2)), n = 40)
print(stats_m |> transmute(cohort, gene, M0 = n_M0, M1 = n_M1, log2FC = round(log2FC_median, 2),
                           delta = round(cliffs_delta, 2), p_adj = signif(p_adj, 2)), n = 40)
cat("\nSaídas em", pe(OUT_DIR), "\n")
