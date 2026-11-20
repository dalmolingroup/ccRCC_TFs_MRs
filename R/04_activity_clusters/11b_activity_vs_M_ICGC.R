### Atividade dos clusters vs metástase ao diagnóstico no ICGC-RECA (e lado a lado com o TCGA) ###
# Revisão CSBJ-D-26-01039 - Revisor 2, Comentário 3, bullet 2:
#   "Report the associations between the two cluster-activity scores and metastatic status directly
#    in ICGC-RECA, with effect sizes and confidence intervals."
#
# Mesmo escore e mesmas métricas da Fig. 4C (TCGA):
#   - escore do cluster = média dos dES dos TFs do cluster
#   - Cliff's delta (M1 vs M0) + IC95% por bootstrap, Mann-Whitney
#   - OR de metástase por DP (regressão logística; escore padronizado dentro da coorte), IC95% de
#     Wald; com poucos M1, também a regressão de Firth (logistf), se instalada
#   - sensibilidade exploratória: Cluster 1 sem PARP12 (não reportada no artigo)
#
# ENTRADAS
#   ACTIVITY_CSV: atividade do ICGC SÓ com os 91 tumores (RT), regulons do TCGA + tni.replace.samples +
#                 tni.gsea2 - a mesma fonte usada no Cox (21_survival_ICGC.R) e na Fig. S5
#   ICGC_RDATA:   donor_data_final (tissue_type, submitted_donor_id, pathologic_m de todas as amostras)
#   TCGA_FILE:    (opcional) dados_por_paciente.csv do script 11 (TCGA) para a tabela lado a lado
#
# SAÍDAS (OUT_DIR)
#   clusters_M1_vs_M0_TCGA_ICGC.csv  - Supplementary Table S28 (OR de Firth no texto do ICGC)
#   escores_clusters_M1_vs_M0_ICGC.(pdf|png|tiff) - versão isolada do painel D da Supplementary Figure S5
#                                                   (a figura final é montada pelo 12_FigS5_ICGC_activity.R)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

set.seed(2026)

## ---------------------------------------------------------------- CONFIG
ACTIVITY_CSV <- "~/Kidney/0_data/dES_ICGC_RT91.csv"   # 91 tumores RT (gerado pelo 21_survival_ICGC.R)
ICGC_RDATA   <- "~/Kidney/0_data/ICGC_DEA_Data.RData"   # donor_data_final: metadados de TODAS as amostras
TCGA_FILE    <- "~/Kidney/results_activity_clinical/dados_por_paciente.csv"   # NULL para pular o TCGA
TUMOR_COL <- "tissue_type"; TUMOR_VALUE <- "RT"
DONOR_COL <- "submitted_donor_id"
M_COL     <- "pathologic_m"
OUT_DIR   <- "results_activity_clinical_ICGC"
N_BOOT    <- 2000

TF14 <- data.frame(
  symbol  = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A",
              "IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  ensembl = c("ENSG00000111206", "ENSG00000101057", "ENSG00000050344", "ENSG00000059378", "ENSG00000137203",
              "ENSG00000117595", "ENSG00000116604", "ENSG00000151623", "ENSG00000185630", "ENSG00000077092",
              "ENSG00000143365", "ENSG00000182568", "ENSG00000162367", "ENSG00000196628"),
  cluster = rep(c("Cluster 1", "Cluster 2"), c(5, 9)))

# escores analisados (o 3º é a sensibilidade sem PARP12)
SCORES <- list(
  "Cluster 1"                = TF14$symbol[TF14$cluster == "Cluster 1"],
  "Cluster 2"                = TF14$symbol[TF14$cluster == "Cluster 2"],
  "Cluster 1 (excl. PARP12)" = setdiff(TF14$symbol[TF14$cluster == "Cluster 1"], "PARP12"))
CL_COL <- c("Cluster 1" = "#542788", "Cluster 2" = "#1a9850")
INK2   <- "#52514e"

pe <- path.expand
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- 1. atividade do ICGC (45 RT)
dES <- as.matrix(read.csv(pe(ACTIVITY_CSV), row.names = 1, check.names = FALSE))
# colunas: aceita símbolo ou ENSG (com ou sem versão)
key <- sub("\\..*", "", colnames(dES))
colnames(dES) <- ifelse(key %in% TF14$ensembl, TF14$symbol[match(key, TF14$ensembl)], colnames(dES))
miss <- setdiff(TF14$symbol, colnames(dES))
if (length(miss)) stop("TFs sem atividade no ICGC: ", paste(miss, collapse = ", "))

# metadados de TODAS as amostras do ICGC (o metadata_ICGC_full.csv só tem as 90 amostras pareadas)
load(pe(ICGC_RDATA))
meta <- as.data.frame(donor_data_final)
rownames(meta) <- paste0(meta$submitted_donor_id, meta$tissue_type)
s <- intersect(rownames(dES), rownames(meta))
cat(sprintf("Atividade: %d amostras | com metadados: %d\n", nrow(dES), length(s)))   # esperado 91 | 91
meta <- meta[s, ]

# trava: só tumores, 1 por doador (o erro corrigido na revisão)
stopifnot(all(meta[[TUMOR_COL]] == TUMOR_VALUE), !anyDuplicated(meta[[DONOR_COL]]))

icgc <- data.frame(sample = s, M = as.character(meta[[M_COL]]))
for (sc in names(SCORES)) icgc[[sc]] <- rowMeans(dES[s, SCORES[[sc]], drop = FALSE])
icgc <- filter(icgc, M %in% c("M0", "M1"))
cat(sprintf("ICGC-RECA: %d tumores com M conhecido (M0 = %d, M1 = %d)\n",
            nrow(icgc), sum(icgc$M == "M0"), sum(icgc$M == "M1")))

## ---------------------------------------------------------------- 2. TCGA (opcional, mesmo escore)
cohorts <- list(`ICGC-RECA` = icgc)
if (!is.null(TCGA_FILE) && file.exists(TCGA_FILE)) {
  tc <- read.csv(TCGA_FILE, check.names = FALSE)
  # no script 11 os escores se chamam cluster_1 / cluster_2 e os TFs estão nas colunas (símbolo ou ENSG)
  key <- sub("\\..*", "", names(tc))
  names(tc) <- ifelse(key %in% TF14$ensembl, TF14$symbol[match(key, TF14$ensembl)], names(tc))
  tc$M <- ifelse(tc$M %in% c("M0", "0"), "M0", ifelse(tc$M %in% c("M1", "1"), "M1", NA))
  for (sc in names(SCORES)) tc[[sc]] <- rowMeans(tc[, SCORES[[sc]], drop = FALSE])
  cohorts <- c(list(`TCGA-KIRC` = filter(tc, M %in% c("M0", "M1"))), cohorts)
  cat(sprintf("TCGA-KIRC: %d tumores (M0 = %d, M1 = %d)\n", nrow(cohorts$`TCGA-KIRC`),
              sum(cohorts$`TCGA-KIRC`$M == "M0"), sum(cohorts$`TCGA-KIRC`$M == "M1")))
}

## ---------------------------------------------------------------- 3. métricas
cliff <- function(x, y) (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
boot_cliff <- function(x, y) {
  b <- replicate(N_BOOT, cliff(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  quantile(b, c(0.025, 0.975), names = FALSE)
}
has_firth <- requireNamespace("logistf", quietly = TRUE)

assoc <- function(d, sc, cohort) {
  x <- d[[sc]][d$M == "M1"]; y <- d[[sc]][d$M == "M0"]
  dd <- data.frame(m1 = as.integer(d$M == "M1"), z = as.numeric(scale(d[[sc]])))  # OR por DP da coorte
  fit <- glm(m1 ~ z, family = binomial, data = dd)
  ci  <- confint.default(fit)["z", ]                                               # IC de Wald, como no manuscrito
  cd_ci <- boot_cliff(x, y)
  out <- data.frame(
    cohort = cohort, score = sc, n_M0 = length(y), n_M1 = length(x),
    median_M0 = median(y), median_M1 = median(x),
    cliffs_delta = cliff(x, y), delta_CI_low = cd_ci[1], delta_CI_high = cd_ci[2],
    p_mannwhitney = wilcox.test(x, y, exact = FALSE)$p.value,
    OR_per_SD = exp(coef(fit)["z"]), OR_CI_low = exp(ci[1]), OR_CI_high = exp(ci[2]),
    p_logistic = summary(fit)$coefficients["z", 4],
    OR_firth = NA_real_, OR_firth_CI_low = NA_real_, OR_firth_CI_high = NA_real_, p_firth = NA_real_)
  if (has_firth) {
    ff <- logistf::logistf(m1 ~ z, data = dd)
    j  <- which(names(coef(ff)) == "z")
    out$OR_firth <- exp(coef(ff)[j]); out$OR_firth_CI_low <- exp(ff$ci.lower[j])
    out$OR_firth_CI_high <- exp(ff$ci.upper[j]); out$p_firth <- ff$prob[j]
  }
  out
}

res <- bind_rows(lapply(names(cohorts), function(co)
  bind_rows(lapply(names(SCORES), function(sc) assoc(cohorts[[co]], sc, co)))))
rownames(res) <- NULL
write.csv(res, file.path(OUT_DIR, "clusters_M1_vs_M0_TCGA_ICGC.csv"), row.names = FALSE)

cat("\n=== Escores dos clusters vs metástase (M1 vs M0) ===\n")
print(res |> transmute(cohort, score, `M0/M1` = paste0(n_M0, "/", n_M1),
                       delta = sprintf("%.2f (%.2f to %.2f)", cliffs_delta, delta_CI_low, delta_CI_high),
                       p_MW = signif(p_mannwhitney, 2),
                       OR_SD = sprintf("%.2f (%.2f-%.2f)", OR_per_SD, OR_CI_low, OR_CI_high),
                       OR_Firth = ifelse(is.na(OR_firth), "-",
                                         sprintf("%.2f (%.2f-%.2f)", OR_firth, OR_firth_CI_low, OR_firth_CI_high))),
      row.names = FALSE)
if (!has_firth) cat("\n[aviso] pacote logistf ausente: com poucos M1, instale para a estimativa de Firth\n")

## ---------------------------------------------------------------- 4. figura (análoga à Fig. 4C, ICGC)
d <- cohorts$`ICGC-RECA`
long <- bind_rows(lapply(c("Cluster 1", "Cluster 2"), function(sc)
  data.frame(cluster = sc, M = d$M, score = d[[sc]])))
ann <- res |> filter(cohort == "ICGC-RECA", score %in% c("Cluster 1", "Cluster 2")) |>
  transmute(cluster = score,
            txt = sprintf("δ = %.2f (%.2f to %.2f)\nOR/SD = %.2f (%.2f-%.2f)",
                          cliffs_delta, delta_CI_low, delta_CI_high, OR_per_SD, OR_CI_low, OR_CI_high))

g <- ggplot(long, aes(M, score, fill = cluster)) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.45, color = "#8a8983") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, color = "#0b0b0b", linewidth = 0.35) +
  geom_text(data = ann, aes(x = 1.5, y = Inf, label = txt), vjust = 1.2, size = 2.6, color = INK2,
            inherit.aes = FALSE) +
  scale_fill_manual(values = CL_COL, guide = "none") +
  scale_x_discrete(labels = c(M0 = sprintf("M0\n(n = %d)", sum(d$M == "M0")),
                              M1 = sprintf("M1\n(n = %d)", sum(d$M == "M1")))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  facet_wrap(~cluster, scales = "free_y") +
  labs(x = NULL, y = "Cluster activity score (mean dES)",
       title = sprintf("ICGC-RECA (n = %d primary tumours)", nrow(d))) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(color = INK2), strip.text = element_text(face = "bold"),
        plot.title = element_text(size = 9, hjust = 0.5))

for (ext in c("pdf", "png", "tiff")) {
  f <- file.path(OUT_DIR, paste0("escores_clusters_M1_vs_M0_ICGC.", ext))
  if (ext == "pdf") ggsave(f, g, width = 5, height = 3.4, device = cairo_pdf)
  else if (ext == "png") ggsave(f, g, width = 5, height = 3.4, dpi = 300, bg = "white")
  else ggsave(f, g, width = 5, height = 3.4, dpi = 300, compression = "lzw", bg = "white")
}
cat("\nSaídas em", OUT_DIR, "\n")
