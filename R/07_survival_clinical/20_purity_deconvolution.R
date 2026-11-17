### Pureza tumoral e composição imune vs atividade dos clusters de regulons ###
# Revisão CSBJ-D-26-01039
#   R3.7  : o sinal imune pode refletir infiltração/pureza e não regulação nas células tumorais
#   R2.M6 : nome "immune-related cluster"
#   R2.a14: pureza como limitação
#
# Reaproveita a deconvolução já feita para o MOFA (IOBR, TCGA-KIRC STAR counts -> TPM):
#   resultados/checkpoints/10_tme_combine.rds  (ou os .rds de cada método)
#
# Perguntas:
#   1. Quanto o escore de cada cluster se correlaciona com pureza, ImmuneScore/StromalScore e
#      com cada tipo celular (quanTIseq como principal, CIBERSORT como sensibilidade)?
#   2. A diferença M1 vs M0 se mantém ajustando pela pureza / ImmuneScore? (regressão logística)
#   3. O HR do cluster no Cox se mantém ajustando pela pureza / ImmuneScore?
#
# A atividade é lida exatamente como em survival_cox_cdr.R (mesmo objeto, mesmo escore).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(survival); library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG
TME_FILE     <- "~/Multiomic/Novos Códigos/results_mofa/checkpoints/10_tme_combine.rds"
CPE_FILE     <- NULL     # opcional: pureza de consenso (Aran et al., Nat Commun 2015, Supp. Data 1),
# CSV com colunas Sample.ID e CPE; NULL = só ESTIMATE
RTN_FILE     <- "0_data/rtni_TF.RData"
RTN_OBJECT   <- "rtni_TF"
ACTIVITY_RDS <- NULL
CDR_FILE     <- "0_data/TCGA-CDR-SupplementalTableS1.xlsx"
META_FILE    <- "0_data/metadata_KIRC_full.csv"
META_ID_COL  <- "sample"
M_COL        <- "ajcc_pathologic_m"
SYMBOL_FILE  <- "0_data/ensembl_to_symbol.csv"
OUT_DIR      <- "results_purity_deconv"

# Clusters por SÍMBOLO (convertidos para o formato dos regulons do RTN automaticamente)
CLUSTERS <- list(
  differentiation_associated = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  proliferative = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A")
)

ENDPOINTS <- c("OS", "PFI", "DSS")
COVARS    <- c("age", "sex", "stage", "grade")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
strip <- function(x) sub("\\..*", "", x)

## ---------------------------------------------------------------- 1. atividade e escores dos clusters
symbols <- read.csv(SYMBOL_FILE)
sym2ens <- setNames(strip(symbols$ensembl), toupper(symbols$symbol))
ens2sym <- setNames(symbols$symbol, strip(symbols$ensembl))

get_activity <- function(tfs) {
  if (!is.null(ACTIVITY_RDS)) return(readRDS(ACTIVITY_RDS))
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(RTN_FILE, envir = env)
  rtni <- get(RTN_OBJECT, envir = env)
  act <- tryCatch(tni.get(rtni, what = "regulonActivity"), error = function(e) NULL)
  have <- if (!is.null(act$differential)) colnames(act$differential) else character(0)
  if (!all(tfs %in% have)) {
    message("Rodando tni.gsea2 para ", length(tfs), " regulons...")
    rtni <- tni.gsea2(rtni, regulatoryElements = tfs)
    act <- tni.get(rtni, what = "regulonActivity")
  }
  act$differential
}

# descobre se os regulons do RTN estão em ENSG (com ou sem versão) ou em símbolo
probe_names <- function() {
  if (!is.null(ACTIVITY_RDS)) return(colnames(readRDS(ACTIVITY_RDS)))
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(RTN_FILE, envir = env)
  tni.get(get(RTN_OBJECT, envir = env), what = "regulatoryElements")
}

reg_names <- probe_names()
to_reg <- function(sym) {
  nm <- names(reg_names)                                                    # RTN: nomes = símbolos, valores = IDs
  if (!is.null(nm) && all(toupper(sym) %in% toupper(nm))) return(unname(reg_names[match(toupper(sym), toupper(nm))]))
  reg_names <- unname(reg_names)
  hit <- reg_names[toupper(reg_names) %in% toupper(sym)]                    # símbolo
  if (length(hit) == length(sym)) return(hit[match(toupper(sym), toupper(hit))])
  ens <- sym2ens[toupper(sym)]
  out <- reg_names[match(ens, strip(reg_names))]                           # ENSG
  if (any(is.na(out))) stop("Regulons não encontrados no RTN: ", paste(sym[is.na(out)], collapse = ", "))
  out
}
CLUSTER_REG <- lapply(CLUSTERS, to_reg)
all_reg <- unlist(CLUSTER_REG)

dES <- get_activity(all_reg)
act <- as.data.frame(dES[, unlist(CLUSTERS)  , drop = FALSE])
colnames(act) <- unlist(CLUSTERS)                                          # colunas por símbolo
act$sample <- rownames(dES)
act <- act |> filter(substr(sample, 14, 15) == "01") |>
  mutate(patient = substr(sample, 1, 12)) |> distinct(patient, .keep_all = TRUE)
for (cl in names(CLUSTERS)) act[[paste0(cl, "_cluster")]] <- rowMeans(act[, CLUSTERS[[cl]], drop = FALSE])
SCORES <- paste0(names(CLUSTERS), "_cluster")   # escore = média dos dES do cluster

## ---------------------------------------------------------------- 2. deconvolução (IOBR)
tme <- readRDS(TME_FILE)
tme <- tme |>
  mutate(sample_type = substr(ID, 14, 15), patient = substr(ID, 1, 12)) |>
  filter(sample_type == "01") |>
  arrange(patient, ID) |> distinct(patient, .keep_all = TRUE)             # 1 tumor primário por paciente

est_cols <- grep("_estimate$", names(tme), value = TRUE)
qts_cols <- setdiff(grep("_quantiseq$", names(tme), value = TRUE), "Other_quantiseq")
cib_cols <- grep("_CIBERSORT$", names(tme), value = TRUE)
cib_cols <- cib_cols[!grepl("^(P-value|Correlation|RMSE)_", cib_cols)]
message(sprintf("ESTIMATE: %d colunas; quanTIseq: %d; CIBERSORT: %d", length(est_cols), length(qts_cols), length(cib_cols)))

purity_col <- grep("^TumorPurity", est_cols, value = TRUE)[1]
immune_col <- grep("^ImmuneScore", est_cols, value = TRUE)[1]
stroma_col <- grep("^StromalScore", est_cols, value = TRUE)[1]
tme$purity <- tme[[purity_col]]; tme$immune_score <- tme[[immune_col]]; tme$stromal_score <- tme[[stroma_col]]
PURITY_VARS <- c(purity = "ESTIMATE tumour purity", immune_score = "ESTIMATE ImmuneScore",
                 stromal_score = "ESTIMATE StromalScore")

if (!is.null(CPE_FILE) && file.exists(CPE_FILE)) {
  cpe <- read.csv(CPE_FILE, check.names = FALSE)
  cpe <- data.frame(patient = substr(cpe$Sample.ID, 1, 12), cpe = suppressWarnings(as.numeric(cpe$CPE))) |>
    filter(!is.na(cpe)) |> distinct(patient, .keep_all = TRUE)
  tme <- left_join(tme, cpe, by = "patient")
  PURITY_VARS <- c(cpe = "Consensus purity (CPE)", PURITY_VARS)
}

## ---------------------------------------------------------------- 3. clínico (CDR) e M
cdr <- read_excel(CDR_FILE, sheet = "TCGA-CDR", na = c("", "NA", "#N/A", "[Not Available]", "[Not Applicable]",
                                                       "[Unknown]", "[Discrepancy]", "[Not Evaluated]"))
clin <- cdr |> filter(type == "KIRC") |>
  transmute(patient = bcr_patient_barcode,
            age = as.numeric(age_at_initial_pathologic_diagnosis),
            sex = factor(tolower(gender), levels = c("female", "male")),
            stage = case_when(grepl("Stage IV", ajcc_pathologic_tumor_stage) ~ "IV",
                              grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
                              grepl("Stage I", ajcc_pathologic_tumor_stage) ~ "I-II"),
            grade = case_when(histological_grade %in% c("G1", "G2") ~ "G1-G2",
                              histological_grade == "G3" ~ "G3", histological_grade == "G4" ~ "G4"),
            across(all_of(as.vector(rbind(ENDPOINTS, paste0(ENDPOINTS, ".time")))), as.numeric)) |>
  mutate(stage = factor(stage, levels = c("I-II", "III", "IV")), grade = factor(grade, levels = c("G1-G2", "G3", "G4")))
meta <- read.csv(META_FILE, check.names = FALSE)
meta <- data.frame(patient = substr(meta[[META_ID_COL]], 1, 12), M = meta[[M_COL]]) |> distinct(patient, .keep_all = TRUE)

keep_tme <- c("patient", names(PURITY_VARS), qts_cols, cib_cols)
df <- act |> inner_join(tme[, keep_tme], by = "patient") |> left_join(clin, by = "patient") |> left_join(meta, by = "patient")
df_ml <- df |> filter(M %in% c("M0", "M1"))                                # coorte do ML (M0/M1 definidos)
flow <- c(activity = nrow(act), with_deconvolution = nrow(df), M0 = sum(df$M == "M0", na.rm = TRUE),
          M1 = sum(df$M == "M1", na.rm = TRUE))
print(flow)
write.csv(data.frame(step = names(flow), n = flow), file.path(OUT_DIR, "fluxo.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 4. correlações
lab_cell <- function(x) gsub("_", " ", sub("_(quantiseq|CIBERSORT|estimate)$", "", x))
cor_block <- function(vars, cells, method_name) {
  out <- expand.grid(variable = vars, cell = cells, stringsAsFactors = FALSE)
  res <- t(mapply(function(v, c) {
    ok <- complete.cases(df_ml[[v]], df_ml[[c]])
    if (sd(df_ml[[c]][ok]) == 0) return(c(NA, NA))
    ct <- suppressWarnings(cor.test(df_ml[[v]][ok], df_ml[[c]][ok], method = "spearman"))
    c(unname(ct$estimate), ct$p.value)
  }, out$variable, out$cell))
  out$rho <- res[, 1]; out$p <- res[, 2]; out$method <- method_name
  out
}
vars_all <- c(SCORES, unlist(CLUSTERS))
cors <- bind_rows(
  cor_block(vars_all, names(PURITY_VARS), "ESTIMATE/CPE"),
  cor_block(vars_all, qts_cols, "quanTIseq"),
  cor_block(vars_all, cib_cols, "CIBERSORT")
) |> group_by(method) |> mutate(FDR = p.adjust(p, "BH")) |> ungroup() |>
  mutate(cell_label = ifelse(cell %in% names(PURITY_VARS), PURITY_VARS[cell], lab_cell(cell)))
write.csv(cors, file.path(OUT_DIR, "correlacoes_atividade_vs_TME.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 5. M1 vs M0 ajustado (logística)
z <- function(x) as.numeric(scale(x))
logit_tab <- list()
adj_sets <- c(list(none = character(0)), as.list(setNames(names(PURITY_VARS), names(PURITY_VARS))),
              list(purity_immune = intersect(c("purity", "immune_score"), names(PURITY_VARS))))
for (s in SCORES) for (a in names(adj_sets)) {
  d <- df_ml[complete.cases(df_ml[, c(s, adj_sets[[a]])]), ]
  d$y <- as.integer(d$M == "M1"); d$x <- z(d[[s]])
  for (v in adj_sets[[a]]) d[[v]] <- z(d[[v]])
  fit <- glm(as.formula(paste("y ~", paste(c("x", adj_sets[[a]]), collapse = " + "))), data = d, family = binomial)
  ci <- suppressMessages(confint.default(fit)["x", ])
  logit_tab[[length(logit_tab) + 1]] <- data.frame(score = s, adjusted_for = a, n = nrow(d), M1 = sum(d$y),
                                                   OR_per_SD = exp(coef(fit)["x"]), CI_low = exp(ci[1]), CI_high = exp(ci[2]),
                                                   p = summary(fit)$coefficients["x", 4])
}
logit_tab <- bind_rows(logit_tab)
write.csv(logit_tab, file.path(OUT_DIR, "M1_vs_M0_logistica_ajustada.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 6. Cox ajustado por pureza
cohorts <- list(all = df_ml, M0_only = filter(df_ml, M == "M0"))
cox_tab <- list()
for (coh in names(cohorts)) for (e in ENDPOINTS) for (s in SCORES) for (a in names(adj_sets)) {
  tcol <- paste0(e, ".time"); extra <- adj_sets[[a]]
  d <- cohorts[[coh]]
  d <- d[complete.cases(d[, c(tcol, e, s, COVARS, extra)]) & d[[tcol]] > 0, ]
  d$stage <- droplevels(d$stage); d$x <- z(d[[s]])
  for (v in extra) d[[v]] <- z(d[[v]])
  fit <- coxph(as.formula(paste0("Surv(", tcol, ", ", e, ") ~ ", paste(c("x", COVARS, extra), collapse = " + "))), data = d)
  sm <- summary(fit)
  cox_tab[[length(cox_tab) + 1]] <- data.frame(cohort = coh, endpoint = e, score = s, adjusted_for = a,
                                               n = sm$n, events = sm$nevent, HR = sm$conf.int["x", 1], CI_low = sm$conf.int["x", 3],
                                               CI_high = sm$conf.int["x", 4], p = sm$coefficients["x", 5])
}
cox_tab <- bind_rows(cox_tab)
write.csv(cox_tab, file.path(OUT_DIR, "cox_ajustado_pureza.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 7. figuras
theme_set(theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), axis.text = element_text(color = "#52514e")))
# (a) heatmap: escores e TFs x pureza e células (quanTIseq)
hm <- cors |> filter(method %in% c("ESTIMATE/CPE", "quanTIseq")) |>
  mutate(variable = factor(variable, levels = rev(vars_all)),
         cell_label = factor(cell_label, levels = unique(cell_label)),
         star = ifelse(!is.na(FDR) & FDR < 0.05, "*", ""))
g1 <- ggplot(hm, aes(cell_label, variable, fill = rho)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f%s", rho, star)), size = 2.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "Spearman ρ") +
  labs(x = NULL, y = NULL, caption = "* BH-adjusted p < 0.05") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "heatmap_atividade_vs_TME.pdf"), g1, width = 8, height = 5.5)
ggsave(file.path(OUT_DIR, "heatmap_atividade_vs_TME.png"), g1, width = 8, height = 5.5, dpi = 300)

# (b) forest: HR do escore com e sem ajuste por pureza
fp <- cox_tab |> mutate(cohort = recode(cohort, all = "All patients", M0_only = "M0 only"),
                        adjusted_for = factor(adjusted_for, levels = rev(names(adj_sets))))
g2 <- ggplot(fp, aes(HR, adjusted_for, color = score)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, position = position_dodge(0.5)) +
  geom_point(position = position_dodge(0.5), size = 1.8) + scale_x_log10() +
  scale_color_manual(values = c(differentiation_associated_cluster = "#1a9850", proliferative_cluster = "#542788"),
                     labels = c("Differentiation-associated cluster", "Proliferative cluster"), name = NULL) +
  facet_grid(cohort ~ endpoint) +
  labs(x = "HR per SD (log scale), adjusted for age, sex, stage, grade + listed variable", y = "Additional adjustment") +
  theme(legend.position = "bottom", panel.border = element_rect(fill = NA, color = "grey85"))
ggsave(file.path(OUT_DIR, "forest_cox_ajuste_pureza.pdf"), g2, width = 8, height = 4.5)
ggsave(file.path(OUT_DIR, "forest_cox_ajuste_pureza.png"), g2, width = 8, height = 4.5, dpi = 300)

## ---------------------------------------------------------------- 8. resumo
cat("\n=== Escores x pureza/ImmuneScore (Spearman) ===\n")
print(cors |> filter(variable %in% SCORES, method == "ESTIMATE/CPE") |>
        transmute(variable, cell_label, rho = round(rho, 2), FDR = signif(FDR, 2)))
cat("\n=== Células mais correlacionadas com o escore imune (quanTIseq) ===\n")
print(cors |> filter(variable == "differentiation_associated_cluster", method == "quanTIseq") |> arrange(desc(abs(rho))) |>
        transmute(cell_label, rho = round(rho, 2), FDR = signif(FDR, 2)))
cat("\n=== M1 vs M0: OR por DP, sem e com ajuste ===\n")
print(logit_tab |> mutate(across(c(OR_per_SD, CI_low, CI_high), ~ round(.x, 2)), p = signif(p, 2)))
cat("\n=== Cox: escore imune, sem e com ajuste por pureza ===\n")
print(cox_tab |> filter(score == "differentiation_associated_cluster") |>
        transmute(cohort, endpoint, adjusted_for, n, events, HR = round(HR, 2),
                  CI = sprintf("%.2f-%.2f", CI_low, CI_high), p = signif(p, 3)))
cat("\nSaídas em", OUT_DIR, "\n")
