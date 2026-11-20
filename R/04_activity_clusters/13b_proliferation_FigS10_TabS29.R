### Comentário 4 (R2): os dois programas são independentes do eixo de proliferação? ###
# Revisão CSBJ-D-26-01039 - Revisor 2, Comentário 4. Responde aos 3 bullets analíticos:
#   (a) "Comparing both cluster scores with established proliferation-related expression scores,
#        such as Hallmark E2F targets or G2M checkpoint activity"
#   (b) "Testing whether the differentiation-associated score retains independent associations with
#        metastasis or survival after adjustment for a conventional proliferation score"
#   (c) "Repeating the principal cluster-level analyses after excluding the [...] regulators whose
#        cluster assignments [...] were inconsistent across cohorts"
#        -> após a correção do Comentário 1, nenhum TF é discordante na especificação principal;
#           o único discordante na sensibilidade (Ward) é o RORC. Reporta-se: sem RORC (principal)
#           e sem os 4 TFs citados pelo revisor (TFAP2A, MEF2D, TAL1, TCF4), por completude.
#
# Escores de proliferação: média dos z-scores (por gene, entre os tumores) dos genes dos conjuntos
#   HALLMARK_E2F_TARGETS e HALLMARK_G2M_CHECKPOINT (MSigDB, msigdbr) na expressão voom do TCGA-KIRC.
#   "Proliferation score" = média dos dois escores (padronizados).
#
# ENTRADAS
#   results_survival/dados_sobrevivencia.csv  (gerado pelo 18_survival_cox_TCGA.R: atividade por TF,
#                                              escores dos clusters, CDR, M)
#   0_data/rtni_TF.RData                      (expressão voom dos tumores, a mesma da rede)
# SAÍDAS (OUT_DIR): Table_prolif_*.csv, TableS_comment4.xlsx (= Supplementary Table S29),
#                   FigS_prolif_correlation.(pdf|png) (= Supplementary Figure S10)

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(ggplot2); library(openxlsx)
})

## ---------------------------------------------------------------- CONFIG
SURV_FILE <- "results_survival/dados_sobrevivencia.csv"
RTN_FILE  <- "0_data/rtni_TF.RData"; RTN_OBJECT <- "rtni_TF"
OUT_DIR   <- "results_comment4_proliferation"
ENDPOINTS <- c("OS", "PFI", "DSS")
COVARS    <- c("age", "sex", "stage", "grade")          # as mesmas do 18_survival_cox_TCGA.R
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TF14 <- data.frame(
  symbol  = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A",
              "IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  ensembl = c("ENSG00000111206", "ENSG00000101057", "ENSG00000050344", "ENSG00000059378", "ENSG00000137203",
              "ENSG00000117595", "ENSG00000116604", "ENSG00000151623", "ENSG00000185630", "ENSG00000077092",
              "ENSG00000143365", "ENSG00000182568", "ENSG00000162367", "ENSG00000196628"),
  cluster = rep(c("C1", "C2"), c(5, 9)))
C1 <- TF14$symbol[TF14$cluster == "C1"]; C2 <- TF14$symbol[TF14$cluster == "C2"]
VARIANTS <- list(
  "All 14 regulators"                  = list(C1 = C1, C2 = C2),
  "Excluding RORC"                     = list(C1 = C1, C2 = setdiff(C2, "RORC")),
  "Excluding TFAP2A, MEF2D, TAL1, TCF4" = list(C1 = setdiff(C1, "TFAP2A"), C2 = setdiff(C2, c("MEF2D", "TAL1", "TCF4"))))
LAB <- c(C1 = "Proliferative cluster", C2 = "Differentiation-associated cluster")

strip_v <- function(x) sub("\\..*", "", x)
zs <- function(x) as.numeric(scale(x))

## ---------------------------------------------------------------- 1. dados por paciente
df <- read.csv(SURV_FILE, check.names = FALSE)
key <- strip_v(names(df))
names(df) <- ifelse(key %in% TF14$ensembl, TF14$symbol[match(key, TF14$ensembl)], names(df))
miss <- setdiff(TF14$symbol, names(df))
if (length(miss)) stop("TFs ausentes em ", SURV_FILE, ": ", paste(miss, collapse = ", "))
df <- df |> mutate(sex = factor(sex, levels = c("female", "male")),
                   stage = factor(stage, levels = c("I-II", "III", "IV")),
                   grade = factor(grade, levels = c("G1-G2", "G3", "G4")))
cat(sprintf("TCGA-KIRC: %d pacientes (M0 = %d, M1 = %d)\n", nrow(df),
            sum(df$M == "M0", na.rm = TRUE), sum(df$M == "M1", na.rm = TRUE)))

## ---------------------------------------------------------------- 2. escores Hallmark (E2F, G2M)
suppressPackageStartupMessages({ library(RTN); library(msigdbr) })
env <- new.env(); load(RTN_FILE, envir = env)
gexp <- tni.get(get(RTN_OBJECT, envir = env), "gexp")            # genes (ENSG) x amostras
rownames(gexp) <- strip_v(rownames(gexp))

hm <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
ens_col <- intersect(c("ensembl_gene", "db_ensembl_gene"), names(hm))[1]
sets <- list(E2F = unique(hm[[ens_col]][hm$gs_name == "HALLMARK_E2F_TARGETS"]),
             G2M = unique(hm[[ens_col]][hm$gs_name == "HALLMARK_G2M_CHECKPOINT"]))

# só os tumores do df (1 por paciente)
smp <- colnames(gexp)[substr(colnames(gexp), 1, 12) %in% df$patient & substr(colnames(gexp), 14, 15) == "01"]
smp <- smp[!duplicated(substr(smp, 1, 12))]
Z <- t(scale(t(gexp[, smp])))                                   # z-score por gene entre os tumores
hall <- data.frame(patient = substr(smp, 1, 12))
for (s in names(sets)) {
  g <- intersect(sets[[s]], rownames(Z))
  cat(sprintf("%s: %d de %d genes encontrados\n", s, length(g), length(sets[[s]])))
  hall[[s]] <- colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
}
hall$Prolif <- rowMeans(cbind(zs(hall$E2F), zs(hall$G2M)))
df <- inner_join(df, hall, by = "patient")

# escores dos clusters em cada variante
for (v in names(VARIANTS)) for (cl in c("C1", "C2"))
  df[[paste(cl, v, sep = "|")]] <- rowMeans(df[, VARIANTS[[v]][[cl]], drop = FALSE])
s_main <- function(cl) paste(cl, "All 14 regulators", sep = "|")

## ---------------------------------------------------------------- 3. (a) correlações com proliferação
sp <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman")); c(rho = unname(ct$estimate), p = ct$p.value) }
pcor_s <- function(x, y, z) {                                   # Spearman parcial (resíduos dos postos)
  rx <- resid(lm(rank(x) ~ rank(z))); ry <- resid(lm(rank(y) ~ rank(z))); sp(rx, ry)
}
tab_cor <- bind_rows(
  lapply(c("C1", "C2"), function(cl) bind_rows(lapply(c("E2F", "G2M", "Prolif"), function(p) {
    r <- sp(df[[s_main(cl)]], df[[p]])
    data.frame(comparison = paste(LAB[cl], "vs", p), rho = r["rho"], p = r["p"])
  }))),
  { r <- sp(df[[s_main("C1")]], df[[s_main("C2")]])
    data.frame(comparison = "Proliferative vs differentiation-associated cluster", rho = r["rho"], p = r["p"]) },
  { r <- pcor_s(df[[s_main("C1")]], df[[s_main("C2")]], df$Prolif)
    data.frame(comparison = "Proliferative vs differentiation-associated cluster | proliferation score (partial)",
               rho = r["rho"], p = r["p"]) })
rownames(tab_cor) <- NULL
write.csv(tab_cor, file.path(OUT_DIR, "Table_prolif_correlations.csv"), row.names = FALSE)
cat("\n=== (a) Correlações de Spearman (n =", nrow(df), ") ===\n"); print(tab_cor |> mutate(rho = round(rho, 2), p = signif(p, 2)))

## ---------------------------------------------------------------- 4. (b) metástase ajustada por proliferação
dm <- df |> filter(M %in% c("M0", "M1")) |> mutate(m1 = as.integer(M == "M1"))
logit_or <- function(d, terms) {
  fit <- glm(as.formula(paste("m1 ~", paste(terms, collapse = " + "))), family = binomial, data = d)
  ci <- confint.default(fit)                                    # Wald, como no manuscrito
  list(fit = fit, or = exp(coef(fit)[terms[1]]), lo = exp(ci[terms[1], 1]), hi = exp(ci[terms[1], 2]),
       p = summary(fit)$coefficients[terms[1], 4])
}
tab_met <- bind_rows(lapply(c("C1", "C2"), function(cl) {
  d <- dm |> mutate(score = zs(.data[[s_main(cl)]]), prolif = zs(Prolif))
  m0 <- glm(m1 ~ prolif, family = binomial, data = d)
  bind_rows(lapply(list("Cluster score" = "score", "+ proliferation score" = c("score", "prolif")), function(t) {
    r <- logit_or(d, t)
    lrt <- if (length(t) > 1) anova(m0, r$fit, test = "LRT")$`Pr(>Chi)`[2] else NA
    data.frame(cluster = LAB[cl], model = if (length(t) > 1) "Adjusted for proliferation" else "Unadjusted",
               n = nrow(d), M1 = sum(d$m1), OR_per_SD = r$or, CI_low = r$lo, CI_high = r$hi, p = r$p,
               LRT_p_added_to_proliferation = lrt)
  }))
}))
write.csv(tab_met, file.path(OUT_DIR, "Table_prolif_metastasis.csv"), row.names = FALSE)
cat("\n=== (b) Metástase (M1 vs M0), OR por DP ===\n")
print(tab_met |> transmute(cluster, model, OR = sprintf("%.2f (%.2f-%.2f)", OR_per_SD, CI_low, CI_high),
                           p = signif(p, 2), LRT_p = signif(LRT_p_added_to_proliferation, 2)), row.names = FALSE)

## ---------------------------------------------------------------- 5. (b) sobrevida ajustada por proliferação
cox_hr <- function(d, e, terms) {
  tc <- paste0(e, ".time")
  d <- d[complete.cases(d[, c(tc, e, terms)]) & d[[tc]] > 0, ]
  if ("stage" %in% terms) d$stage <- droplevels(d$stage)
  d$score <- zs(d$score); if ("prolif" %in% terms) d$prolif <- zs(d$prolif)
  s <- summary(coxph(as.formula(paste0("Surv(", tc, ", ", e, ") ~ ", paste(terms, collapse = " + "))), data = d))
  data.frame(n = s$n, events = s$nevent, HR = s$conf.int["score", 1], CI_low = s$conf.int["score", 3],
             CI_high = s$conf.int["score", 4], p = s$coefficients["score", 5])
}
cohorts <- list("All patients" = df, "M0 only" = df[!is.na(df$M) & df$M == "M0", ])
models  <- list("Unadjusted" = "score", "+ proliferation" = c("score", "prolif"),
                "+ clinical covariates" = c("score", COVARS),
                "+ proliferation + clinical covariates" = c("score", "prolif", COVARS))
tab_surv <- bind_rows(lapply(c("C1", "C2"), function(cl) bind_rows(lapply(names(cohorts), function(co)
  bind_rows(lapply(ENDPOINTS, function(e) bind_rows(lapply(names(models), function(m) {
    d <- cohorts[[co]] |> mutate(score = .data[[s_main(cl)]], prolif = Prolif)
    cbind(cluster = LAB[cl], cohort = co, endpoint = e, model = m, cox_hr(d, e, models[[m]]))
  }))))))))
write.csv(tab_surv, file.path(OUT_DIR, "Table_prolif_survival.csv"), row.names = FALSE)
cat("\n=== (b) Cox, HR por DP do escore do cluster ===\n")
print(tab_surv |> transmute(cluster = sub(" cluster", "", cluster), cohort, endpoint, model, n, events,
                            HR = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high), p = signif(p, 2)), row.names = FALSE)

## ---------------------------------------------------------------- 6. (c) análises principais sem os TFs discordantes
# análises principais do manuscrito: OR de metástase por DP (não ajustado) e Cox multivariável (COVARS)
tab_excl <- bind_rows(lapply(names(VARIANTS), function(v) bind_rows(lapply(c("C1", "C2"), function(cl) {
  col <- paste(cl, v, sep = "|")
  d <- dm |> mutate(score = zs(.data[[col]])); r <- logit_or(d, "score")
  met <- data.frame(variant = v, cluster = LAB[cl], analysis = "Metastasis (M1 vs M0), OR per SD",
                    cohort = "All patients", n = nrow(d), events = sum(d$m1),
                    estimate = r$or, CI_low = r$lo, CI_high = r$hi, p = r$p)
  sv <- bind_rows(lapply(names(cohorts), function(co) bind_rows(lapply(ENDPOINTS, function(e) {
    h <- cox_hr(cohorts[[co]] |> mutate(score = .data[[col]]), e, c("score", COVARS))
    data.frame(variant = v, cluster = LAB[cl], analysis = paste0(e, ", HR per SD (multivariable)"),
               cohort = co, n = h$n, events = h$events, estimate = h$HR, CI_low = h$CI_low,
               CI_high = h$CI_high, p = h$p)
  }))))
  bind_rows(met, sv)
}))))
write.csv(tab_excl, file.path(OUT_DIR, "Table_exclusion_sensitivity.csv"), row.names = FALSE)
cat("\n=== (c) Sensibilidade: exclusão de reguladores ===\n")
print(tab_excl |> transmute(variant, cluster = sub(" cluster", "", cluster), analysis, cohort,
                            est = sprintf("%.2f (%.2f-%.2f)", estimate, CI_low, CI_high), p = signif(p, 2)),
      row.names = FALSE)

## ---------------------------------------------------------------- 7. tabela suplementar + figura
wb <- createWorkbook()
for (sh in list(list("Correlations", tab_cor), list("Metastasis", tab_met),
                list("Survival", tab_surv), list("Exclusion_sensitivity", tab_excl))) {
  addWorksheet(wb, sh[[1]]); writeData(wb, sh[[1]], sh[[2]])
  setColWidths(wb, sh[[1]], cols = seq_len(ncol(sh[[2]])), widths = "auto")
}
saveWorkbook(wb, file.path(OUT_DIR, "TableS_comment4.xlsx"), overwrite = TRUE)

long <- bind_rows(lapply(c("C1", "C2"), function(cl)
  data.frame(cluster = LAB[cl], score = df[[s_main(cl)]], Prolif = df$Prolif, M = df$M)))
ann <- tab_cor |> filter(grepl("vs Prolif$", comparison)) |>
  mutate(cluster = sub(" vs Prolif$", "", comparison), txt = sprintf("ρ = %.2f", rho))
g <- ggplot(long, aes(Prolif, score)) +
  geom_point(aes(color = M), size = 0.8, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "#0b0b0b", linewidth = 0.5, formula = y ~ x) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt), hjust = -0.2, vjust = 1.5, size = 3,
            inherit.aes = FALSE) +
  scale_color_manual(values = c(M0 = "#2b6cb0", M1 = "#d73027"), na.value = "grey70", name = NULL) +
  facet_wrap(~cluster, scales = "free_y") +
  labs(x = "Proliferation score (Hallmark E2F targets + G2M checkpoint)",
       y = "Cluster activity score (mean dES)") +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "FigS_prolif_correlation.pdf"), g, width = 6.5, height = 3.2)
ggsave(file.path(OUT_DIR, "FigS_prolif_correlation.png"), g, width = 6.5, height = 3.2, dpi = 300, bg = "white")
cat("\nSaídas em", OUT_DIR, "\n")
