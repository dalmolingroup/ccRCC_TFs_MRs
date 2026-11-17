### Sobrevivência com atividade de regulons (RTN) e desfechos do TCGA-CDR ###
# Revisão CSBJ-D-26-01039 (Bloco A)
#   R2.M4 / R3.6 : Cox multivariável, HR [IC 95%], FDR, teste de riscos proporcionais
#   R1.M5        : PFI além de OS; análise só nos pacientes M0 (evita a circularidade M1/M0)
#   R2 Fig 6     : Kaplan-Meier com número em risco e eventos na legenda
#
# Entradas
#   - Objeto TNI do RTN (o mesmo do artigo). A atividade de cada regulon por amostra
#     (dES) vem de tni.gsea2(); se o objeto já tiver rodado tni.gsea2, ela é reaproveitada.
#   - TCGA-CDR (Liu et al., Cell 2018), Tabela S1, aba "TCGA-CDR":
#     https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
#   - metadata_KIRC_full.csv, só para o status M (análise nos M0)
#
# Abordagem: o escore de atividade entra no Cox como variável CONTÍNUA (HR por
# desvio padrão), sem ponto de corte. A divisão em alta/baixa atividade (mediana)
# serve só para as curvas de Kaplan-Meier.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(survival)
  library(survminer)
  library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG
RTN_FILE    <- "0_data/rtni_TF.RData"
RTN_OBJECT  <- "rtni_TF"          # nome do objeto dentro do .RData
ACTIVITY_RDS <- NULL                           # opcional: matriz amostras x regulons (dES) já salva;
# se preenchido, o RTN não é usado
CDR_FILE    <- "0_data/TCGA-CDR-SupplementalTableS1.xlsx"
META_FILE   <- "0_data/metadata_KIRC_full.csv"
META_ID_COL <- "sample"
M_COL       <- "ajcc_pathologic_m"
OUT_DIR     <- "results_survival"

# Clusters de regulons por SÍMBOLO (convertidos para os nomes do RTN automaticamente).
# O primeiro é o teste primário; os escores dos clusters ficam fora da correção de BH dos TFs.
CLUSTERS <- list(
  "Differentiation-associated cluster"        = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  "Proliferative cluster" = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A")
)
CLUSTER_COLORS <- c(differentiation_associated_cluster = "#1a9850", proliferative_cluster = "#542788")

# Conjuntos de MRs testados um a um (forest plot + tabela com FDR); arquivos ausentes são ignorados
MR_FILES <- c(signature = "results_ml_revision/importancia_features_modelo_final_class_weight.csv")   # assinatura única do artigo (14 TFs): .txt (um ID por linha) ou .csv (IDs na 1ª coluna)

ENDPOINTS <- c("OS", "PFI", "DSS")             # CDR recomenda OS, PFI e DSS para KIRC; DFI tem poucos eventos
SYMBOL_FILE <- "0_data/ensembl_to_symbol.csv"  # opcional: colunas ensembl,symbol (rótulos das figuras)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- 1. atividade dos regulons
# ---- identificadores: tudo é convertido para os nomes dos regulons no RTN ----
# Aceita ENSG (com ou sem versão) ou símbolo, e funciona se o RTN usar ENSG ou símbolo.
reg_universe <- function() {
  if (!is.null(ACTIVITY_RDS)) return(colnames(readRDS(ACTIVITY_RDS)))
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(RTN_FILE, envir = env)
  tni.get(get(RTN_OBJECT, envir = env), what = "regulatoryElements")
}
id_key <- function(x) toupper(sub("\\..*", "", trimws(gsub('"', "", as.character(x)))))
to_reg <- function(ids, regs) {
  vals <- unname(regs); nm <- names(regs)
  s <- if (file.exists(SYMBOL_FILE)) read.csv(SYMBOL_FILE) else data.frame(ensembl = character(), symbol = character())
  ens2sym <- setNames(toupper(s$symbol), id_key(s$ensembl))
  sym2ens <- setNames(id_key(s$ensembl), toupper(s$symbol))
  out <- vapply(ids, function(id) {
    k <- id_key(id)
    cand <- unique(na.omit(c(k, unname(ens2sym[k]), unname(sym2ens[k]),     # o próprio ID, o símbolo e o ENSG
                             unname(REG_V2N[k]), unname(REG_N2V[k]))))      # e o par ID/nome do próprio RTN
    for (pool in list(vals, nm)) {
      if (is.null(pool)) next
      hit <- match(cand, id_key(pool)); hit <- hit[!is.na(hit)]
      if (length(hit)) return(vals[hit[1]])
    }
    NA_character_
  }, character(1))
  if (any(is.na(out))) stop("Não encontrados no RTN: ", paste(ids[is.na(out)], collapse = ", "),
                            "\nExemplos de nomes no RTN: ", paste(head(vals), collapse = ", "))
  unname(out)
}
# lê a assinatura: .txt (um ID por linha) ou .csv (IDs na 1ª coluna, ex.: importância do modelo final)
read_tf_file <- function(f) {
  if (!file.exists(f)) return(NULL)
  ids <- if (grepl("\\.csv$", f, ignore.case = TRUE)) as.character(read.csv(f, check.names = FALSE)[[1]]) else readLines(f, warn = FALSE)
  ids <- trimws(gsub('"', "", ids))
  ids[ids != "" & !tolower(ids) %in% c("x", "ensembl", "gene", "feature", "id")]
}
CLUSTER_VARS  <- gsub("[^a-z]+", "_", tolower(names(CLUSTERS)))          # differentiation_associated_cluster, proliferative_cluster
CLUSTER_LABEL <- setNames(names(CLUSTERS), CLUSTER_VARS)
regs_all    <- reg_universe()
# o RTN guarda ID (ex.: ENSG) e nome (ex.: símbolo); a atividade (dES) usa os nomes
REG_V2N <- if (!is.null(names(regs_all))) setNames(id_key(names(regs_all)), id_key(regs_all)) else character(0)
REG_N2V <- if (!is.null(names(regs_all))) setNames(id_key(regs_all), id_key(names(regs_all))) else character(0)
CLUSTER_REG <- setNames(lapply(CLUSTERS, to_reg, regs = regs_all), CLUSTER_VARS)
MR_SETS <- lapply(MR_FILES, read_tf_file)
MR_SETS <- MR_SETS[!vapply(MR_SETS, is.null, logical(1))]
MR_SETS <- lapply(MR_SETS, to_reg, regs = regs_all)
message("Assinatura: ", paste(lengths(MR_SETS), "TFs", collapse = "; "), " | clusters: ",
        paste(lengths(CLUSTER_REG), collapse = " + "), " TFs")
tfs_needed <- unique(c(unlist(CLUSTER_REG), unlist(MR_SETS)))

get_activity <- function() {
  if (!is.null(ACTIVITY_RDS)) return(readRDS(ACTIVITY_RDS))
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(RTN_FILE, envir = env)
  rtni <- get(RTN_OBJECT, envir = env)
  act <- tryCatch(tni.get(rtni, what = "regulonActivity"), error = function(e) NULL)
  have <- if (!is.null(act$differential)) colnames(act$differential) else character(0)
  ok <- length(have) > 0 && !inherits(try(to_reg(tfs_needed, have), silent = TRUE), "try-error")
  if (!ok) {
    message("Rodando tni.gsea2 para ", length(tfs_needed), " regulons (pode levar alguns minutos)...")
    rtni <- tni.gsea2(rtni, regulatoryElements = tfs_needed)
    act <- tni.get(rtni, what = "regulonActivity")
  }
  act$differential                           # amostras x regulons (dES)
}

dES <- get_activity()
# converte clusters e assinatura para os nomes das colunas da atividade (ex.: ENSG -> símbolo)
CLUSTER_REG <- lapply(CLUSTER_REG, to_reg, regs = colnames(dES))
MR_SETS     <- lapply(MR_SETS, to_reg, regs = colnames(dES))
tfs_needed  <- unique(c(unlist(CLUSTER_REG), unlist(MR_SETS)))
missing_tfs <- setdiff(tfs_needed, colnames(dES))
if (length(missing_tfs))
  stop("Regulons sem atividade: ", paste(missing_tfs, collapse = ", "),
       "\nConfira se os nomes batem com os do RTN, ex.: ", paste(head(colnames(dES)), collapse = ", "))

symbols <- if (file.exists(SYMBOL_FILE)) {
  s <- read.csv(SYMBOL_FILE); setNames(s$symbol, sub("\\..*", "", s$ensembl))
} else character(0)
lab <- function(g) ifelse(sub("\\..*", "", g) %in% names(symbols), symbols[sub("\\..*", "", g)], g)

# 1 tumor primário por paciente
act <- as.data.frame(dES[, tfs_needed, drop = FALSE])
act$sample <- rownames(dES)
act <- act |>
  filter(substr(sample, 14, 15) == "01") |>
  mutate(patient = substr(sample, 1, 12)) |>
  distinct(patient, .keep_all = TRUE)
# escore de cada cluster = média dos dES dos TFs do cluster
for (cv in CLUSTER_VARS) act[[cv]] <- rowMeans(act[, CLUSTER_REG[[cv]], drop = FALSE])

## ---------------------------------------------------------------- 2. desfechos e covariáveis (CDR)
cdr <- read_excel(CDR_FILE, sheet = "TCGA-CDR", na = c("", "NA", "#N/A", "[Not Available]",
                                                       "[Not Applicable]", "[Unknown]",
                                                       "[Discrepancy]", "[Not Evaluated]"))
need <- c("bcr_patient_barcode", "type", "age_at_initial_pathologic_diagnosis", "gender",
          "ajcc_pathologic_tumor_stage", "histological_grade",
          as.vector(rbind(ENDPOINTS, paste0(ENDPOINTS, ".time"))))
if (length(setdiff(need, names(cdr))))
  stop("Colunas ausentes no CDR: ", paste(setdiff(need, names(cdr)), collapse = ", "))

clin <- cdr |>
  filter(type == "KIRC") |>
  transmute(
    patient = bcr_patient_barcode,
    age   = as.numeric(age_at_initial_pathologic_diagnosis),
    sex   = factor(tolower(gender), levels = c("female", "male")),
    stage = case_when(grepl("Stage IV", ajcc_pathologic_tumor_stage) ~ "IV",
                      grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
                      grepl("Stage I", ajcc_pathologic_tumor_stage) ~ "I-II"),   # I e II
    grade = case_when(histological_grade %in% c("G1", "G2") ~ "G1-G2",
                      histological_grade == "G3" ~ "G3",
                      histological_grade == "G4" ~ "G4"),
    across(all_of(as.vector(rbind(ENDPOINTS, paste0(ENDPOINTS, ".time")))), as.numeric)
  ) |>
  mutate(stage = factor(stage, levels = c("I-II", "III", "IV")),
         grade = factor(grade, levels = c("G1-G2", "G3", "G4")))

meta <- read.csv(META_FILE, check.names = FALSE)
meta <- data.frame(patient = substr(meta[[META_ID_COL]], 1, 12), M = meta[[M_COL]]) |>
  distinct(patient, .keep_all = TRUE)

df <- act |> inner_join(clin, by = "patient") |> left_join(meta, by = "patient")
write.csv(df, file.path(OUT_DIR, "dados_sobrevivencia.csv"), row.names = FALSE)


flow <- c(regulon_activity_primary_tumor = nrow(act), with_CDR = nrow(df),
          M0 = sum(df$M == "M0", na.rm = TRUE), M1 = sum(df$M == "M1", na.rm = TRUE))
for (e in ENDPOINTS) {
  ok <- !is.na(df[[e]]) & !is.na(df[[paste0(e, ".time")]]) & df[[paste0(e, ".time")]] > 0
  flow[paste0(e, "_n")] <- sum(ok); flow[paste0(e, "_events")] <- sum(df[[e]][ok] == 1)
}
print(flow)
write.csv(data.frame(step = names(flow), n = flow), file.path(OUT_DIR, "fluxo_sobrevivencia.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------- 3. Cox
COVARS <- c("age", "sex", "stage", "grade")

fit_cox <- function(data, endpoint, var, covars = character(0)) {
  tcol <- paste0(endpoint, ".time")
  d <- data[, c(tcol, endpoint, var, covars)]
  d <- d[complete.cases(d) & d[[tcol]] > 0, ]
  d$z <- as.numeric(scale(d[[var]]))                         # HR por 1 desvio padrão
  f <- as.formula(paste0("Surv(", tcol, ", ", endpoint, ") ~ ",
                         paste(c("z", covars), collapse = " + ")))
  fit <- coxph(f, data = d)
  s <- summary(fit)
  zph <- tryCatch(cox.zph(fit)$table, error = function(e) NULL)
  data.frame(endpoint = endpoint, variable = var,
             model = if (length(covars)) "multivariable" else "univariable",
             n = s$n, events = s$nevent,
             HR = s$conf.int["z", "exp(coef)"],
             CI_low = s$conf.int["z", "lower .95"], CI_high = s$conf.int["z", "upper .95"],
             p = s$coefficients["z", "Pr(>|z|)"],
             PH_p_variable = if (!is.null(zph)) zph["z", "p"] else NA,
             PH_p_global = if (!is.null(zph)) zph["GLOBAL", "p"] else NA)
}

vars <- c(CLUSTER_VARS, unique(unlist(MR_SETS)))
cohorts <- list(all = df, M0_only = df[!is.na(df$M) & df$M == "M0", ])

res <- list()
for (coh in names(cohorts)) for (e in ENDPOINTS) for (v in vars) {
  d <- cohorts[[coh]]
  # nos M0 o estágio IV quase não existe: ajusta por idade, sexo, grau e estágio I-III
  cv <- COVARS
  if (coh == "M0_only") d$stage <- droplevels(d$stage)
  res[[length(res) + 1]] <- cbind(cohort = coh, fit_cox(d, e, v))
  res[[length(res) + 1]] <- cbind(cohort = coh, fit_cox(d, e, v, cv))
}
cox <- bind_rows(res) |>
  mutate(label = ifelse(variable %in% CLUSTER_VARS, CLUSTER_LABEL[variable], lab(variable)),
         in_set = vapply(variable, function(v) paste(names(MR_SETS)[vapply(MR_SETS, \(s) v %in% s, TRUE)],
                                                     collapse = ";"), ""),
         # FDR entre os TFs testados, dentro de cada coorte x desfecho x modelo (o escore do cluster é o teste primário, fora da correção)
         FDR = NA_real_)
tf_rows <- !cox$variable %in% CLUSTER_VARS
cox$FDR[tf_rows] <- ave(cox$p[tf_rows], cox$cohort[tf_rows], cox$endpoint[tf_rows], cox$model[tf_rows],
                        FUN = function(p) p.adjust(p, "BH"))
write.csv(cox, file.path(OUT_DIR, "cox_resultados.csv"), row.names = FALSE)

# modelo completo do escore de cada cluster (todas as covariáveis), para tabela do texto
full_tab <- list()
for (cv in CLUSTER_VARS) for (coh in names(cohorts)) for (e in ENDPOINTS) {
  d <- cohorts[[coh]]; tcol <- paste0(e, ".time")
  d <- d[complete.cases(d[, c(tcol, e, cv, COVARS)]) & d[[tcol]] > 0, ]
  d$stage <- droplevels(d$stage)
  d$score_SD <- as.numeric(scale(d[[cv]]))
  fit <- coxph(as.formula(paste0("Surv(", tcol, ", ", e, ") ~ score_SD + ",
                                 paste(COVARS, collapse = " + "))), data = d)
  s <- summary(fit)
  full_tab[[length(full_tab) + 1]] <- data.frame(cluster = CLUSTER_LABEL[[cv]], cohort = coh, endpoint = e,
                                                 term = rownames(s$conf.int),
                                                 HR = s$conf.int[, 1], CI_low = s$conf.int[, 3],
                                                 CI_high = s$conf.int[, 4], p = s$coefficients[, 5],
                                                 n = s$n, events = s$nevent)
}
write.csv(bind_rows(full_tab), file.path(OUT_DIR, "cox_cluster_modelo_completo.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 4. figuras
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GRAY <- "#b8b6af"; INK2 <- "#52514e"
theme_set(theme_minimal(base_size = 9) +
            theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
                  axis.text = element_text(color = INK2), strip.text = element_text(face = "bold")))

# Forest plot: multivariável, todos os pacientes e só M0
fp <- cox |> filter(model == "multivariable") |>
  mutate(sig = ifelse(variable %in% CLUSTER_VARS, ifelse(p < 0.05, "p < 0.05", "n.s."),
                      ifelse(FDR < 0.05, "p < 0.05", "n.s.")),
         cohort = recode(cohort, all = "All patients", M0_only = "M0 only"),
         endpoint = factor(endpoint, levels = ENDPOINTS),
         label = factor(label, levels = rev(unique(c(unname(CLUSTER_LABEL), lab(unique(unlist(MR_SETS))))))))
g <- ggplot(fp, aes(HR, label, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_x_log10() +
  scale_color_manual(values = c("p < 0.05" = BLUE, "n.s." = GRAY), name = NULL) +
  facet_grid(cohort ~ endpoint) +
  labs(x = "Hazard ratio per SD of regulon activity (log scale)", y = NULL,
       caption = "Adjusted for age, sex, stage and grade. TFs: BH-adjusted p; cluster scores: nominal p.") +
  theme(legend.position = "bottom", plot.caption = element_text(color = INK2, hjust = 0))
ggsave(file.path(OUT_DIR, "forest_cox_multivariavel.pdf"), g, width = 7.5, height = 1.2 + 0.32 * nlevels(fp$label) * 2)
ggsave(file.path(OUT_DIR, "forest_cox_multivariavel.png"), g, width = 7.5, height = 1.2 + 0.32 * nlevels(fp$label) * 2, dpi = 300)

# Kaplan-Meier do escore do cluster (mediana), com número em risco
km_plot <- function(d, endpoint, title, file, var) {
  tcol <- paste0(endpoint, ".time")
  d <- d[!is.na(d[[endpoint]]) & !is.na(d[[tcol]]) & d[[tcol]] > 0, ]
  d$time_months <- d[[tcol]] / 30.44
  d$event <- d[[endpoint]]
  d$group <- factor(ifelse(d[[var]] > median(d[[var]]), "High", "Low"),
                    levels = c("Low", "High"))
  ev <- tapply(d$event, d$group, sum); nn <- table(d$group)
  legend_labs <- sprintf("%s (n = %d; %d events)", levels(d$group), nn, ev)
  fit <- survfit(Surv(time_months, event) ~ group, data = d)
  p <- ggsurvplot(fit, data = d, pval = TRUE, pval.size = 3, risk.table = TRUE,
                  risk.table.height = 0.25, palette = c(ORANGE, BLUE), legend.labs = legend_labs,
                  legend.title = "Activity", legend = "top", xlab = "Months", ylab = paste(endpoint, "probability"),
                  title = title, censor.size = 2, size = 0.7, break.time.by = 24,
                  ggtheme = theme_classic(base_size = 9), tables.theme = theme_cleantable(base_size = 8),
                  fontsize = 3, tables.y.text = FALSE)
  for (ext in c("pdf", "png")) {
    f <- file.path(OUT_DIR, paste0(file, ".", ext))
    if (ext == "pdf") pdf(f, width = 5, height = 4.6, onefile = FALSE) else png(f, width = 5, height = 4.6, units = "in", res = 300)
    print(p); dev.off()
  }
}
for (cv in CLUSTER_VARS) for (e in ENDPOINTS) {
  km_plot(df, e, paste0(CLUSTER_LABEL[[cv]], ", ", e, ", all patients"), paste0("km_", cv, "_", e, "_todos"), cv)
  km_plot(cohorts$M0_only, e, paste0(CLUSTER_LABEL[[cv]], ", ", e, ", M0 patients"), paste0("km_", cv, "_", e, "_M0"), cv)
}

## ---------------------------------------------------------------- 5. resumo no console
cat("\n=== Escores dos clusters: Cox multivariável (HR por DP) ===\n")
print(cox |> filter(variable %in% CLUSTER_VARS, model == "multivariable") |>
        transmute(cluster = label, cohort, endpoint, n, events, HR = round(HR, 2),
                  CI = sprintf("%.2f-%.2f", CI_low, CI_high), p = signif(p, 3),
                  PH_p = signif(PH_p_variable, 2)))
cat("\nPH_p < 0.05 indica violação de riscos proporcionais: reporte e considere estratificar a covariável.\n")
cat("Saídas em", OUT_DIR, "\n")
