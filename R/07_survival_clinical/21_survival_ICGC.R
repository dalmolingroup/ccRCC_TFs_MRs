### Sobrevivência no ICGC-RECA, dois clusters — R2.M4 ###
# "Independent validation of these survival associations would further strengthen the findings."
#
# 1. Calcula a atividade dos 14 regulons em TODOS os tumores primários (RT) do ICGC-RECA, 1 por doador
#    (91 tumores; o pareamento com RA só é necessário nas comparações tumor x normal): regulons da rede do TCGA + tni.replace.samples (voom) + tni.gsea2.
#    As amostras de tecido adjacente (RA) são excluídas ANTES da normalização e do tni.gsea2, porque
#    a atividade é relativa à média da coorte (correção da revisão: a versão anterior usava 45 RT + 45 RA).
#    Salva dES_ICGC_RT91.csv e surv_ICGC_RT91.csv - fonte ÚNICA da atividade do ICGC, usada também
#    na Fig. S5 (12_FigS5_ICGC_activity.R) e na associação com metástase (11b_activity_vs_M_ICGC.R).
# 2. Para CADA cluster: escore = média dos dES; Cox de OS em camadas (o ICGC não tem grau nem estágio):
#    não ajustado | idade + sexo | idade + pT | idade + pT + M | só M0
# 3. KM (mediana) para visualização
#
# Saídas (OUT_DIR): cox_ICGC_clusters.csv (lido pelo 22_Fig7_FigS7-S9_TabS13-S25.R), KMs por cluster

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(survminer); library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG
RTN_FILE     <- "~/Kidney/0_data/rtni_TF.RData"; RTN_OBJECT <- "rtni_TF"
ICGC_RDATA   <- "~/Kidney/0_data/ICGC_DEA_Data.RData"           # count_final_ICGC, donor_data_final
ACTIVITY_CSV <- "~/Kidney/0_data/dES_ICGC_RT91.csv"             # amostras (91 tumores RT) x regulons
SURVIVAL_CSV <- "~/Kidney/0_data/surv_ICGC_RT91.csv"
META_FILE    <- "~/Kidney/0_data/metadata_ICGC_full.csv"        # covariáveis que faltem
FORCE_RECOMPUTE <- TRUE      # TRUE: recalcula a atividade mesmo se os CSVs existirem (evita cache antigo)
N_EXPECTED   <- 91           # nº de tumores RT com expressão no ICGC-RECA (trava de segurança)
TUMOR_VALUE  <- "RT"
TIME_COL <- "obs.time"; EVENT_COL <- "status"; EVENT_VALUE <- "deceased"
PT_COL <- "pathologic_t"; M_COL <- "pathologic_m"
TIME_UNIT <- "days"
OUT_DIR  <- "~/Kidney/results_survival_ICGC"

CLUSTERS <- list(
  "Proliferative cluster"              = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A"),
  "Differentiation-associated cluster" = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"))
CL_COL <- c("Proliferative cluster" = "#542788", "Differentiation-associated cluster" = "#1a9850")

pe <- path.expand
dir.create(pe(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
strip_v <- function(x) sub("\\..*", "", x)

# Checagem: só tumores (sufixo RT no ID = doador + tipo de tecido) e 1 amostra por doador
check_samples <- function(ids, where) {
  tissue <- substr(ids, nchar(ids) - 1, nchar(ids))
  donor  <- substr(ids, 1, nchar(ids) - 2)
  if (!all(tissue == TUMOR_VALUE))
    stop(where, ": há amostras que não são tumor (", paste(unique(tissue[tissue != TUMOR_VALUE]), collapse = ", "), ")")
  if (anyDuplicated(donor)) stop(where, ": há doadores duplicados")
  if (length(ids) != N_EXPECTED)
    warning(sprintf("%s: %d amostras (esperado %d) - confira", where, length(ids), N_EXPECTED))
  invisible(TRUE)
}

## ---------------------------------------------------------------- 1. atividade no ICGC (só RT)
if (FORCE_RECOMPUTE || !file.exists(pe(ACTIVITY_CSV)) || !file.exists(pe(SURVIVAL_CSV))) {
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(pe(RTN_FILE), envir = env); rtni <- get(RTN_OBJECT, envir = env)
  load(pe(ICGC_RDATA))
  don <- donor_data_final
  don$sample <- paste0(don$submitted_donor_id, don$tissue_type)
  cat("Amostras no ICGC por tipo de tecido:\n"); print(table(don$tissue_type, useNA = "ifany"))

  # todos os tumores primários (RT) com expressão, 1 por doador (91); os RA ficam fora
  # (só servem para as comparações tumor x normal)
  don <- don[don$tissue_type == TUMOR_VALUE & don$sample %in% colnames(count_final_ICGC), ]
  don <- don[!duplicated(don$submitted_donor_id), ]
  check_samples(don$sample, "seleção de amostras")

  # contagens SÓ dos tumores -> voom -> regulons do TCGA -> atividade relativa à coorte de tumores
  cnt <- count_final_ICGC; rownames(cnt) <- strip_v(cnt$gene_id); cnt$gene_id <- NULL
  cnt <- as.matrix(cnt[, don$sample])
  expr <- GDCRNATools::gdcVoomNormalization(cnt, filter = FALSE)
  expr <- expr[rownames(expr) %in% rownames(tni.get(rtni, "gexp")), ]
  rtni_icgc <- tni.replace.samples(rtni, expData = expr)
  ids <- tni.get(rtni, "regulatoryElements")                        # nomes = símbolos, valores = IDs
  tfs <- unname(ids[names(ids) %in% unlist(CLUSTERS)])
  rtni_icgc <- tni.gsea2(rtni_icgc, regulatoryElements = tfs)
  write.csv(tni.get(rtni_icgc, what = "regulonActivity")$differential, pe(ACTIVITY_CSV))
  surv <- don[, c("sample", TIME_COL, EVENT_COL, PT_COL, M_COL, "gender", grep("age", names(don), value = TRUE)[1])]
  rownames(surv) <- surv$sample; surv$sample <- NULL       # (subset herdava rownames: column_to_rownames falhava)
  write.csv(surv, pe(SURVIVAL_CSV))
  message("Atividade salva em ", ACTIVITY_CSV)
}

## ---------------------------------------------------------------- 2. dados
dES  <- as.matrix(read.csv(pe(ACTIVITY_CSV), row.names = 1, check.names = FALSE))
surv <- read.csv(pe(SURVIVAL_CSV), row.names = 1, check.names = FALSE); surv$sample <- rownames(surv)
check_samples(rownames(dES), ACTIVITY_CSV)       # trava também para CSVs reaproveitados
check_samples(surv$sample, SURVIVAL_CSV)

# colunas da atividade -> símbolo (aceita símbolo ou ENSG)
sym_file <- "~/Kidney/0_data/ensembl_to_symbol.csv"
if (!all(unlist(CLUSTERS) %in% colnames(dES)) && file.exists(pe(sym_file))) {
  s <- read.csv(pe(sym_file)); m <- setNames(s$symbol, strip_v(s$ensembl))
  colnames(dES) <- ifelse(strip_v(colnames(dES)) %in% names(m), m[strip_v(colnames(dES))], colnames(dES))
}
miss <- setdiff(unlist(CLUSTERS), colnames(dES))
if (length(miss)) stop("TFs sem atividade no ICGC: ", paste(miss, collapse = ", "),
                       "\nColunas: ", paste(head(colnames(dES), 10), collapse = ", "))

base <- inner_join(data.frame(sample = rownames(dES)), surv, by = "sample")
if (file.exists(pe(META_FILE))) {
  meta <- read.csv(pe(META_FILE), check.names = FALSE, row.names = 1); meta$sample <- rownames(meta)
  base <- left_join(base, meta[, c("sample", setdiff(names(meta), names(base)))], by = "sample")
}
age_col <- names(base)[grepl("age", names(base), ignore.case = TRUE)][1]
sex_col <- names(base)[grepl("gender|^sex$", names(base), ignore.case = TRUE)][1]
ev <- trimws(tolower(as.character(base[[EVENT_COL]])))
base$event <- if (is.numeric(base[[EVENT_COL]]) && all(base[[EVENT_COL]] %in% c(1, 2, NA))) as.integer(base[[EVENT_COL]] == 2) else
  if (is.numeric(base[[EVENT_COL]])) as.integer(base[[EVENT_COL]] == 1) else
    as.integer(ev %in% c(tolower(EVENT_VALUE), "deceased", "dead", "died"))
base$time <- as.numeric(base[[TIME_COL]])
base$age  <- as.numeric(base[[age_col]])
base$sex  <- factor(tolower(base[[sex_col]]))
base$pT   <- factor(ifelse(grepl("T[12]", base[[PT_COL]]), "T1-T2", ifelse(grepl("T[34]", base[[PT_COL]]), "T3-T4", NA)),
                    levels = c("T1-T2", "T3-T4"))
base$M    <- factor(ifelse(base[[M_COL]] %in% c("M0", "M1"), base[[M_COL]], NA), levels = c("M0", "M1"))
n_before <- nrow(base)
base <- base[!is.na(base$time) & base$time > 0 & !is.na(base$event), ]
cat(sprintf("ICGC: %d tumores (%d excluídos por falta de tempo/evento), %d óbitos (M0 = %d, M1 = %d, M ausente = %d)\n",
            nrow(base), n_before - nrow(base), sum(base$event),
            sum(base$M == "M0", na.rm = TRUE), sum(base$M == "M1", na.rm = TRUE), sum(is.na(base$M))))

## ---------------------------------------------------------------- 3. Cox em camadas, por cluster
fit <- function(d, covars, label, cohort) {
  d <- d[complete.cases(d[, c("time", "event", "z", covars)]), ]
  m <- coxph(as.formula(paste("Surv(time, event) ~", paste(c("z", covars), collapse = " + "))), data = d)
  s <- summary(m)
  data.frame(cohort = cohort, model = label, n = s$n, events = s$nevent,
             HR = s$conf.int["z", 1], CI_low = s$conf.int["z", 3], CI_high = s$conf.int["z", 4],
             p = s$coefficients["z", 5], PH_p = tryCatch(cox.zph(m)$table["z", "p"], error = function(e) NA),
             events_per_variable = s$nevent / (1 + length(covars)))
}
models <- list("Unadjusted" = character(0), "Age + sex" = c("age", "sex"),
               "Age + pT" = c("age", "pT"), "Age + pT + M" = c("age", "pT", "M"))

res <- list()
for (cl in names(CLUSTERS)) {
  df <- base
  df$score <- rowMeans(dES[df$sample, CLUSTERS[[cl]], drop = FALSE])
  df$z <- as.numeric(scale(df$score))                       # HR por DP (padronizado nos tumores do ICGC)
  r <- bind_rows(lapply(names(models), function(n) fit(df, models[[n]], n, "All patients")))
  m0 <- filter(df, M == "M0")
  r <- bind_rows(r, fit(m0, character(0), "Unadjusted", "M0 only"), fit(m0, "age", "Age", "M0 only"))
  res[[cl]] <- cbind(cluster = cl, r)

  # KM (mediana), só para visualização
  km <- function(d, title, file) {
    d$t <- if (TIME_UNIT == "days") d$time / 30.44 else d$time
    d$group <- factor(ifelse(d$score > median(d$score), "High", "Low"), levels = c("High", "Low"))
    nn <- table(d$group); ee <- tapply(d$event, d$group, sum)
    p <- ggsurvplot(survfit(Surv(t, event) ~ group, data = d), data = d, pval = TRUE, pval.size = 3,
                    risk.table = TRUE, risk.table.height = 0.25, palette = c(unname(CL_COL[cl]), "#9a9893"),
                    legend.labs = sprintf("%s (n = %d; %d deaths)", levels(d$group), nn, ee),
                    legend.title = "Activity", legend = c(0.7, 0.88), xlab = "Months",
                    ylab = "Overall survival probability", title = title, censor.size = 2, size = 0.7,
                    ggtheme = theme_classic(base_size = 9), tables.theme = theme_cleantable(base_size = 8),
                    fontsize = 3, tables.y.text = FALSE)
    for (ext in c("pdf", "png")) {
      f <- file.path(pe(OUT_DIR), paste0(file, ".", ext))
      if (ext == "pdf") pdf(f, width = 5, height = 4.6, onefile = FALSE) else png(f, width = 5, height = 4.6, units = "in", res = 300)
      print(p); dev.off()
    }
  }
  tag <- gsub("[^a-z]+", "_", tolower(cl))
  km(df, paste0(cl, ", OS, ICGC-RECA"), paste0("km_", tag, "_OS_ICGC"))
  if (sum(m0$event) >= 5) km(m0, paste0(cl, ", OS, ICGC-RECA M0"), paste0("km_", tag, "_OS_ICGC_M0"))
}
res <- bind_rows(res)
write.csv(res, file.path(pe(OUT_DIR), "cox_ICGC_clusters.csv"), row.names = FALSE)

cat("\n=== Escores dos clusters no ICGC (Cox de OS, HR por DP; 91 tumores RT) ===\n")
print(res |> transmute(cluster, cohort, model, n, events, HR = round(HR, 2),
                       CI = sprintf("%.2f-%.2f", CI_low, CI_high), p = signif(p, 3),
                       PH_p = signif(PH_p, 2), EPV = round(events_per_variable, 1)), row.names = FALSE)
cat("\nEPV < 10: poucos eventos por variável; interpretar com cautela.\n")
