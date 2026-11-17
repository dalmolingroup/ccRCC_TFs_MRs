### Tabela 2 do texto principal: TFs da assinatura (R2.a15) ###
# Revisão CSBJ-D-26-01039
# Junta as saídas que já existem; não refaz nenhuma análise.
#
# Entradas (rode antes os scripts correspondentes, com a assinatura de 14 TFs em MR_FILES):
#   figure5_expression.R    -> SuppTable_expression_tumor_vs_normal.csv, SuppTable_expression_M1_vs_M0.csv
#   activity_vs_clinical.R  -> atividade_vs_clinico.csv
#   survival_cox_cdr.R      -> cox_resultados.csv
#   (opcional) tabela da MRA do RTN (tna.get(rtna, what = "mra")), colunas Regulon e Adjusted.Pvalue
#
# Saídas: Table2_signature_TFs.csv (números completos) e Table2_signature_TFs.tex (formatada)

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

## ---------------------------------------------------------------- CONFIG
FIG5_DIR  <- "~/Kidney/Figures/Figure5"
ACT_FILE  <- "results_activity_clinical/atividade_vs_clinico.csv"
COX_FILE  <- "results_survival/cox_resultados.csv"
MRA_FILE  <- "~/Kidney/0_data/MRA.RData"   # .RData com o objeto TNA, ou o CSV; NULL = coluna omitida
OUT_DIR   <- "results_table2"
SYMBOL_FILE <- "~/Kidney/0_data/ensembl_to_symbol.csv"   # gerado pelo figure5_expression.R (colunas ensembl, symbol)

GRADE_FEATURE <- "Grade"     # associação de atividade mostrada (Grade, pT ou Stage)
COX_COLS <- list(            # colunas de Cox: rótulo = c(coorte, desfecho)
  "DSS, all patients" = c("all", "DSS"),
  "PFI, M0 patients"  = c("M0_only", "PFI")
)
CLUSTER_ORDER <- c("Proliferative cluster", "Differentiation-associated cluster")
# nomes das variáveis dos escores nas saídas de survival_cox_cdr.R e activity_vs_clinical.R
CLUSTER_ROWS_LABEL <- c(proliferative_cluster = "Proliferative cluster",
                        differentiation_associated_cluster = "Differentiation-associated cluster")
# nomes antigos/numéricos dos escores nas saídas -> nomes acima
CLUSTER_ALIAS <- c(cluster_1 = "proliferative_cluster", cluster_2 = "differentiation_associated_cluster",
                   immune_cluster = "differentiation_associated_cluster")
fix_cl <- function(v) ifelse(v %in% names(CLUSTER_ALIAS), CLUSTER_ALIAS[v], v)
CLUSTER_ROWS <- names(CLUSTER_ROWS_LABEL)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
pe <- path.expand
strip <- function(x) sub("\\..*", "", x)
fp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                              formatC(p, format = "f", digits = 3)))

## ---------------------------------------------------------------- expressão (Fig 5)
tn <- read.csv(file.path(pe(FIG5_DIR), "SuppTable_expression_tumor_vs_normal.csv"))
mm <- read.csv(file.path(pe(FIG5_DIR), "SuppTable_expression_M1_vs_M0.csv"))

tab <- tn |> filter(cohort == "TCGA-KIRC") |>
  transmute(ensembl = strip(ensembl), TF = gene, cluster,
            TvN_log2FC = median_log2FC_paired, TvN_r = rank_biserial_paired, TvN_padj = p_adj_paired) |>
  left_join(mm |> filter(cohort == "TCGA-KIRC") |>
              transmute(ensembl = strip(ensembl), M1vM0_delta = cliffs_delta, M1vM0_padj = p_adj),
            by = "ensembl") |>
  left_join(mm |> filter(cohort == "ICGC-RECA") |>
              transmute(ensembl = strip(ensembl), ICGC_delta = cliffs_delta),
            by = "ensembl") |>
  mutate(ICGC_concordant = ifelse(is.na(ICGC_delta), NA, sign(ICGC_delta) == sign(M1vM0_delta)))

## ---------------------------------------------------------------- identificadores
# Todas as fontes são convertidas para ENSG sem versão antes de qualquer junção.
# Aceita ENSG (com ou sem versão) ou símbolo. Para os TFs da assinatura, o símbolo é
# convertido no ENSG usado na Fig 5 (evita símbolos com mais de um ENSG).
sym2ens <- setNames(tab$ensembl, toupper(tab$TF))
if (file.exists(pe(SYMBOL_FILE))) {
  s <- read.csv(pe(SYMBOL_FILE))
  s <- s[!is.na(s$symbol) & s$symbol != "" & !toupper(s$symbol) %in% names(sym2ens), ]
  s <- s[!duplicated(toupper(s$symbol)), ]
  sym2ens <- c(sym2ens, setNames(strip(s$ensembl), toupper(s$symbol)))
}
to_ens <- function(x) {
  x <- trimws(as.character(x))
  out <- ifelse(grepl("^ENSG", x, ignore.case = TRUE), strip(toupper(x)), unname(sym2ens[toupper(x)]))
  out[grepl("_cluster$|^cluster_score$", x)] <- x[grepl("_cluster$|^cluster_score$", x)]
  out
}
check_ids <- function(ids, fonte) {
  bad <- unique(ids[is.na(to_ens(ids))])
  if (length(bad)) message("[aviso] ", fonte, ": ", length(bad), " identificador(es) sem ENSG: ",
                           paste(head(bad, 10), collapse = ", "))
}

## ---------------------------------------------------------------- atividade vs grau
act <- read.csv(ACT_FILE) |> filter(cohort == "all", feature == GRADE_FEATURE) |> mutate(variable = fix_cl(variable))
check_ids(act$variable, "atividade_vs_clinico")
act_tf <- act |> filter(!variable %in% CLUSTER_ROWS, variable != "cluster_score") |>
  transmute(ensembl = to_ens(variable), act_rho = rho, act_FDR = FDR_spearman) |>
  filter(!is.na(ensembl)) |> distinct(ensembl, .keep_all = TRUE)
tab <- left_join(tab, act_tf, by = "ensembl")

## ---------------------------------------------------------------- Cox
cox <- read.csv(COX_FILE) |> filter(model == "multivariable") |> mutate(variable = fix_cl(variable))
check_ids(cox$variable, "cox_resultados")
for (nm in names(COX_COLS)) {
  sel <- cox |> filter(cohort == COX_COLS[[nm]][1], endpoint == COX_COLS[[nm]][2])
  x <- sel |> filter(!variable %in% CLUSTER_ROWS, variable != "cluster_score") |>
    transmute(ensembl = to_ens(variable), HR = HR, CI_low = CI_low, CI_high = CI_high, FDR = FDR) |>
    filter(!is.na(ensembl)) |> distinct(ensembl, .keep_all = TRUE)
  names(x)[-1] <- paste0(nm, "|", names(x)[-1])
  tab <- left_join(tab, x, by = "ensembl")
}

## ---------------------------------------------------------------- MRA (opcional)
read_mra <- function(f) {
  if (grepl("\\.rdata$|\\.rda$", tolower(f))) {
    suppressPackageStartupMessages(library(RTN))
    env <- new.env(); load(pe(f), envir = env)
    obj <- Filter(function(o) is(o, "TNA"), mget(ls(env), envir = env))
    if (!length(obj)) stop("Nenhum objeto TNA em ", f)
    message("MRA a partir do objeto: ", names(obj)[1])
    m <- tna.get(obj[[1]], what = "mra", ntop = -1)
    m$id <- rownames(m)
    m
  } else {
    m <- read.csv(pe(f), check.names = FALSE)
    if (!"id" %in% names(m)) m$id <- if (!is.null(m[[1]])) as.character(m[[1]]) else NA
    m
  }
}
if (!is.null(MRA_FILE) && file.exists(pe(MRA_FILE))) {
  mra <- read_mra(MRA_FILE)
  # id (rownames) e Regulon podem ser ENSG ou símbolo: converte os dois e fica com o que casar
  ens <- to_ens(mra$id)
  if ("Regulon" %in% names(mra)) ens <- ifelse(is.na(ens), to_ens(mra$Regulon), ens)
  look <- data.frame(ensembl = ens, MRA_padj = mra$Adjusted.Pvalue) |>
    filter(!is.na(ensembl)) |> distinct(ensembl, .keep_all = TRUE)
  tab <- left_join(tab, look, by = "ensembl")
  if (any(is.na(tab$MRA_padj)))
    message("[aviso] TFs da assinatura fora da tabela da MRA: ", paste(tab$TF[is.na(tab$MRA_padj)], collapse = ", "))
}

miss <- tab |> filter(is.na(act_rho) | is.na(.data[[paste0(names(COX_COLS)[1], "|HR")]]))
if (nrow(miss)) message("[aviso] TFs sem atividade/Cox (confira MR_FILES nos scripts de origem): ",
                        paste(miss$TF, collapse = ", "))

tab <- tab |> mutate(cluster = factor(cluster, levels = CLUSTER_ORDER)) |> arrange(cluster, TF)
write.csv(tab, file.path(OUT_DIR, "Table2_signature_TFs.csv"), row.names = FALSE)


## ---------------------------------------------------------------- tabela formatada
fmt <- tab |> transmute(
  TF, Cluster = sub(" cluster", "", as.character(cluster)),
  `Tumour vs normal, log2FC (adj. p)` = sprintf("%.2f (%s)", TvN_log2FC, fp(TvN_padj)),
  `M1 vs M0, δ (adj. p)` = sprintf("%.2f (%s)", M1vM0_delta, fp(M1vM0_padj)),
  `ICGC direction` = ifelse(is.na(ICGC_concordant), "NA", ifelse(ICGC_concordant, "Concordant", "Discordant")),
  !!paste0("Activity vs ", tolower(GRADE_FEATURE), ", ρ (adj. p)") := sprintf("%.2f (%s)", act_rho, fp(act_FDR)))
for (nm in names(COX_COLS)) {
  fmt[[paste0("HR ", nm, " (adj. p)")]] <- sprintf("%.2f (%.2f–%.2f; %s)",
                                                   tab[[paste0(nm, "|HR")]], tab[[paste0(nm, "|CI_low")]], tab[[paste0(nm, "|CI_high")]],
                                                   fp(tab[[paste0(nm, "|FDR")]]))
}
if ("MRA_padj" %in% names(tab)) fmt <- fmt |> mutate(`MRA adj. p` = fp(tab$MRA_padj), .after = Cluster)

# linhas finais: escores dos clusters (p nominal)
for (cv in names(CLUSTER_ROWS_LABEL)) {
  clu_act <- act |> filter(variable == cv)
  clu_row <- setNames(as.list(rep("", ncol(fmt))), names(fmt))
  clu_row$TF <- paste(CLUSTER_ROWS_LABEL[[cv]], "score"); clu_row$Cluster <- sub(" cluster", "", CLUSTER_ROWS_LABEL[[cv]])
  if (nrow(clu_act)) clu_row[[grep("^Activity", names(fmt), value = TRUE)]] <-
    sprintf("%.2f (%s)", clu_act$rho, fp(clu_act$p_spearman))
  for (nm in names(COX_COLS)) {
    cc <- cox |> filter(variable == cv, cohort == COX_COLS[[nm]][1], endpoint == COX_COLS[[nm]][2])
    if (nrow(cc)) clu_row[[paste0("HR ", nm, " (adj. p)")]] <- sprintf("%.2f (%.2f–%.2f; %s)",
                                                                       cc$HR, cc$CI_low, cc$CI_high, fp(cc$p))
  }
  if (nrow(clu_act) || any(cox$variable == cv)) fmt <- bind_rows(fmt, as.data.frame(clu_row, check.names = FALSE))
  else message("[aviso] escore ", cv, " ausente nas saídas (rode de novo o Cox e a atividade)")
}
write.csv(fmt, file.path(OUT_DIR, "Table2_signature_TFs_formatted.csv"), row.names = FALSE)

if (requireNamespace("knitr", quietly = TRUE)) {
  # cabeçalho em dois níveis: Expressão | Atividade do regulon | Sobrevivência
  n_id  <- which(names(fmt) == "Cluster") + ("MRA adj. p" %in% names(fmt))
  hdr <- c(" " = n_id, "Expression (TCGA-KIRC)" = 3, "Regulon activity" = 1, "Cox, HR per SD (95% CI; adj. p)" = length(COX_COLS))
  short <- names(fmt)
  short <- sub("^Tumour vs normal, log2FC \\(adj. p\\)$", "T vs N, log2FC (adj. p)", short)
  short <- sub("^ICGC direction$", "ICGC", short)
  short <- sub("^HR (.*) \\(adj. p\\)$", "\\1", short)
  tex <- knitr::kable(setNames(fmt, short), format = "latex", booktabs = TRUE, escape = TRUE,
                      caption = "Transcription factors of the metastasis-associated signature.")
  if (requireNamespace("kableExtra", quietly = TRUE)) {
    tex <- kableExtra::add_header_above(tex, hdr)
    tex <- kableExtra::kable_styling(tex, latex_options = "scale_down", font_size = 7)
    tex <- kableExtra::landscape(tex)
  }
  for (r in list(c("\u03b4", "$\\delta$"), c("\u03c1", "$\\rho$"), c("\u2013", "--")))
    tex <- gsub(enc2utf8(r[1]), r[2], enc2utf8(as.character(tex)), fixed = TRUE, useBytes = TRUE)
  writeLines(tex, file.path(OUT_DIR, "Table2_signature_TFs.tex"))
}
print(fmt, right = FALSE)
cat("\nSaídas em", OUT_DIR, "\n")