### Tabela 1 - Características clínicas por status metastático ###
# Revisão CSBJ-D-26-01039: R2.M3, R2.adic ("table summarizing the clinical
# characteristics of both cohorts"), R1.M2, R3.2
#
# Ajuste os nomes das colunas no bloco CONFIG. Colunas que não existirem
# no arquivo são apenas avisadas e puladas.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(gtsummary)
  library(flextable)
})

## ---------------------------------------------------------------- CONFIG
OUT_DIR <- "results_table1"
dir.create(OUT_DIR, showWarnings = FALSE)

# Valores que significam "não informado" em qualquer coluna categórica
UNKNOWN <- c("", "NA", "Not Reported", "not reported", "[Not Available]",
             "[Not Evaluated]", "[Unknown]", "Unknown", "[Discrepancy]",
             "GX", "TX", "NX", "MX")

TCGA_CFG <- list(
  file      = "0_data/metadata_KIRC_full.csv",
  id_col    = "sample",                 # coluna com o barcode (vroom/readr nomeia o índice assim)
  m_col     = "ajcc_pathologic_m",
  is_tcga   = TRUE,                   # aplica filtro de tumor primário e 1 amostra por paciente
  vars = c(                           # nome no arquivo = rótulo na tabela
    age_at_index               = "Age at diagnosis (years)",
    gender                     = "Sex",
    race                       = "Race",
    neoplasm.histologic.grade  = "Histological grade",
    ajcc_pathologic_t          = "Pathological T",
    ajcc_pathologic_n          = "Pathological N",
    ajcc_pathologic_stage      = "AJCC stage",
    longest_dimension              = "Tumor size, longest dimension (cm)",     # ajuste ou remova
    VHL_Status                 = "VHL status",          # ajuste ou remova
    vital_status               = "Vital status"
  ),
  continuous = c("age_at_index", "longest_dimension"),
  no_test    = c("ajcc_pathologic_stage"),   # estágio IV = M1 por definição: descrever sem p-valor
  label      = "TCGA-KIRC"
)

ICGC_CFG <- list(
  file      = "0_data/meta_icgc_data.csv",
  id_col    = "patient",
  m_col     = "pathologic_m",
  is_tcga   = FALSE,
  vars = c(
    age = "Age at diagnosis (years)",
    gender              = "Sex",
    tumour_grade           = "Histological grade",
    pathologic_t           = "Pathological T",
    pathologic_n           = "Pathological N",
    tumour_stage           = "Stage",
    status     = "Vital status"
  ),
  continuous = c("age"),
  no_test    = c("tumour_stage"),
  label      = "ICGC-RECA"
)

## ---------------------------------------------------------------- funções
load_cohort <- function(cfg) {
  meta <- read_csv(cfg$file, show_col_types = FALSE)
  flow <- c(inicial = nrow(meta))
  
  stopifnot(cfg$id_col %in% names(meta), cfg$m_col %in% names(meta))
  meta <- meta %>% rename(id = all_of(cfg$id_col), M = all_of(cfg$m_col))
  
  meta <- meta %>% filter(M %in% c("M0", "M1"))
  flow["M0_M1"] <- nrow(meta)
  
  if (cfg$is_tcga) {
    meta <- meta %>% filter(str_sub(id, 14, 15) == "01")          # tumor primário
    flow["tumor_primario"] <- nrow(meta)
    meta <- meta %>% distinct(paciente = str_sub(id, 1, 12), .keep_all = TRUE)
    flow["1_por_paciente"] <- nrow(meta)
  }
  flow["M0"] <- sum(meta$M == "M0"); flow["M1"] <- sum(meta$M == "M1")
  message(cfg$label, " | fluxo: ", paste(names(flow), flow, sep = "=", collapse = ", "))
  list(meta = meta, flow = flow)
}

prepare <- function(meta, cfg) {
  present <- intersect(names(cfg$vars), names(meta))
  missing <- setdiff(names(cfg$vars), names(meta))
  if (length(missing))
    warning(cfg$label, ": colunas não encontradas e puladas: ",
            paste(missing, collapse = ", "),
            "\n  Colunas disponíveis: ", paste(names(meta), collapse = ", "))
  
  df <- meta %>% select(M, all_of(present))
  
  cont <- intersect(cfg$continuous, present)
  cat_vars <- setdiff(present, cont)
  
  df <- df %>%
    mutate(across(all_of(cont), ~ suppressWarnings(as.numeric(.x))),
           across(all_of(cat_vars), ~ {
             x <- as.character(.x)
             x[x %in% UNKNOWN] <- NA          # "não informado" vira NA e sai dos testes
             factor(x)
           }),
           M = factor(M, levels = c("M0", "M1")))
  
  # TCGA às vezes traz idade em dias negativos (days_to_birth)
  for (v in cont) if (median(df[[v]], na.rm = TRUE) < 0) df[[v]] <- abs(df[[v]]) / 365.25
  
  list(df = df, present = present, cont = cont, cat_vars = cat_vars)
}

make_table1 <- function(cfg) {
  ld <- load_cohort(cfg)
  pr <- prepare(ld$meta, cfg)
  labels <- as.list(cfg$vars[pr$present])
  
  tbl <- pr$df %>%
    tbl_summary(
      by = M,
      label = labels,
      type = list(any_of(pr$cont) ~ "continuous"),
      statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                       all_categorical() ~ "{n} ({p}%)"),
      missing = "ifany",
      missing_text = "Not reported"
    ) %>%
    add_overall() %>%
    add_p(
      include = -any_of(cfg$no_test),
      test = list(all_continuous() ~ "wilcox.test"),
      # categóricas: gtsummary usa chisq.test, e fisher.test quando há
      # células esperadas < 5; simulação evita falha em tabelas grandes
      test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE, B = 1e5),
      pvalue_fun = ~ style_pvalue(.x, digits = 3)
    ) %>%
    add_q(method = "fdr") %>%
    bold_labels() %>%
    modify_caption(paste0("Table 1. Clinical characteristics by metastatic status at diagnosis, ",
                          cfg$label, ". Stage is shown without a test because M1 implies stage IV."))
  
  ft <- as_flex_table(tbl)
  save_as_docx(ft, path = file.path(OUT_DIR, paste0("Table1_", cfg$label, ".docx")))
  write_csv(as_tibble(tbl, col_labels = TRUE),
            file.path(OUT_DIR, paste0("Table1_", cfg$label, ".csv")))
  list(tbl = tbl, flow = ld$flow)
}

## ---------------------------------------------------------------- execução
res <- list()
res$tcga <- make_table1(TCGA_CFG)
if (file.exists(ICGC_CFG$file)) res$icgc <- make_table1(ICGC_CFG)

# Tabela combinada lado a lado (TCGA | ICGC), útil para o texto principal
if (!is.null(res$icgc)) {
  merged <- tbl_merge(list(res$tcga$tbl, res$icgc$tbl),
                      tab_spanner = c("**TCGA-KIRC**", "**ICGC-RECA**"))
  save_as_docx(as_flex_table(merged), path = file.path(OUT_DIR, "Table1_TCGA_ICGC.docx"))
}

# Fluxo de pacientes (para o fluxograma, R3.2)
flow_df <- bind_rows(lapply(res, function(r) as.data.frame(t(r$flow))), .id = "coorte")
write_csv(flow_df, file.path(OUT_DIR, "fluxo_pacientes.csv"))
print(flow_df)
res$tcga$tbl
res$icgc$tbl
