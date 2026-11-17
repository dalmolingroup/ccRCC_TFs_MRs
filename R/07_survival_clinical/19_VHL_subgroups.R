### Subgrupo VHL selvagem — R1.M6 ###
# "The findings may reflect VHL/HIF-driven biology [...] Please test whether your findings pan out
#  in VHL-wild-type ccRCC."
#
# Status do VHL montado DIRETO do MAF (não da coluna VHL_Status do metadata):
#   Mutado      = mutação inativante no VHL (Missense, Nonsense, Frame_Shift, Splice_Site, In_Frame,
#                 Translation_Start_Site, Nonstop)
#   Selvagem    = paciente presente no MAF (foi sequenciado) sem mutação inativante no VHL
#   Desconhecido= paciente ausente do MAF (sem sequenciamento) -> fora da análise do subgrupo
# Limitação a declarar: a metilação do promotor do VHL (~5-10%) não foi considerada.
#
# Análises (mesma atividade do RTN e mesmos desfechos de survival_cox_cdr.R):
#   1. Contagens por status do VHL e M
#   2. Atividade de cada cluster: VHL mutado vs selvagem (o cluster é só reflexo do VHL?)
#   3. Nos VHL selvagem: atividade M1 vs M0 e vs grau/pT
#   4. Cox (OS, PFI, DSS) por subgrupo + teste de interação cluster x VHL em todos os pacientes

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(survival)
  library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG
RTN_FILE     <- "0_data/rtni_TF.RData"
RTN_OBJECT   <- "rtni_TF"
ACTIVITY_RDS <- NULL
CDR_FILE     <- "0_data/TCGA-CDR-SupplementalTableS1.xlsx"
META_FILE    <- "0_data/metadata_KIRC_full.csv"       # coorte do ML: define pacientes, M e pT
META_ID_COL  <- "sample"
M_COL        <- "ajcc_pathologic_m"
PT_COL       <- "ajcc_pathologic_t"
MAF_FILE     <- "0_data/muta.rds"                      # .rds, .RData (objeto MAF_OBJECT) ou .maf/.maf.gz/.tsv
MAF_OBJECT   <- "muta"
OUT_DIR      <- "results_vhl"

# Clusters por SÍMBOLO (os mesmos de survival_cox_cdr.R); convertidos para os nomes do RTN
CLUSTERS <- list(
  "Differentiation-associated cluster" = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  "Proliferative cluster"              = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A")
)
SYMBOL_FILE  <- "0_data/ensembl_to_symbol.csv"
ENDPOINTS    <- c("OS", "PFI", "DSS")
INACTIVATING <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                  "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site",
                  "Nonstop_Mutation")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- 1. status do VHL a partir do MAF
read_maf <- function(f) {
  if (grepl("\\.rds$", f, ignore.case = TRUE)) return(readRDS(f))
  if (grepl("\\.rdata$|\\.rda$", f, ignore.case = TRUE)) {
    e <- new.env(); load(f, envir = e); return(get(MAF_OBJECT, envir = e))
  }
  as.data.frame(data.table::fread(f, sep = "\t", quote = "", skip = "Hugo_Symbol"))
}
maf <- as.data.frame(read_maf(MAF_FILE))
maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)
sequenced <- unique(maf$patient)
vhl_rows <- maf[maf$Hugo_Symbol == "VHL", c("patient", "Variant_Classification", "HGVSp_Short")]
cat("Classes de variante no VHL (MAF):\n"); print(table(vhl_rows$Variant_Classification))
vhl_mut <- unique(vhl_rows$patient[vhl_rows$Variant_Classification %in% INACTIVATING])
write.csv(vhl_rows, file.path(OUT_DIR, "VHL_variantes_MAF.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 2. atividade (igual aos outros scripts)
id_key <- function(x) toupper(sub("\\..*", "", trimws(as.character(x))))
to_reg <- function(ids, regs) {                      # símbolo/ENSG -> nome usado em regs
  vals <- unname(regs); nm <- names(regs)
  s <- if (file.exists(SYMBOL_FILE)) read.csv(SYMBOL_FILE) else data.frame(ensembl = character(), symbol = character())
  e2s <- setNames(toupper(s$symbol), id_key(s$ensembl)); s2e <- setNames(id_key(s$ensembl), toupper(s$symbol))
  out <- vapply(ids, function(id) {
    k <- id_key(id); cand <- unique(na.omit(c(k, unname(e2s[k]), unname(s2e[k]),
                                              unname(REG_V2N[k]), unname(REG_N2V[k]))))
    for (pool in list(vals, nm)) {
      if (is.null(pool)) next
      hit <- match(cand, id_key(pool)); hit <- hit[!is.na(hit)]
      if (length(hit)) return(vals[hit[1]])
    }
    NA_character_
  }, character(1))
  if (any(is.na(out))) stop("Não encontrados: ", paste(ids[is.na(out)], collapse = ", "))
  unname(out)
}
if (is.null(ACTIVITY_RDS)) {
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(RTN_FILE, envir = env); rtni <- get(RTN_OBJECT, envir = env)
  regs_all <- tni.get(rtni, what = "regulatoryElements")
} else regs_all <- colnames(readRDS(ACTIVITY_RDS))
REG_V2N <- if (!is.null(names(regs_all))) setNames(id_key(names(regs_all)), id_key(regs_all)) else character(0)
REG_N2V <- if (!is.null(names(regs_all))) setNames(id_key(regs_all), id_key(names(regs_all))) else character(0)
CLUSTER_VARS  <- gsub("[^a-z]+", "_", tolower(names(CLUSTERS)))
CLUSTER_LABEL <- setNames(names(CLUSTERS), CLUSTER_VARS)
CLUSTER_REG   <- setNames(lapply(CLUSTERS, to_reg, regs = regs_all), CLUSTER_VARS)

get_activity <- function() {
  if (!is.null(ACTIVITY_RDS)) return(readRDS(ACTIVITY_RDS))
  act <- tryCatch(tni.get(rtni, what = "regulonActivity"), error = function(e) NULL)
  have <- if (!is.null(act$differential)) colnames(act$differential) else character(0)
  ok <- length(have) > 0 && !inherits(try(to_reg(unlist(CLUSTER_REG), have), silent = TRUE), "try-error")
  if (!ok) {
    rtni <- tni.gsea2(rtni, regulatoryElements = unlist(CLUSTER_REG))
    act <- tni.get(rtni, what = "regulonActivity")
  }
  act$differential
}
dES <- get_activity()
CLUSTER_REG <- lapply(CLUSTER_REG, to_reg, regs = colnames(dES))     # nomes das colunas da atividade
act <- data.frame(sample = rownames(dES))
for (cv in CLUSTER_VARS) act[[cv]] <- rowMeans(dES[, CLUSTER_REG[[cv]], drop = FALSE])
act <- act |> filter(substr(sample, 14, 15) == "01") |>
  mutate(patient = substr(sample, 1, 12)) |> distinct(patient, .keep_all = TRUE)

## ---------------------------------------------------------------- 3. clínico
meta <- read.csv(META_FILE, check.names = FALSE)
find_col <- function(d, name) {
  if (name %in% names(d)) return(name)
  hit <- names(d)[tolower(gsub("[._]", "", names(d))) == tolower(gsub("[._]", "", name))]
  if (length(hit)) hit[1] else stop("Coluna não encontrada: ", name)
}
meta <- data.frame(patient = substr(meta[[META_ID_COL]], 1, 12),
                   M = meta[[find_col(meta, M_COL)]],
                   pT_ord = suppressWarnings(as.numeric(sub(".*T([1-4]).*", "\\1", meta[[find_col(meta, PT_COL)]])))) |>
  distinct(patient, .keep_all = TRUE)

cdr <- read_excel(CDR_FILE, sheet = "TCGA-CDR", na = c("", "NA", "#N/A", "[Not Available]", "[Not Applicable]",
                                                       "[Unknown]", "[Discrepancy]", "[Not Evaluated]")) |>
  filter(type == "KIRC") |>
  transmute(patient = bcr_patient_barcode,
            age = as.numeric(age_at_initial_pathologic_diagnosis),
            sex = factor(tolower(gender), levels = c("female", "male")),
            stage = factor(case_when(grepl("Stage IV", ajcc_pathologic_tumor_stage) ~ "IV",
                                     grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "III",
                                     grepl("Stage I", ajcc_pathologic_tumor_stage) ~ "I-II"),
                           levels = c("I-II", "III", "IV")),
            grade_ord = suppressWarnings(as.numeric(sub("G", "", histological_grade))),
            grade = factor(case_when(histological_grade %in% c("G1", "G2") ~ "G1-G2",
                                     histological_grade == "G3" ~ "G3", histological_grade == "G4" ~ "G4"),
                           levels = c("G1-G2", "G3", "G4")),
            across(all_of(as.vector(rbind(ENDPOINTS, paste0(ENDPOINTS, ".time")))), as.numeric))

df <- meta |> inner_join(act, by = "patient") |> left_join(cdr, by = "patient") |>
  mutate(VHL = factor(case_when(patient %in% vhl_mut ~ "Mutated",
                                patient %in% sequenced ~ "Wild-type",
                                TRUE ~ NA_character_), levels = c("Wild-type", "Mutated")))

counts <- df |> count(VHL, M) |> pivot_wider(names_from = M, values_from = n, values_fill = 0)
cat("\nPacientes da coorte por status do VHL e M:\n"); print(counts)
write.csv(counts, file.path(OUT_DIR, "contagens_VHL_M.csv"), row.names = FALSE)
dv <- filter(df, !is.na(VHL))
wt <- filter(dv, VHL == "Wild-type")
write.csv(dv, file.path(OUT_DIR, "dados_VHL_por_paciente.csv"), row.names = FALSE)

ALL <- list()
for (clv in CLUSTER_VARS) {
  dv$score <- dv[[clv]]; wt <- filter(dv, VHL == "Wild-type")
  ## ---------------------------------------------------------------- 4. testes de atividade
  cliff <- function(x, y) unname(2 * suppressWarnings(wilcox.test(x, y, exact = FALSE))$statistic /
                                   (length(x) * length(y)) - 1)
  act_tests <- list()
  add_test <- function(label, x, y, gx, gy) {
    act_tests[[length(act_tests) + 1]] <<- data.frame(
      comparison = label, group1 = gx, n1 = length(x), median1 = median(x),
      group2 = gy, n2 = length(y), median2 = median(y), cliffs_delta = cliff(x, y),
      p = suppressWarnings(wilcox.test(x, y, exact = FALSE))$p.value)
  }
  add_test("VHL (all patients)", dv$score[dv$VHL == "Mutated"], dv$score[dv$VHL == "Wild-type"], "Mutated", "Wild-type")
  for (s in levels(dv$VHL)) {
    d <- filter(dv, VHL == s, M %in% c("M0", "M1"))
    add_test(paste0("M1 vs M0 (VHL ", s, ")"), d$score[d$M == "M1"], d$score[d$M == "M0"], "M1", "M0")
  }
  act_tab <- bind_rows(act_tests)
  cor_tab <- bind_rows(lapply(levels(dv$VHL), function(s) {
    d <- filter(dv, VHL == s)
    bind_rows(lapply(c(grade = "grade_ord", pT = "pT_ord"), function(v) {
      ok <- !is.na(d[[v]])
      ct <- suppressWarnings(cor.test(d$score[ok], d[[v]][ok], method = "spearman"))
      data.frame(VHL = s, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
    }), .id = "feature")
  }))
  write.csv(act_tab, file.path(OUT_DIR, paste0("atividade_testes_", clv, ".csv")), row.names = FALSE)
  write.csv(cor_tab, file.path(OUT_DIR, paste0("atividade_vs_grau_pT_por_VHL_", clv, ".csv")), row.names = FALSE)
  
  ## ---------------------------------------------------------------- 5. Cox por subgrupo + interação
  COVARS <- c("age", "sex", "stage", "grade")
  prep <- function(d, e) {
    tc <- paste0(e, ".time")
    d <- d[complete.cases(d[, c(tc, e, "score", COVARS)]) & d[[tc]] > 0, ]
    d$z <- (d$score - mean(dv$score)) / sd(dv$score)      # mesma escala (DP de todos) em todos os modelos
    d$stage <- droplevels(d$stage); d$grade <- droplevels(d$grade)
    d
  }
  fit_one <- function(d, e, subgroup, cohort) {
    d <- prep(d, e)
    if (sum(d[[e]]) < 10) return(NULL)                     # poucos eventos: não ajusta
    cv <- COVARS[vapply(COVARS, function(v) length(unique(d[[v]])) > 1, TRUE)]
    f <- as.formula(paste0("Surv(", e, ".time, ", e, ") ~ z + ", paste(cv, collapse = " + ")))
    s <- summary(coxph(f, data = d))
    data.frame(cohort = cohort, subgroup = subgroup, endpoint = e, n = s$n, events = s$nevent,
               HR = s$conf.int["z", 1], CI_low = s$conf.int["z", 3], CI_high = s$conf.int["z", 4],
               p = s$coefficients["z", 5])
  }
  cox_rows <- list(); inter_rows <- list()
  for (coh in c("all", "M0_only")) {
    base <- if (coh == "all") dv else filter(dv, M == "M0")
    for (e in ENDPOINTS) {
      for (s in levels(dv$VHL)) cox_rows[[length(cox_rows) + 1]] <- fit_one(filter(base, VHL == s), e, s, coh)
      d <- prep(base, e)
      cv <- COVARS[vapply(COVARS, function(v) length(unique(d[[v]])) > 1, TRUE)]
      rhs <- paste(c("z", "VHL", cv), collapse = " + ")
      f0 <- as.formula(paste0("Surv(", e, ".time, ", e, ") ~ ", rhs))
      f1 <- as.formula(paste0("Surv(", e, ".time, ", e, ") ~ ", rhs, " + z:VHL"))
      m0 <- coxph(f0, data = d); m1 <- coxph(f1, data = d)
      inter_rows[[length(inter_rows) + 1]] <- data.frame(
        cohort = coh, endpoint = e, n = m1$n, events = m1$nevent,
        p_interaction = anova(m0, m1)$`Pr(>|Chi|)`[2])
    }
  }
  cox_tab <- bind_rows(cox_rows); inter_tab <- bind_rows(inter_rows)
  write.csv(cox_tab, file.path(OUT_DIR, paste0("cox_por_VHL_", clv, ".csv")), row.names = FALSE)
  write.csv(inter_tab, file.path(OUT_DIR, paste0("cox_interacao_cluster_VHL_", clv, ".csv")), row.names = FALSE)
  
  ## ---------------------------------------------------------------- 6. figuras
  BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK2 <- "#52514e"
  fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  
  bx <- bind_rows(
    dv |> transmute(panel = "By VHL status (all patients)", group = as.character(VHL), score),
    wt |> filter(M %in% c("M0", "M1")) |> transmute(panel = "VHL wild-type: by metastatic status", group = M, score)
  ) |> mutate(group = factor(group, levels = c("Wild-type", "Mutated", "M0", "M1")))
  ann <- data.frame(panel = c("By VHL status (all patients)", "VHL wild-type: by metastatic status"),
                    label = c(sprintf("Cliff's delta = %.2f, %s", act_tab$cliffs_delta[1], fmt_p(act_tab$p[1])),
                              with(act_tab[act_tab$comparison == "M1 vs M0 (VHL Wild-type)", ],
                                   sprintf("Cliff's delta = %.2f, %s", cliffs_delta, fmt_p(p)))))
  ns <- bx |> count(panel, group)
  g1 <- ggplot(bx, aes(group, score)) +
    geom_jitter(width = 0.18, size = 0.7, alpha = 0.35, color = "#8a8983") +
    geom_boxplot(aes(fill = group), width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
    geom_text(data = ns, aes(group, -Inf, label = paste0("n = ", n)), vjust = -0.6, size = 2.6, color = INK2) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5,
              size = 2.8, color = INK2, inherit.aes = FALSE) +
    scale_fill_manual(values = c("Wild-type" = "#86b6ef", "Mutated" = "#1c5cab", M0 = "#2166ac", M1 = "#d73027"),
                      guide = "none") +
    facet_wrap(~panel, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(x = NULL, y = paste(CLUSTER_LABEL[[clv]], "activity (mean dES)")) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"), axis.text = element_text(color = INK2))
  ggsave(file.path(OUT_DIR, paste0("atividade_cluster_VHL_", clv, ".pdf")), g1, width = 6.4, height = 3.2)
  ggsave(file.path(OUT_DIR, paste0("atividade_cluster_VHL_", clv, ".png")), g1, width = 6.4, height = 3.2, dpi = 300)
  
  if (nrow(cox_tab)) {
    fp <- cox_tab |>
      left_join(inter_tab |> select(cohort, endpoint, p_interaction), by = c("cohort", "endpoint")) |>
      mutate(cohort = recode(cohort, all = "All patients", M0_only = "M0 only"),
             endpoint = factor(endpoint, levels = ENDPOINTS),
             subgroup = factor(subgroup, levels = c("Mutated", "Wild-type")),
             lab = sprintf("%s (n = %d; %d events)", subgroup, n, events))
    pint <- fp |> distinct(cohort, endpoint, p_interaction) |>
      mutate(txt = paste0("P interaction = ", sprintf("%.2f", p_interaction)))
    g2 <- ggplot(fp, aes(HR, subgroup, color = subgroup)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
      geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.6) +
      geom_point(size = 2.2) +
      geom_text(data = pint, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.2,
                size = 2.5, color = INK2, inherit.aes = FALSE) +
      scale_x_log10() +
      scale_color_manual(values = c("Wild-type" = BLUE, "Mutated" = ORANGE), guide = "none") +
      facet_grid(cohort ~ endpoint) +
      labs(x = paste("Hazard ratio per SD of", tolower(CLUSTER_LABEL[[clv]]), "activity (log scale)"), y = NULL,
           caption = "Multivariable Cox models adjusted for age, sex, stage and grade, fitted within each VHL subgroup.") +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
            strip.text = element_text(face = "bold"), plot.caption = element_text(color = INK2, hjust = 0))
    ggsave(file.path(OUT_DIR, paste0("forest_cox_por_VHL_", clv, ".pdf")), g2, width = 7, height = 3.6)
    ggsave(file.path(OUT_DIR, paste0("forest_cox_por_VHL_", clv, ".png")), g2, width = 7, height = 3.6, dpi = 300)
  }
  
  
  ALL[[clv]] <- list(act = cbind(cluster = CLUSTER_LABEL[[clv]], act_tab), cor = cbind(cluster = CLUSTER_LABEL[[clv]], cor_tab),
                     cox = if (nrow(cox_tab)) cbind(cluster = CLUSTER_LABEL[[clv]], cox_tab) else NULL,
                     inter = cbind(cluster = CLUSTER_LABEL[[clv]], inter_tab))
}
act_tab <- bind_rows(lapply(ALL, `[[`, "act")); cor_tab <- bind_rows(lapply(ALL, `[[`, "cor"))
cox_tab <- bind_rows(lapply(ALL, `[[`, "cox")); inter_tab <- bind_rows(lapply(ALL, `[[`, "inter"))
for (nm in c("act_tab", "cor_tab", "cox_tab", "inter_tab"))
  write.csv(get(nm), file.path(OUT_DIR, paste0(nm, "_dois_clusters.csv")), row.names = FALSE)

## ---------------------------------------------------------------- resumo
cat("\n=== Atividade dos clusters ===\n"); print(act_tab |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
cat("\n=== Atividade vs grau/pT por VHL ===\n"); print(cor_tab |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
cat("\n=== Cox por subgrupo (HR por DP) ===\n"); print(cox_tab |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
cat("\n=== Interação cluster x VHL ===\n"); print(inter_tab |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
cat("\nSaídas em", OUT_DIR, "\n")