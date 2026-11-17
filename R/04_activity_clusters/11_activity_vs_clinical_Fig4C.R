### Atividade dos regulons (RTN) vs grau, pT e estágio ###
# Revisão CSBJ-D-26-01039, R1.M2: "consider examining MR activity in relation to grade, stage,
# or subsequent progression" (progressão: ver survival_cox_cdr.R)
#
# Para o escore do cluster imune e para cada TF das assinaturas:
#   - Spearman com a variável ordinal (grau 1-4, pT 1-4, estágio I-IV): tendência monotônica
#   - Kruskal-Wallis entre os grupos agrupados (G1-G2 / G3 / G4; T1-T2 / T3-T4; I-II / III / IV)
#   - em todos os pacientes e só nos M0 (mostra que a associação não vem só dos M1)
# A atividade é lida exatamente como em survival_cox_cdr.R (mesmo objeto, mesmo escore).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## ---------------------------------------------------------------- CONFIG (igual ao script de sobrevivência)
RTN_FILE     <- "0_data/rtni_TF.RData"
RTN_OBJECT   <- "rtni_TF"
ACTIVITY_RDS <- NULL                  # opcional: matriz amostras x regulons (dES) já salva
META_FILE    <- "0_data/metadata_KIRC_full.csv"
META_ID_COL  <- "sample"
M_COL        <- "ajcc_pathologic_m"
GRADE_COL    <- "neoplasm.histologic.grade"
PT_COL       <- "ajcc_pathologic_T"
STAGE_COL    <- "ajcc_pathologic_stage"   # se não existir, o estágio é pulado
OUT_DIR      <- "results_activity_clinical"

# Clusters de regulons por SÍMBOLO (convertidos para os nomes do RTN automaticamente).
# O primeiro é o teste primário; os escores dos clusters ficam fora da correção de BH dos TFs.
CLUSTERS <- list(
  "Cluster_1" = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A"),
  "Cluster_2" = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4")
)
CLUSTER_COLORS <- c(cluster_2 = "#1a9850", cluster_1 = "#542788")
MR_FILES <- c(signature = "results_ml_revision/importancia_features_modelo_final_class_weight.csv")   # assinatura única do artigo (14 TFs)
SYMBOL_FILE <- "0_data/ensembl_to_symbol.csv"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- atividade (mesma lógica do script de sobrevivência)
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
CLUSTER_VARS  <- gsub("[^a-z0-9]+", "_", tolower(names(CLUSTERS)))          # cluster 2, cluster 1
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
    message("Rodando tni.gsea2 para ", length(tfs_needed), " regulons...")
    rtni <- tni.gsea2(rtni, regulatoryElements = tfs_needed)
    act <- tni.get(rtni, what = "regulonActivity")
  }
  act$differential
}
dES <- get_activity()
# converte clusters e assinatura para os nomes das colunas da atividade (ex.: ENSG -> símbolo)
CLUSTER_REG <- lapply(CLUSTER_REG, to_reg, regs = colnames(dES))
MR_SETS     <- lapply(MR_SETS, to_reg, regs = colnames(dES))
tfs_needed  <- unique(c(unlist(CLUSTER_REG), unlist(MR_SETS)))
miss <- setdiff(tfs_needed, colnames(dES))
if (length(miss)) stop("Regulons sem atividade: ", paste(miss, collapse = ", "))

symbols <- if (file.exists(SYMBOL_FILE)) {
  s <- read.csv(SYMBOL_FILE); setNames(s$symbol, sub("\\..*", "", s$ensembl))
} else character(0)
lab <- function(g) ifelse(sub("\\..*", "", g) %in% names(symbols), symbols[sub("\\..*", "", g)], g)

act <- as.data.frame(dES[, tfs_needed, drop = FALSE])
act$sample <- rownames(dES)
act <- act |> filter(substr(sample, 14, 15) == "01") |>
  mutate(patient = substr(sample, 1, 12)) |> distinct(patient, .keep_all = TRUE)
# escore de cada cluster = média dos dES dos TFs do cluster
for (cv in CLUSTER_VARS) act[[cv]] <- rowMeans(act[, CLUSTER_REG[[cv]], drop = FALSE])

## ---------------------------------------------------------------- clínico
meta <- read.csv(META_FILE, check.names = FALSE)
find_col <- function(name) {
  if (name %in% names(meta)) return(name)
  norm <- function(s) tolower(gsub("[._]", "", s))
  hit <- names(meta)[norm(names(meta)) == norm(name)]
  if (length(hit)) hit[1] else NA_character_
}
gc <- find_col(GRADE_COL); tc <- find_col(PT_COL); sc <- find_col(STAGE_COL); mc <- find_col(M_COL)
if (is.na(gc) || is.na(tc) || is.na(mc)) stop("Confira GRADE_COL, PT_COL e M_COL. Colunas: ", paste(names(meta), collapse = ", "))
if (is.na(sc)) message("[aviso] coluna de estágio não encontrada; estágio pulado")

clin <- data.frame(
  patient = substr(meta[[META_ID_COL]], 1, 12),
  M = meta[[mc]],
  grade_ord = as.numeric(sub(".*G([1-4]).*", "\\1", ifelse(grepl("G[1-4]", meta[[gc]]), meta[[gc]], NA))),
  pT_ord    = as.numeric(sub(".*T([1-4]).*", "\\1", ifelse(grepl("T[1-4]", meta[[tc]]), meta[[tc]], NA)))
)
if (!is.na(sc)) {
  st <- toupper(meta[[sc]])
  clin$stage_ord <- case_when(grepl("STAGE IV", st) ~ 4, grepl("STAGE III", st) ~ 3,
                              grepl("STAGE II", st) ~ 2, grepl("STAGE I", st) ~ 1)
}
clin <- clin |> distinct(patient, .keep_all = TRUE) |>
  mutate(Grade = factor(c("G1-G2", "G1-G2", "G3", "G4")[grade_ord], levels = c("G1-G2", "G3", "G4")),
         pT = factor(ifelse(pT_ord <= 2, "T1-T2", "T3-T4"), levels = c("T1-T2", "T3-T4")))
if ("stage_ord" %in% names(clin))
  clin$Stage <- factor(c("I-II", "I-II", "III", "IV")[clin$stage_ord], levels = c("I-II", "III", "IV"))

df <- inner_join(act, clin, by = "patient")
cohorts <- list(all = df, M0_only = filter(df, M == "M0"))
write.csv(df, file.path(OUT_DIR, "dados_por_paciente.csv"), row.names = FALSE)

features <- list(Grade = "grade_ord", pT = "pT_ord")
if ("Stage" %in% names(df)) features$Stage <- "stage_ord"

## ---------------------------------------------------------------- testes
vars <- c(CLUSTER_VARS, unique(unlist(MR_SETS)))
rows <- list()
for (coh in names(cohorts)) for (f in names(features)) for (v in vars) {
  d <- cohorts[[coh]]
  if (coh == "M0_only" && f == "Stage") d <- filter(d, stage_ord < 4)   # estágio IV nos M0 é raro (T4)
  d <- d[!is.na(d[[features[[f]]]]) & !is.na(d[[f]]), ]
  if (length(unique(d[[f]])) < 2) next
  sp <- suppressWarnings(cor.test(d[[v]], d[[features[[f]]]], method = "spearman"))
  kw <- kruskal.test(d[[v]] ~ d[[f]])
  rows[[length(rows) + 1]] <- data.frame(
    cohort = coh, feature = f, variable = v,
    label = ifelse(v %in% CLUSTER_VARS, CLUSTER_LABEL[v], lab(v)),
    n = nrow(d), rho = unname(sp$estimate), p_spearman = sp$p.value, p_kruskal = kw$p.value)
}
res <- bind_rows(rows)
tf <- !res$variable %in% CLUSTER_VARS     # FDR entre os TFs; os escores dos clusters ficam fora
res$FDR_spearman <- NA_real_
res$FDR_spearman[tf] <- ave(res$p_spearman[tf], res$cohort[tf], res$feature[tf],
                            FUN = function(p) p.adjust(p, "BH"))
write.csv(res, file.path(OUT_DIR, "atividade_vs_clinico.csv"), row.names = FALSE)

## ---------------------------------------------------------------- figuras
SEQ <- c("#86b6ef", "#2a78d6", "#1c5cab"); INK2 <- "#52514e"
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

plot_cohort <- function(coh, title, cv) {
  d <- cohorts[[coh]]
  long <- bind_rows(lapply(names(features), function(f) {
    dd <- d[!is.na(d[[f]]), ]
    if (coh == "M0_only" && f == "Stage") dd <- filter(dd, stage_ord < 4) |> mutate(Stage = droplevels(Stage))
    data.frame(feature = f, group = as.character(dd[[f]]), score = dd[[cv]])
  }))
  long$feature <- factor(long$feature, levels = names(features))
  lev <- unlist(lapply(names(features), function(f) levels(droplevels(factor(d[[f]][!is.na(d[[f]])])))))
  long$group <- factor(long$group, levels = unique(lev))
  long <- long |> group_by(feature) |>
    mutate(shade = SEQ[pmin(as.integer(droplevels(group)), 3)]) |> ungroup()
  ann <- res |> filter(cohort == coh, variable == cv) |>
    mutate(feature = factor(feature, levels = names(features)),
           txt = sprintf("rho = %.2f, %s", rho, fmt_p(p_spearman)))
  ns <- long |> count(feature, group)
  ymax <- max(long$score); ymin <- min(long$score); span <- ymax - ymin
  ggplot(long, aes(group, score)) +
    geom_jitter(width = 0.18, size = 0.7, alpha = 0.35, color = "#8a8983") +
    geom_boxplot(aes(fill = shade), width = 0.55, outlier.shape = NA, alpha = 0.85,
                 color = "#0b0b0b", linewidth = 0.35) +
    scale_fill_identity() +
    geom_text(data = ns, aes(group, ymin - 0.08 * span, label = paste0("n = ", n)),
              size = 2.6, color = INK2) +
    geom_text(data = ann, aes(x = Inf, y = ymax + 0.12 * span, label = txt),
              hjust = 1.05, size = 2.8, color = INK2, inherit.aes = FALSE) +
    facet_wrap(~feature, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.12))) +
    labs(x = NULL, y = paste(CLUSTER_LABEL[[cv]], "activity (mean dES)"), title = title) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(color = INK2), strip.text = element_text(face = "bold"))
}
for (cv in CLUSTER_VARS) for (coh in names(cohorts)) {
  title <- paste0(CLUSTER_LABEL[[cv]], ", ", if (coh == "all") "all patients" else "M0 patients (stage IV excluded)")
  g <- plot_cohort(coh, title, cv)
  w <- 2.3 * length(features) + 0.6
  ggsave(file.path(OUT_DIR, paste0("atividade_", cv, "_vs_clinico_", coh, ".pdf")), g, width = w, height = 3.2)
  ggsave(file.path(OUT_DIR, paste0("atividade_", cv, "_vs_clinico_", coh, ".png")), g, width = w, height = 3.2, dpi = 300)
}

## ---------------------------------------------------------------- escores dos clusters: M1 vs M0 (Fig 4)
cliff <- function(x, y) (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
dm <- filter(df, M %in% c("M0", "M1"))
m_tab <- bind_rows(lapply(CLUSTER_VARS, function(cv) {
  x <- dm[[cv]][dm$M == "M1"]; y <- dm[[cv]][dm$M == "M0"]
  data.frame(cluster = CLUSTER_LABEL[[cv]], n_M1 = length(x), n_M0 = length(y),
             median_M1 = median(x), median_M0 = median(y), cliffs_delta = cliff(x, y),
             p_mannwhitney = wilcox.test(x, y)$p.value)
}))
write.csv(m_tab, file.path(OUT_DIR, "escores_clusters_M1_vs_M0.csv"), row.names = FALSE)

long_m <- bind_rows(lapply(CLUSTER_VARS, function(cv)
  data.frame(cluster = CLUSTER_LABEL[[cv]], M = dm$M, score = dm[[cv]])))
long_m$cluster <- factor(long_m$cluster, levels = unname(CLUSTER_LABEL))
ann_m <- m_tab |> mutate(cluster = factor(cluster, levels = unname(CLUSTER_LABEL)),
                         txt = sprintf("δ = %.2f\np = %s", cliffs_delta,
                                       ifelse(p_mannwhitney < 0.001, formatC(p_mannwhitney, format = "e", digits = 1),
                                              sprintf("%.3f", p_mannwhitney))))
CL_COL <- c(cluster_2 = "#1a9850", cluster_1 = "#542788")
fill_map <- setNames(CL_COL[CLUSTER_VARS], unname(CLUSTER_LABEL))
gm <- ggplot(long_m, aes(M, score, fill = cluster)) +
  geom_jitter(width = 0.18, size = 0.6, alpha = 0.3, color = "#8a8983") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8, color = "#0b0b0b", linewidth = 0.35) +
  geom_text(data = ann_m, aes(x = 1.5, y = Inf, label = txt), vjust = 1.3, size = 2.8, color = INK2,
            inherit.aes = FALSE) +
  scale_fill_manual(values = fill_map, guide = "none") +
  scale_x_discrete(labels = c(M0 = sprintf("M0\n(n = %d)", sum(dm$M == "M0")),
                              M1 = sprintf("M1\n(n = %d)", sum(dm$M == "M1")))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  facet_wrap(~cluster, scales = "free_y") +
  labs(x = NULL, y = "Cluster activity score (mean dES)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(color = INK2), strip.text = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "escores_clusters_M1_vs_M0.pdf"), gm, width = 5, height = 3.2)
ggsave(file.path(OUT_DIR, "escores_clusters_M1_vs_M0.png"), gm, width = 5, height = 3.2, dpi = 300)
ggsave(file.path(OUT_DIR, "escores_clusters_M1_vs_M0.tiff"), gm, width = 5, height = 3.2, dpi = 300, compression = "lzw")
cat("\n=== Escores dos clusters: M1 vs M0 (Mann-Whitney) ===\n")
print(m_tab |> mutate(across(where(is.numeric), ~ signif(.x, 3))))

## ---------------------------------------------------------------- resumo
cat("\n=== Escores dos clusters vs variáveis clínicas (Spearman com a variável ordinal) ===\n")
print(res |> filter(variable %in% CLUSTER_VARS) |>
        transmute(cluster = label, cohort, feature, n, rho = round(rho, 2), p_spearman = signif(p_spearman, 3),
                  p_kruskal = signif(p_kruskal, 3)))
cat("\nSaídas em", OUT_DIR, "\n")