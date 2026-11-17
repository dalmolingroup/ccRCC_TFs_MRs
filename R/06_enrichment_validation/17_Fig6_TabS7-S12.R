### Seção 3.4 — Tabelas suplementares e Figura 6 a partir das saídas do regulon_validation_enrichment.R ###
# Revisão CSBJ-D-26-01039 (R2.M6, R2.a8, R3.4, R3.7, R3.13)
#
# NÃO refaz o enriquecimento: só lê os CSVs já salvos em results_regulon_validation/ e organiza.
# (A única coisa lida do RTN é a lista de alvos com o modo, para a Tabela S de membros dos regulons.)
#
# Saídas em OUT_DIR:
#   TableS_ORA_per_regulon.csv        GO BP + KEGG por regulon (universo = rede)      [Tabela S "SX"]
#   TableS_ORA_universe_summary.csv   nº de termos significativos: universo rede x padrão  [SX, 2ª parte]
#   TableS_ORA_cluster_by_sign.csv    união dos alvos por cluster e sinal (A+/A-)    [SY]
#   TableS_ORA_recurrence.csv         termos recorrentes entre os regulons do cluster  [SZ]
#   TableS_regulon_validation.csv     ChIP-Atlas + DoRothEA/CollecTRI lado a lado     [SW]
#   TableS_regulon_members.csv        alvos de cada regulon com o modo (R3.13)        [SV]
#   SupplementaryTables_section3.4.xlsx  todas acima, uma aba cada (se o openxlsx estiver instalado)
#   Figure6.(tiff|png|pdf)            dotplot: clusters x alvos A+/A-
#   + no console: os números que vão no texto da 3.4

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

## ---------------------------------------------------------------- CONFIG
IN_DIR   <- "~/Kidney/results_regulon_validation"
OUT_DIR  <- "~/Kidney/results_regulon_validation/section3.4"
RTN_FILE <- "~/Kidney/0_data/rtni_TF.RData"; RTN_OBJECT <- "rtni_TF"   # NULL para pular a tabela de membros
SYMBOL_FILE <- "~/Kidney/0_data/ensembl_to_symbol.csv"
PADJ     <- 0.05
TOP_N    <- 6            # termos por grupo (cluster x sinal) na Figura 6
FORCE_TERMS <- c("lymphocyte homeostasis")   # termos citados no texto que devem aparecer na figura
# (se forem significativos); c() para nenhum

# Nomes que aparecem a partir da seção 3.4
CL_NAME <- c(Cluster_1 = "Proliferative\n(Cluster 1)", Cluster_2 = "Differentiation-associated\n(Cluster 2)")
CL_TXT  <- c(Cluster_1 = "Proliferative (Cluster 1)", Cluster_2 = "Differentiation-associated (Cluster 2)")
CL_COL  <- c(Cluster_1 = "#542788", Cluster_2 = "#1a9850")
TF14 <- list(Cluster_1 = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A"),
             Cluster_2 = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"))

pe <- path.expand
rd <- function(f) { p <- file.path(pe(IN_DIR), f); if (file.exists(p)) read.csv(p, check.names = FALSE) else {
  message("[aviso] não encontrado: ", f); NULL } }
dir.create(pe(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
wr <- function(d, f) if (!is.null(d)) write.csv(d, file.path(pe(OUT_DIR), f), row.names = FALSE)
fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
# aceita os rótulos antigos e novos dos CSVs ("Immune cluster", "Proliferative cluster (+)", "Cluster_2 (-)", "C1"...)
norm_cl <- function(x) {
  y <- trimws(sub("\\s*\\((\\+|-|\u2212)\\)\\s*$", "", as.character(x)))
  out <- ifelse(grepl("immune|differentiation|cluster[_ ]?2|^C2$", y, ignore.case = TRUE), "Cluster_2",
                ifelse(grepl("prolif|cluster[_ ]?1|^C1$", y, ignore.case = TRUE), "Cluster_1", NA))
  if (anyNA(out)) stop("Rótulo de cluster não reconhecido: ", paste(unique(x[is.na(out)]), collapse = ", "))
  out
}
cl_lab <- function(x) unname(CL_TXT[norm_cl(x)])
tabs <- list()

## ---------------------------------------------------------------- SX: ORA por regulon
ora <- rd("ORA_regulons_universo_rede.csv")
if (!is.null(ora)) {
  sx <- ora |> filter(p.adjust < PADJ) |>
    transmute(Regulon = tf, Cluster = cl_lab(cluster), Database = sub("GO_BP", "GO BP", database),
              ID, Term = Description, GeneRatio, BgRatio, p = pvalue, `Adjusted p (BH)` = p.adjust,
              Genes = geneID) |>
    arrange(Cluster, Regulon, Database, `Adjusted p (BH)`)
  tabs$SX_ORA_per_regulon <- sx; wr(sx, "TableS_ORA_per_regulon.csv")
}
us <- rd("ORA_resumo_universos.csv")
if (!is.null(us)) {
  us <- us |> transmute(Regulon = tf, Cluster = cl_lab(cluster), Database = sub("GO_BP", "GO BP", database),
                        `Significant terms, network background` = terms_sig_universe_network,
                        `Significant terms, genome-wide background` = terms_sig_universe_default)
  tabs$SX_universe_summary <- us; wr(us, "TableS_ORA_universe_summary.csv")
}

## ---------------------------------------------------------------- SY: união dos alvos por sinal
sg <- rd("ORA_por_cluster_sinal.csv")
if (!is.null(sg)) {
  sg <- sg |> mutate(cl = norm_cl(Cluster),
                     sign = ifelse(grepl("\\(\\+\\)\\s*$", Cluster), "A+", "A-"),
                     group = paste0(cl, "_", sign))
  sy <- sg |> filter(p.adjust < PADJ) |>
    transmute(Cluster = cl_lab(cl), Targets = ifelse(sign == "A+", "Positive (A+)", "Negative (A-)"),
              ID, Term = Description, GeneRatio, BgRatio, Count, p = pvalue, `Adjusted p (BH)` = p.adjust,
              Genes = geneID) |>
    arrange(Cluster, desc(Targets), `Adjusted p (BH)`)
  tabs$SY_cluster_by_sign <- sy; wr(sy, "TableS_ORA_cluster_by_sign.csv")
}

## ---------------------------------------------------------------- SZ: recorrência
rec <- rd("ORA_recorrencia_por_cluster.csv")
if (!is.null(rec)) {
  sz <- rec |> transmute(Cluster = cl_lab(cluster), Database = sub("GO_BP", "GO BP", database), ID,
                         Term = Description, `Regulons with the term` = n_regulons,
                         `Fraction of cluster regulons` = round(fraction, 2), Regulons = regulons,
                         `Lowest adjusted p` = min_padj) |>
    arrange(Cluster, desc(`Regulons with the term`), `Lowest adjusted p`)
  tabs$SZ_recurrence <- sz; wr(sz, "TableS_ORA_recurrence.csv")
}

## ---------------------------------------------------------------- SW: ChIP-Atlas + DoRothEA/CollecTRI
chip <- rd("ChIPAtlas_sobreposicao.csv"); ct <- rd("CollecTRI_sobreposicao.csv")
if (!is.null(chip) || !is.null(ct)) {
  pick <- function(d, pref, extra = NULL) {
    if (is.null(d)) return(NULL)
    keep <- c("tf", extra, "ref_targets_in_universe", "overlap", "pct_regulon_supported", "odds_ratio", "p", "FDR")
    d <- d[, intersect(keep, names(d))]
    names(d)[-1] <- paste(pref, names(d)[-1]); d
  }
  base <- bind_rows(lapply(names(TF14), function(k) data.frame(tf = TF14[[k]], Cluster = cl_lab(k))))
  size_from <- function(d) if (!is.null(d) && all(c("tf", "regulon_size") %in% names(d))) d[, c("tf", "regulon_size")] else NULL
  regsize <- size_from(chip); if (is.null(regsize)) regsize <- size_from(ct)
  if (is.null(regsize)) regsize <- data.frame(tf = base$tf, regulon_size = NA_integer_)
  regsize <- regsize[!duplicated(regsize$tf), ]
  src <- if (!is.null(ct) && "fonte" %in% names(ct)) ct$fonte[1] else "DoRothEA/CollecTRI"
  sw <- base |> left_join(regsize, by = "tf")
  if (!is.null(chip)) sw <- left_join(sw, pick(chip, "ChIP-Atlas:", "chip_experiments"), by = "tf")
  if (!is.null(ct))   sw <- left_join(sw, pick(ct, paste0(src, ":"), "sign_concordance"), by = "tf")
  names(sw)[names(sw) == "tf"] <- "Regulon"
  names(sw)[names(sw) == "regulon_size"] <- "Regulon size"
  names(sw) <- sub("ref_targets_in_universe", "reference targets in network", names(sw))
  names(sw) <- sub("pct_regulon_supported", "% of regulon supported", names(sw))
  names(sw) <- sub("chip_experiments", "experiments", names(sw))
  names(sw) <- sub("sign_concordance", "mode concordance", names(sw))
  tabs$SW_regulon_validation <- sw; wr(sw, "TableS_regulon_validation.csv")
}

## ---------------------------------------------------------------- SV: membros dos regulons (R3.13)
if (!is.null(RTN_FILE) && file.exists(pe(RTN_FILE))) {
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(pe(RTN_FILE), envir = env); rtni <- get(RTN_OBJECT, envir = env)
  ids <- tni.get(rtni, what = "regulatoryElements")          # nomes = símbolos, valores = IDs
  reg <- tni.get(rtni, what = "regulons.and.mode")
  smap <- if (file.exists(pe(SYMBOL_FILE))) { s <- read.csv(pe(SYMBOL_FILE)); setNames(s$symbol, sub("\\..*", "", s$ensembl)) } else c()
  sv <- bind_rows(lapply(names(TF14), function(k) bind_rows(lapply(TF14[[k]], function(tf) {
    v <- reg[[unname(ids[tf])]]
    e <- sub("\\..*", "", names(v))
    data.frame(Regulon = tf, Cluster = cl_lab(k), `Target (Ensembl)` = e,
               `Target (symbol)` = unname(ifelse(e %in% names(smap), smap[e], NA)),
               Mode = ifelse(v > 0, "Positive (A+)", "Negative (A-)"),
               `Association weight (signed MI)` = signif(unname(v), 3), check.names = FALSE)
  }))))
  tabs$SV_regulon_members <- sv; wr(sv, "TableS_regulon_members.csv")
  cat("\nAlvos por regulon (A+ / A-):\n")
  print(sv |> count(Regulon, Mode) |> tidyr::pivot_wider(names_from = Mode, values_from = n), n = Inf)
}

if (requireNamespace("openxlsx", quietly = TRUE) && length(tabs)) {
  openxlsx::write.xlsx(tabs, file.path(pe(OUT_DIR), "SupplementaryTables_section3.4.xlsx"),
                       firstRow = TRUE, colWidths = "auto")
}

## ---------------------------------------------------------------- Figura 6
if (!is.null(sg)) {
  gr <- function(x) sapply(strsplit(x, "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]))
  s <- sg |> filter(p.adjust < PADJ) |> mutate(ratio = gr(GeneRatio))
  top <- s |> group_by(group) |> slice_min(p.adjust, n = TOP_N, with_ties = FALSE) |> ungroup()
  forced <- s |> filter(tolower(Description) %in% tolower(FORCE_TERMS))
  terms <- unique(c(top$Description, forced$Description))
  d <- s |> filter(Description %in% terms)
  # cada termo fica no bloco do grupo em que é mais significativo
  best <- d |> group_by(Description) |> slice_min(p.adjust, n = 1, with_ties = FALSE) |> ungroup()
  grp_lvl <- c("Cluster_1_A+", "Cluster_1_A-", "Cluster_2_A+", "Cluster_2_A-")
  best <- best |> mutate(g = factor(group, grp_lvl)) |> arrange(g, desc(p.adjust))
  wrap <- function(x) sapply(x, function(z) paste(strwrap(z, 45), collapse = "\n"))
  d <- d |> mutate(term = factor(wrap(Description), unique(wrap(rev(best$Description)))),
                   cl = factor(CL_NAME[cl], CL_NAME),
                   sign = factor(ifelse(sign == "A+", "Positive\ntargets (A+)", "Negative\ntargets (A-)"),
                                 c("Positive\ntargets (A+)", "Negative\ntargets (A-)")))
  fig6 <- ggplot(d, aes(sign, term)) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradient(low = "#9ecae1", high = "#08306b", name = expression(-log[10]~adj.~p)) +
    scale_size_area(max_size = 5, name = "Genes") +
    facet_grid(~ cl) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 7) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(color = "white", face = "bold"),
          legend.key.size = unit(3, "mm"))
  g6 <- ggplotGrob(fig6)                                        # faixa do facet com a cor do cluster
  for (i in grep("^strip-t", g6$layout$name)) {
    txt <- tryCatch(g6$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label, error = function(e) "")
    k <- names(CL_NAME)[CL_NAME == txt]
    bg <- grep("rect", sapply(g6$grobs[[i]]$grobs[[1]]$children, function(x) class(x)[1]))
    if (length(k) && length(bg)) g6$grobs[[i]]$grobs[[1]]$children[[bg[1]]]$gp$fill <- CL_COL[k]
  }
  W <- 140; H <- max(70, 8 + 4.2 * length(terms))               # mm; 140 = 1,5 coluna
  sv6 <- function(f, ...) ggsave(file.path(pe(OUT_DIR), f), g6, width = W, height = H, units = "mm", bg = "white", ...)
  sv6("Figure6.tiff", dpi = 500, compression = "lzw"); sv6("Figure6.png", dpi = 300); sv6("Figure6.pdf", device = cairo_pdf)
  cat(sprintf("\nFigura 6: %d termos, %.0f x %.0f mm\n", length(terms), W, H))
}

## ---------------------------------------------------------------- números para o texto
cat("\n================ NÚMEROS PARA A SEÇÃO 3.4 ================\n")
if (!is.null(us)) {
  tot <- us |> group_by(Database) |>
    summarise(rede = sum(`Significant terms, network background`),
              padrao = sum(`Significant terms, genome-wide background`), .groups = "drop")
  cat("Termos significativos somados entre regulons (rede vs padrão):\n"); print(tot)
  cat("(O 87% / 80% da carta é a fração dos termos do universo padrão que continuam significativos com o\n",
      " universo da rede; esse cálculo precisa dos dois universos linha a linha, que o CSV salvo não tem.\n",
      " Se quiser confirmar, rode o trecho 'retencao' no fim deste arquivo dentro da sessão do\n",
      " regulon_validation_enrichment.R, onde o objeto 'ora' completo existe.)\n")
}
if (!is.null(sg)) {
  cat("\nTermo principal de cada grupo (união dos alvos):\n")
  print(sg |> group_by(group) |> slice_min(p.adjust, n = 3, with_ties = FALSE) |>
          transmute(group, Description, GeneRatio, padj = fmt_p(p.adjust)) |> as.data.frame(), row.names = FALSE)
  cat("\nTermos imunes no Cluster 2 (A+):\n")
  imm <- "T cell|lymphocyte|leukocyte|immune|myeloid|B cell|natural killer|interferon|cytokine|antigen|inflammat"
  print(sg |> filter(group == "Cluster_2_A+", p.adjust < PADJ, grepl(imm, Description, ignore.case = TRUE),
                     !grepl("cytokinesis", Description)) |>
          transmute(Description, GeneRatio, padj = fmt_p(p.adjust)) |> head(10), row.names = FALSE)
}
if (!is.null(rec)) {
  cat("\nRecorrência: termos significativos em >= metade dos regulons do cluster:\n")
  print(rec |> filter(fraction >= 0.5) |> mutate(cluster = cl_lab(cluster)) |> count(cluster, database, name = "n_termos"))
}
if (!is.null(chip)) {
  ok <- chip |> filter(!is.na(p))
  sig <- ok |> filter(FDR < 0.01)
  cat(sprintf("\nChIP-Atlas: dados para %d/14; enriquecidos (FDR < 0.01): %d -> %s; OR %.1f–%.1f\n",
              nrow(ok), nrow(sig), paste(sig$tf, collapse = ", "),
              min(sig$odds_ratio), max(sig$odds_ratio)))
}
if (!is.null(ct)) {
  ok <- ct |> filter(!is.na(p))
  cat(sprintf("%s: dados para %d/14 (%s); FDR < 0.05: %s; concordância do modo: %s\n", src, nrow(ok),
              paste(ok$tf, collapse = ", "),
              paste(sprintf("%s (FDR = %s)", ok$tf[ok$FDR < 0.05], fmt_p(ok$FDR[ok$FDR < 0.05])), collapse = ", "),
              paste(sprintf("%s %.0f%%", ok$tf, 100 * ok$sign_concordance), collapse = ", ")))
}

## ---------------------------------------------------------------- (opcional) retencao
# Rode só na sessão do regulon_validation_enrichment.R (precisa do 'ora' com os dois universos):
# ret <- ora |> filter(p.adjust < 0.05) |> distinct(tf, database, universe, ID) |>
#   group_by(tf, database, ID) |> summarise(def = "default" %in% universe, net = "network" %in% universe, .groups = "drop") |>
#   filter(def) |> group_by(database) |> summarise(retidos = mean(net))
# ret