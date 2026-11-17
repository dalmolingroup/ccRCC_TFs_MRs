### Validação dos regulons e enriquecimento com universo definido — R2.M6, R2.a8 (R3.4, R3.7, R3.13) ###
# R2.M6: "validating key interactions using available ChIP-seq data, transcription-factor binding
#         databases [...] clarify the background gene sets used for enrichment"
# R2.a8: "specify the background gene sets used for GO and KEGG analyses"
#
# Parte A — ORA por regulon (regulador + alvos; GO BP e KEGG), universo = genes da rede do RTN,
#           comparada com o universo padrão do clusterProfiler (retenção dos termos)
# Parte B — sobreposição dos alvos com interações curadas (CollecTRI/DoRothEA) + concordância do modo
# Parte C — sobreposição com alvos de ChIP-seq do ChIP-Atlas (TSS ± 5 kb)
#           (B e C: Fisher unilateral, universo = genes da rede; FDR (BH) entre os TFs)
# Parte D — recorrência dos termos entre os regulons de cada cluster
# Parte E — ORA da união dos alvos de cada cluster: todos e separados por sinal (A+ / A-)
# Parte F — (opcional) cobertura dos alvos no ICGC-RECA
#
# Clusters: Cluster_1 = proliferativo; Cluster_2 = diferenciação. Os nomes por extenso entram
# no section34_tables_fig6.R, que lê as saídas deste script e monta as Tabelas S e a Figura 6.

suppressPackageStartupMessages({
  library(dplyr)
  library(RTN)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## ---------------------------------------------------------------- CONFIG
RTN_FILE    <- "~/Kidney/0_data/rtni_TF.RData"
RTN_OBJECT  <- "rtni_TF"
SYMBOL_FILE <- "~/Kidney/0_data/ensembl_to_symbol.csv"     # gerado pelo figure5_expression.R
SIG_FILE    <- "results_ml_revision/importancia_features_modelo_final_class_weight.csv"  # 1ª coluna = ENSG
OUT_DIR     <- "~/Kidney/results_regulon_validation"

CLUSTERS <- list(
  Cluster_1 = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A"),
  Cluster_2 = c("IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"))

RUN_ORA       <- TRUE
RUN_COLLECTRI <- TRUE
RUN_CHIPATLAS <- TRUE
RUN_CLUSTER   <- TRUE
COLLECTRI_FILE <- NULL            # opcional: CSV com colunas source,target,mor (se não houver internet)
CHIPATLAS_DIR  <- NULL            # opcional: pasta com <TF>.5.tsv já baixados
CHIPATLAS_URL  <- "https://chip-atlas.dbcls.jp/data/hg38/target/%s.5.tsv"   # ±5 kb do TSS
CHIPATLAS_MIN_SCORE <- 0          # alvo = gene com média de score MACS2 > este valor
ORA_PCUT <- 0.05

# Parte F (opcional): .RData com o objeto TNI do ICGC usado no tni.replace.samples
ICGC_RDATA  <- NULL               # ex.: "~/Kidney/0_data/rtni1st_icgcsamples.RData"
ICGC_OBJECT <- "rtni1st_icgcsamples"

pe <- path.expand
dir.create(pe(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
outf <- function(f) file.path(pe(OUT_DIR), f)
strip_v <- function(x) sub("\\..*", "", x)

## ---------------------------------------------------------------- rede, universo e alvos
env <- new.env(); load(pe(RTN_FILE), envir = env)
rtni <- get(RTN_OBJECT, envir = env)
ids  <- tni.get(rtni, what = "regulatoryElements")            # nomes = símbolos, valores = IDs
reg_mode <- tni.get(rtni, what = "regulons.and.mode")         # lista: ID do TF -> vetor nomeado (sinal = modo)
universe_ens <- strip_v(rownames(tni.get(rtni, "gexp")))

sym_map <- read.csv(pe(SYMBOL_FILE))
sym_map <- setNames(sym_map$symbol, strip_v(sym_map$ensembl))
to_sym <- function(x) unname(ifelse(strip_v(x) %in% names(sym_map), sym_map[strip_v(x)], x))

# assinatura: IDs do RTN dos 14 TFs (lidos da tabela de importância e conferidos com os clusters)
sig_ens <- strip_v(read.csv(SIG_FILE, check.names = FALSE)[[1]])
sig <- unname(ids[strip_v(ids) %in% sig_ens])                 # IDs exatamente como no RTN
names(sig) <- names(ids)[strip_v(ids) %in% sig_ens]           # símbolos
if (length(sig) != length(sig_ens)) stop("Nem todos os TFs da assinatura estão no RTN: ",
                                         paste(setdiff(sig_ens, strip_v(sig)), collapse = ", "))
if (!setequal(names(sig), unlist(CLUSTERS))) stop("TFs da assinatura e de CLUSTERS não batem: ",
                                                  paste(sort(names(sig)), collapse = ", "))
cluster_of <- setNames(rep(names(CLUSTERS), lengths(CLUSTERS)), unlist(CLUSTERS))
sig <- sig[unlist(CLUSTERS)]                                  # ordem: Cluster_1, Cluster_2

universe_sym <- unique(to_sym(universe_ens))
targets <- lapply(sig, function(id) {
  v <- reg_mode[[id]]
  data.frame(target = to_sym(names(v)), mode = sign(v))
})
info <- data.frame(tf = names(sig), cluster = unname(cluster_of[names(sig)]),
                   regulon_size = vapply(targets, nrow, 1L), row.names = NULL)
message("Universo (genes da rede): ", length(universe_sym), " | regulons: ",
        paste(sprintf("%s (%d)", info$tf, info$regulon_size), collapse = ", "))

## ---------------------------------------------------------------- utilitário: Fisher contra o universo
overlap_test <- function(reg_targets, ref_targets, universe) {
  reg <- intersect(reg_targets, universe); ref <- intersect(ref_targets, universe)
  a <- length(intersect(reg, ref)); b <- length(setdiff(reg, ref))
  c <- length(setdiff(ref, reg)); d <- length(universe) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater")
  data.frame(ref_targets_in_universe = length(ref), overlap = a,
             pct_regulon_supported = 100 * a / max(length(reg), 1),
             expected = length(reg) * length(ref) / length(universe),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}
to_entrez <- function(ens) unique(suppressMessages(bitr(strip_v(ens), "ENSEMBL", "ENTREZID", org.Hs.eg.db))$ENTREZID)
universe_ez <- to_entrez(universe_ens)

## ---------------------------------------------------------------- Parte A: ORA por regulon
if (RUN_ORA) {
  ora_rows <- list()
  for (tf in names(sig)) {
    genes <- to_entrez(c(sig[[tf]], names(reg_mode[[sig[[tf]]]])))   # regulador + alvos, como no artigo
    run <- function(fun, univ, ...) {
      r <- tryCatch(fun(gene = genes, universe = univ, pvalueCutoff = 1, qvalueCutoff = 1, ...), error = function(e) NULL)
      if (is.null(r) || !nrow(as.data.frame(r))) return(NULL)
      as.data.frame(setReadable(r, org.Hs.eg.db, keyType = "ENTREZID"))
    }
    for (db in c("GO_BP", "KEGG")) for (u in c("network", "default")) {
      univ <- if (u == "network") universe_ez else NULL
      r <- if (db == "GO_BP") run(enrichGO, univ, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
      else run(enrichKEGG, univ, organism = "hsa", pAdjustMethod = "BH")
      if (!is.null(r)) ora_rows[[length(ora_rows) + 1]] <- cbind(tf = tf, cluster = cluster_of[[tf]],
                                                                 database = db, universe = u, r)
    }
  }
  ora <- bind_rows(ora_rows)
  write.csv(ora |> filter(universe == "network"), outf("ORA_regulons_universo_rede.csv"), row.names = FALSE)
  write.csv(ora, outf("ORA_regulons_dois_universos.csv"), row.names = FALSE)
  
  ora_sum <- ora |> group_by(tf, cluster, database, universe) |>
    summarise(n_sig = sum(p.adjust < ORA_PCUT), .groups = "drop") |>
    tidyr::pivot_wider(names_from = universe, values_from = n_sig, values_fill = 0,
                       names_prefix = "terms_sig_universe_")
  write.csv(ora_sum, outf("ORA_resumo_universos.csv"), row.names = FALSE)
  
  # retenção: fração dos termos significativos com o universo padrão que continuam com o da rede
  ret <- ora |> filter(p.adjust < ORA_PCUT) |> distinct(tf, database, universe, ID) |>
    group_by(tf, database, ID) |>
    summarise(def = "default" %in% universe, net = "network" %in% universe, .groups = "drop") |>
    filter(def) |> group_by(database) |>
    summarise(terms_default = n(), retained_with_network = sum(net), pct_retained = round(100 * mean(net), 1))
  write.csv(ret, outf("ORA_retencao_universos.csv"), row.names = FALSE)
  
  writeLines(c(paste("clusterProfiler", packageVersion("clusterProfiler")),
               paste("org.Hs.eg.db", packageVersion("org.Hs.eg.db")),
               paste("GO.db", packageVersion("GO.db")),
               paste("KEGG acessado em", Sys.Date()),
               paste("Universo da rede:", length(universe_ez), "genes com Entrez")),
             outf("ORA_versoes.txt"))
  cat("\n=== ORA: termos significativos (BH <", ORA_PCUT, ") por regulon, universo da rede vs padrão ===\n")
  print(ora_sum, n = 60)
  cat("\n=== Retenção dos termos ao trocar o universo (vai no texto: 87% GO / 80% KEGG?) ===\n"); print(ret)
}

## ---------------------------------------------------------------- Parte B: CollecTRI / DoRothEA
if (RUN_COLLECTRI) {
  # fontes em ordem: arquivo local -> decoupleR -> OmnipathR -> DoRothEA (offline, níveis A-C)
  std <- function(d) {
    d <- as.data.frame(d)
    if (!nrow(d)) return(NULL)
    names(d)[names(d) %in% c("tf", "source_genesymbol")] <- "source"
    names(d)[names(d) %in% c("target_genesymbol")] <- "target"
    if (!"mor" %in% names(d)) {
      d$mor <- if (all(c("is_stimulation", "is_inhibition") %in% names(d)))
        ifelse(d$is_inhibition == 1 & d$is_stimulation == 0, -1, 1) else 1
    }
    if (!all(c("source", "target") %in% names(d))) return(NULL)
    d[, c("source", "target", "mor")]
  }
  try_src <- function(expr, src) {
    out <- tryCatch(std(expr), error = function(e) { message("[", src, "] falhou: ", conditionMessage(e)); NULL })
    if (!is.null(out)) attr(out, "fonte") <- src
    out
  }
  ct <- NULL
  if (!is.null(COLLECTRI_FILE)) ct <- try_src(read.csv(pe(COLLECTRI_FILE)), "arquivo local")
  if (is.null(ct) && requireNamespace("decoupleR", quietly = TRUE))
    ct <- try_src(decoupleR::get_collectri(organism = "human", split_complexes = FALSE), "CollecTRI (decoupleR)")
  if (is.null(ct) && requireNamespace("OmnipathR", quietly = TRUE))
    ct <- try_src(OmnipathR::collectri(organism = 9606, genesymbols = TRUE, loops = TRUE), "CollecTRI (OmnipathR)")
  if (is.null(ct) && requireNamespace("dorothea", quietly = TRUE)) {
    data("dorothea_hs", package = "dorothea", envir = environment())
    ct <- try_src(subset(dorothea_hs, confidence %in% c("A", "B", "C")), "DoRothEA (A-C)")
  }
  if (is.null(ct)) stop("Nenhuma fonte de interações TF-alvo disponível. Instale o pacote 'dorothea' ",
                        "(BiocManager::install('dorothea')) ou baixe o CollecTRI e use COLLECTRI_FILE.")
  fonte <- attr(ct, "fonte")
  message("Fonte das interações TF-alvo: ", fonte, " (", nrow(ct), " interações)")
  ct_tab <- bind_rows(lapply(info$tf, function(tf) {
    ref <- ct[ct$source == tf, ]
    if (!nrow(ref)) return(cbind(info[info$tf == tf, ], data.frame(ref_targets_in_universe = 0, overlap = NA,
                                                                   pct_regulon_supported = NA, expected = NA, odds_ratio = NA, p = NA, sign_concordance = NA)))
    sh <- merge(targets[[tf]], ref, by = "target")
    cbind(info[info$tf == tf, ], overlap_test(targets[[tf]]$target, ref$target, universe_sym),
          sign_concordance = if (nrow(sh)) mean(sign(sh$mode) == sign(sh$mor)) else NA)
  })) |> mutate(FDR = p.adjust(p, "BH"), fonte = fonte)
  write.csv(ct_tab, outf("CollecTRI_sobreposicao.csv"), row.names = FALSE)
  cat("\n=== Interações TF-alvo:", fonte, "===\n"); print(ct_tab |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
}

## ---------------------------------------------------------------- Parte C: ChIP-Atlas
if (RUN_CHIPATLAS) {
  read_chip <- function(tf) {
    f <- if (!is.null(CHIPATLAS_DIR)) file.path(pe(CHIPATLAS_DIR), paste0(tf, ".5.tsv")) else sprintf(CHIPATLAS_URL, tf)
    d <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(d) || ncol(d) < 2) return(NULL)
    list(genes = d[[1]][d[[2]] > CHIPATLAS_MIN_SCORE], n_exp = ncol(d) - 2)
  }
  chip_tab <- bind_rows(lapply(info$tf, function(tf) {
    ch <- read_chip(tf)
    if (is.null(ch)) return(cbind(info[info$tf == tf, ], data.frame(chip_experiments = 0, ref_targets_in_universe = 0,
                                                                    overlap = NA, pct_regulon_supported = NA, expected = NA, odds_ratio = NA, p = NA)))
    cbind(info[info$tf == tf, ], chip_experiments = ch$n_exp, overlap_test(targets[[tf]]$target, ch$genes, universe_sym))
  })) |> mutate(FDR = p.adjust(p, "BH"))
  write.csv(chip_tab, outf("ChIPAtlas_sobreposicao.csv"), row.names = FALSE)
  cat("\n=== ChIP-Atlas (TSS ± 5 kb) ===\n"); print(chip_tab |> mutate(across(where(is.numeric), \(x) signif(x, 3))))
}

## ---------------------------------------------------------------- Parte D: recorrência por cluster
ora_net <- read.csv(outf("ORA_regulons_universo_rede.csv"))
n_tf <- lengths(CLUSTERS)
rec <- ora_net |> filter(p.adjust < ORA_PCUT) |>
  group_by(cluster, database, ID, Description) |>
  summarise(n_regulons = n_distinct(tf), regulons = paste(sort(unique(tf)), collapse = ", "),
            min_padj = min(p.adjust), .groups = "drop") |>
  mutate(fraction = n_regulons / n_tf[cluster]) |>
  arrange(cluster, desc(n_regulons), min_padj)
write.csv(rec, outf("ORA_recorrencia_por_cluster.csv"), row.names = FALSE)
cat("\n=== Termos mais recorrentes por cluster ===\n")
rec |> group_by(cluster) |> slice_head(n = 10) |> print(n = Inf)

## ---------------------------------------------------------------- Parte E: união dos alvos por cluster
if (RUN_CLUSTER) {
  tf_ids_of <- function(cl) unname(ids[names(ids) %in% CLUSTERS[[cl]]])
  # (1) todos os alvos (reguladores + união)
  gene_sets <- lapply(setNames(names(CLUSTERS), names(CLUSTERS)), function(cl) {
    tf_ids <- tf_ids_of(cl)
    to_entrez(unique(c(tf_ids, unlist(lapply(reg_mode[tf_ids], names)))))
  })
  print(lengths(gene_sets))
  cc <- compareCluster(gene_sets, fun = "enrichGO", universe = universe_ez, OrgDb = org.Hs.eg.db,
                       ont = "BP", pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, readable = TRUE)
  write.csv(as.data.frame(cc), outf("ORA_por_cluster_uniao.csv"), row.names = FALSE)
  
  # (2) separados pelo sinal do modo: A+ e A- (unname evita colar o nome do TF no do alvo)
  gene_sets_signed <- unlist(lapply(names(CLUSTERS), function(cl) {
    v <- do.call(c, unname(reg_mode[tf_ids_of(cl)]))
    pos <- unique(names(v)[v > 0]); neg <- unique(names(v)[v < 0])
    setNames(list(to_entrez(pos), to_entrez(neg)), paste(cl, c("(+)", "(-)")))
  }), recursive = FALSE)
  cat("\nGenes por conjunto (devem ser centenas):\n"); print(lengths(gene_sets_signed))
  cc_s <- compareCluster(gene_sets_signed, fun = "enrichGO", universe = universe_ez, OrgDb = org.Hs.eg.db,
                         ont = "BP", pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, readable = TRUE)
  write.csv(as.data.frame(cc_s), outf("ORA_por_cluster_sinal.csv"), row.names = FALSE)
  cat("\n=== União dos alvos por sinal: termos principais ===\n")
  as.data.frame(cc_s) |> group_by(Cluster) |> slice_min(p.adjust, n = 8) |>
    dplyr::select(Cluster, Description, GeneRatio, p.adjust) |> print(n = Inf)
  
  imm <- "T cell|lymphocyte|leukocyte|immune|myeloid|B cell|natural killer|interferon|cytokine|antigen|inflammat"
  cat("\n=== Termos imunes no Cluster_2 (por regulon e por sinal) ===\n")
  ora_net |> filter(cluster == "Cluster_2", grepl(imm, Description, ignore.case = TRUE),
                    !grepl("cytokinesis|giant", Description, ignore.case = TRUE), p.adjust < ORA_PCUT) |>
    arrange(p.adjust) |> dplyr::select(tf, database, Description, p.adjust) |> print()
  as.data.frame(cc_s) |> filter(grepl("^Cluster_2", Cluster), grepl(imm, Description, ignore.case = TRUE)) |>
    dplyr::select(Cluster, Description, GeneRatio, p.adjust) |> print()
}

## ---------------------------------------------------------------- Parte F: cobertura no ICGC (opcional)
if (!is.null(ICGC_RDATA) && file.exists(pe(ICGC_RDATA))) {
  e2 <- new.env(); load(pe(ICGC_RDATA), envir = e2)
  icgc_genes <- strip_v(rownames(tni.get(get(ICGC_OBJECT, envir = e2), "gexp")))
  cov_df <- bind_rows(lapply(names(sig), function(s) {
    t <- strip_v(names(reg_mode[[sig[[s]]]]))
    data.frame(TF = s, targets = length(t), in_ICGC = sum(t %in% icgc_genes),
               pct = round(100 * mean(t %in% icgc_genes), 1))
  }))
  write.csv(cov_df, outf("cobertura_alvos_ICGC.csv"), row.names = FALSE); print(cov_df)
}

cat("\nSaídas em", pe(OUT_DIR), "\n")