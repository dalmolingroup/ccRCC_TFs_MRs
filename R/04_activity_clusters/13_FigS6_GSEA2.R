### Figura suplementar: GSEA de duas caudas (tna.gsea2) dos 14 regulons da assinatura ###
# Revisão CSBJ-D-26-01039 (R2.M5; Resultados 3.3) — substitui a figura da versão original (11 TFs)
#
# Refaz o tna.gsea2 SÓ com os 14 TFs, para a correção de BH ser feita entre eles (é o que o texto
# reporta: 13/14 significativos, TFAP2A dES = -0,32, p aj. = 0,23). O fenótipo é o mesmo da MRA
# (log2FC tumor x normal pareado no ICGC-RECA, DESeq2).
#
# O gráfico é feito em ggplot (não usa tna.plot.gsea2, que salva um PDF por TF): um painel por TF com
# a curva de enriquecimento dos alvos positivos (A+) e negativos (A-), as posições dos alvos no ranking
# e o dES e o p ajustado da tna.gsea2.
#
# Saídas (OUT_DIR): FigS_GSEA2_14TFs.(tiff|png|pdf), GSEA2_14TFs.csv

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

## ---------------------------------------------------------------- CONFIG
MRA_RDATA  <- "~/Kidney/0_data/MRA.RData"   # .RData com o objeto TNA (o primeiro objeto TNA é usado)
OUT_DIR    <- "~/Kidney/Figures/FigS_GSEA2"
SEED       <- 2026
N_PERM     <- 1000
MIN_SIZE   <- 15
SIZE_FILTER <- "posORneg"   # "posANDneg" (padrão do RTN) exige >= 15 alvos positivos E >= 15 negativos e
                            # pode descartar TFs; confira a lista "testados" na saída
PHENO_LAB  <- "Genes ranked by log2FC, tumour vs normal (ICGC-RECA)"

TF14 <- data.frame(
  symbol  = c("FOXM1", "MYBL2", "NFE2L3", "PARP12", "TFAP2A",
              "IRF6", "MEF2D", "NR3C2", "PBX1", "RARB", "RORC", "SATB1", "TAL1", "TCF4"),
  ensembl = c("ENSG00000111206", "ENSG00000101057", "ENSG00000050344", "ENSG00000059378", "ENSG00000137203",
              "ENSG00000117595", "ENSG00000116604", "ENSG00000151623", "ENSG00000185630", "ENSG00000077092",
              "ENSG00000143365", "ENSG00000182568", "ENSG00000162367", "ENSG00000196628"),
  cluster = rep(c("Cluster 1", "Cluster 2"), c(5, 9)))
CL_COL  <- c("Cluster 1" = "#542788", "Cluster 2" = "#1a9850")
TAR_COL <- c("Positive targets (A+)" = "#b2182b", "Negative targets (A-)" = "#2166ac")

pe <- path.expand
strip_v <- function(x) sub("\\.\\d+$", "", x)

## ---------------------------------------------------------------- 1. GSEA2 com os 14 TFs
run_gsea2 <- function() {
  suppressPackageStartupMessages(library(RTN))
  env <- new.env(); load(pe(MRA_RDATA), envir = env)
  obj <- Filter(function(o) is(o, "TNA"), mget(ls(env), envir = env))
  if (!length(obj)) stop("Nenhum objeto TNA em ", MRA_RDATA)
  rtna <- obj[[1]]
  reg <- tna.get(rtna, what = "regulatoryElements")      # nomes e valores: símbolo / ID
  pick <- function(v) {                                   # entrada do regulatoryElements de cada TF
    hit <- which(strip_v(reg) %in% v | names(reg) %in% v | strip_v(names(reg)) %in% v)
    if (length(hit)) hit[1] else NA_integer_
  }
  idx <- sapply(seq_len(nrow(TF14)), function(i) pick(c(TF14$symbol[i], TF14$ensembl[i])))
  if (anyNA(idx)) stop("TFs fora do objeto TNA: ", paste(TF14$symbol[is.na(idx)], collapse = ", "))
  # tna.gsea2 aceita em 'tfs' o mesmo tipo usado em tna.mra; testa as duas formas
  args <- list(object = rtna, tfs = unname(reg[idx]), pValueCutoff = 1, pAdjustMethod = "BH",
               minRegulonSize = MIN_SIZE, sizeFilterMethod = SIZE_FILTER, nPermutations = N_PERM,
               verbose = FALSE)
  args <- args[names(args) %in% c("object", names(formals(tna.gsea2)))]
  set.seed(SEED)
  rtna2 <- tryCatch(do.call(tna.gsea2, args), error = function(e) NULL)
  res <- if (!is.null(rtna2)) tna.get(rtna2, what = "gsea2", ntop = -1) else NULL
  if (is.null(res) || !NROW(res$differential)) {
    args$tfs <- names(reg[idx]); set.seed(SEED)
    rtna2 <- do.call(tna.gsea2, args)
    res <- tna.get(rtna2, what = "gsea2", ntop = -1)
  }
  list(res = res,
       phen = tna.get(rtna2, what = "pheno"),
       regmode = tna.get(rtna2, what = "regulons.and.mode"),
       map = setNames(TF14$symbol, names(reg[idx])), map2 = setNames(TF14$symbol, unname(reg[idx])))
}
g <- if (exists("GSEA_INPUT")) GSEA_INPUT else run_gsea2()   # GSEA_INPUT só para teste

to_sym <- function(x) {
  s <- ifelse(x %in% names(g$map), g$map[x], ifelse(x %in% names(g$map2), g$map2[x], NA))
  s <- ifelse(is.na(s) & strip_v(x) %in% TF14$ensembl, TF14$symbol[match(strip_v(x), TF14$ensembl)], s)
  ifelse(is.na(s) & x %in% TF14$symbol, x, s)
}
col_of <- function(df, pat) df[[grep(pat, names(df), ignore.case = TRUE)[1]]]
tab_from <- function(df, tag) {
  df <- as.data.frame(df); id <- if ("Regulon" %in% names(df)) df$Regulon else rownames(df)
  data.frame(symbol = to_sym(id), score = col_of(df, "Observed|Score"),
             p = col_of(df, "^Pvalue"), padj = col_of(df, "Adjusted")) |> setNames(c("symbol", paste0(tag, c("", "_p", "_padj"))))
}
tab <- tab_from(g$res$differential, "dES") |>
  left_join(tab_from(g$res$positive, "ES_pos"), by = "symbol") |>
  left_join(tab_from(g$res$negative, "ES_neg"), by = "symbol") |>
  right_join(TF14, by = "symbol") |> relocate(symbol, ensembl, cluster)
dir.create(pe(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
write.csv(tab, file.path(pe(OUT_DIR), "GSEA2_14TFs.csv"), row.names = FALSE)

## ---------------------------------------------------------------- 2. curvas de enriquecimento
phen <- sort(g$phen[!is.na(g$phen)], decreasing = TRUE)
names(phen) <- strip_v(names(phen))
N <- length(phen)
running_es <- function(set) {                       # GSEA ponderada (expoente 1), como no RTN
  hit <- names(phen) %in% set
  if (!any(hit)) return(NULL)
  ph <- abs(phen) * hit; ph <- ph / sum(ph)
  ms <- (!hit) / (N - sum(hit))
  data.frame(rank = seq_len(N), es = cumsum(ph - ms), hit = hit)
}
curves <- list(); ticks <- list()
for (r in names(g$regmode)) {
  s <- to_sym(r); if (is.na(s)) next
  m <- g$regmode[[r]]; names(m) <- strip_v(names(m))
  for (sg in c("pos", "neg")) {
    set <- names(m)[if (sg == "pos") m > 0 else m < 0]
    rc <- running_es(set); if (is.null(rc)) next
    lab <- if (sg == "pos") names(TAR_COL)[1] else names(TAR_COL)[2]
    curves[[paste(s, sg)]] <- data.frame(symbol = s, targets = lab, rank = rc$rank, es = rc$es)
    ticks[[paste(s, sg)]]  <- data.frame(symbol = s, targets = lab, rank = rc$rank[rc$hit])
  }
}
cv <- bind_rows(curves); tk <- bind_rows(ticks)
cv <- cv[cv$rank %% 20 == 0 | cv$rank %in% c(1, N), ]      # amostra a curva (arquivo leve; forma igual)

lvl <- TF14$symbol
facet_lab <- setNames(sprintf("%s (%s)\ndES = %.2f, adj. p %s", tab$symbol, sub("Cluster ", "C", tab$cluster),
                              tab$dES, ifelse(tab$dES_padj < 0.001, "< 0.001", sprintf("= %.3f", tab$dES_padj))),
                      tab$symbol)
cv$symbol <- factor(cv$symbol, lvl); tk$symbol <- factor(tk$symbol, lvl)
yr <- range(cv$es); tick_y <- c(pos = yr[1] - 0.06 * diff(yr), neg = yr[1] - 0.14 * diff(yr))
tk$y <- ifelse(tk$targets == names(TAR_COL)[1], tick_y["pos"], tick_y["neg"])

pl <- ggplot(cv, aes(rank, es, color = targets)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_segment(data = tk, aes(x = rank, xend = rank, y = y - 0.035 * diff(yr), yend = y + 0.035 * diff(yr)),
               linewidth = 0.1, alpha = 0.6) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = TAR_COL, name = NULL) +
  scale_x_continuous(breaks = c(1, N), labels = c("Up in tumour", "Down in tumour"), expand = c(0.01, 0)) +
  facet_wrap(~ symbol, ncol = 5, labeller = as_labeller(facet_lab)) +
  labs(x = PHENO_LAB, y = "Running enrichment score") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.text = element_text(size = 6, color = "white"),
        axis.text.x = element_text(hjust = c(0, 1)))

# cor da faixa de cada painel pelo cluster (sem pacote extra)
gt <- ggplotGrob(pl)
strips <- grep("^strip-t", gt$layout$name)
for (i in strips) {
  lab_txt <- tryCatch(gt$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label, error = function(e) NULL)
  if (is.null(lab_txt)) next
  s <- sub(" .*", "", lab_txt)
  col <- CL_COL[TF14$cluster[match(s, TF14$symbol)]]
  bg <- grep("rect", sapply(gt$grobs[[i]]$grobs[[1]]$children, function(x) class(x)[1]))
  if (length(bg) && !is.na(col)) gt$grobs[[i]]$grobs[[1]]$children[[bg[1]]]$gp$fill <- col
}

W <- 190; H <- 135
save_all <- function(f, ...) ggsave(file.path(pe(OUT_DIR), f), gt, width = W, height = H, units = "mm", bg = "white", ...)
save_all("FigS_GSEA2_14TFs.tiff", dpi = 500, compression = "lzw")
save_all("FigS_GSEA2_14TFs.png", dpi = 300)
save_all("FigS_GSEA2_14TFs.pdf", device = cairo_pdf)

## ---------------------------------------------------------------- 3. resumo
cat("TFs testados na tna.gsea2:", sum(!is.na(tab$dES)), "de 14\n")
if (any(is.na(tab$dES))) cat("Não testados (filtro de tamanho?):", paste(tab$symbol[is.na(tab$dES)], collapse = ", "), "\n")
cat(sprintf("Significativos (p aj. < 0.05): %d | maior p aj. entre eles = %.3f\n",
            sum(tab$dES_padj < 0.05, na.rm = TRUE), max(tab$dES_padj[tab$dES_padj < 0.05], na.rm = TRUE)))
print(tab |> transmute(symbol, cluster, dES = round(dES, 2), padj = signif(dES_padj, 2)), row.names = FALSE)
cat("(Esperado pela carta: 13/14 com p aj. < 0,02; TFAP2A dES = -0,32, p aj. = 0,23)\n")
