### Figura 5 - montagem final no formato da CSBJ (largura de 2 colunas = 7 in) ###
# Revisão CSBJ-D-26-01039 - pedido do editor: reduzir o espaço em branco da Figura 5.
# Guia da revista: 1 coluna = 3,4 in, 2 colunas = 7 in; >= 300 dpi no tamanho final; RGB;
# tipos de tamanho próximo (~7 pt); evitar espaço em branco; preferir arquivos vetoriais.
#
# Layout:  A  tumor x normal pareado, 14 TFs em 7 colunas x 2 linhas (largura total)
#          B  M1 x M0 (só tumores), 14 TFs em 5 colunas x 3 linhas; legenda no 15º espaço
#          C  concordância do Cliff's delta TCGA x ICGC, ao lado de B
#
# Entradas (rodar antes): 14_expression_Fig5AB_TabS6.R  -> Figures/Figure5/fig5_data.rds
#                         15_Fig5C_concordance.R          -> Figures/Figure5/concordancia_TCGA_ICGC.csv
# Saídas: Figures/Figure5/Figure5.(pdf|tiff|png)  (o PDF é vetorial: envie este à revista)

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(ggh4x)
})

FIG_DIR <- path.expand("~/Kidney/Figures/Figure5")
dat  <- readRDS(file.path(FIG_DIR, "fig5_data.rds"))
conc <- read.csv(file.path(FIG_DIR, "concordancia_TCGA_ICGC.csv"))

W <- 7; H <- 6.4                 # polegadas (2 colunas)
COH <- c("TCGA-KIRC" = "TCGA", "ICGC-RECA" = "ICGC")
CL_LAB <- c(cluster_1 = "Proliferative", cluster_2 = "Differentiation-associated")

## ---------------------------------------------------------------- painéis A e B (boxplots)
panel_box <- function(df, st, groups, ncol, legend_inside = NULL) {
  lev <- dat$gene_order
  df <- df |> mutate(gene = factor(gene, levels = lev), type = factor(type, levels = groups),
                     coh = factor(COH[cohort], levels = COH))
  st <- st |> mutate(gene = factor(gene, levels = lev), coh = factor(COH[cohort], levels = COH))
  br <- df |> group_by(gene, coh) |>
    summarise(top = max(expression), rng = diff(range(expression)), .groups = "drop") |>
    mutate(y = top + 0.06 * rng, x = as.integer(coh)) |>
    left_join(select(st, gene, coh, signif), by = c("gene", "coh"))
  strip_fill <- unname(dat$CLUSTER_COLORS[dat$cluster_of[lev]])

  p <- ggplot(df, aes(coh, expression, fill = type)) +
    geom_boxplot(outlier.size = 0.25, outlier.stroke = 0, width = 0.65,
                 position = position_dodge(0.75), linewidth = 0.25) +
    geom_segment(data = br, aes(x = x - 0.3, xend = x + 0.3, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_text(data = br, aes(x = x, y = y, label = signif), inherit.aes = FALSE,
              vjust = -0.15, size = 1.9) +
    scale_fill_manual(values = dat$GROUP_COLORS[groups], name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.24)), n.breaks = 4) +
    labs(x = NULL, y = "Expression (log2 normalised counts)") +
    facet_wrap2(~gene, ncol = ncol, scales = "free_y",
                strip = strip_themed(background_x = elem_list_rect(fill = strip_fill, colour = NA))) +
    theme_bw(base_size = 7) +
    theme(strip.text = element_text(face = "bold", colour = "white", size = 6.5,
                                    margin = margin(1.5, 0, 1.5, 0)),
          axis.text = element_text(size = 5.5, colour = "black"),
          axis.title.y = element_text(size = 6.5),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.spacing = unit(2.5, "pt"),
          legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 6.5),
          legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 3, 2, 2))
  if (is.null(legend_inside)) {
    p + theme(legend.position = "top", legend.box.spacing = unit(1, "pt"))
  } else {
    p + theme(legend.position = "inside", legend.position.inside = legend_inside,
              legend.direction = "vertical")
  }
}

pA <- panel_box(dat$plot_tn, dat$stats_tn, c("Normal", "Tumor"), ncol = 7)
pB <- panel_box(dat$plot_m,  dat$stats_m,  c("M0", "M1"),       ncol = 5,
                legend_inside = c(0.90, 0.16))                    # 15º espaço, vazio, da grade 5 x 3

## ---------------------------------------------------------------- painel C (concordância)
conc <- conc |> mutate(cl = CL_LAB[unname(dat$cluster_of[gene])],
                       shape = ifelse(concordant, "Concordant", "Discordant"))
ct  <- cor.test(conc$delta_TCGA, conc$delta_ICGC, method = "spearman", exact = FALSE)
# letras gregas via plotmath (independe do locale/fonte do servidor)
lab1 <- sprintf("Same direction: %d/%d", sum(conc$concordant), nrow(conc))
lab2 <- sprintf("'Spearman'~rho == '%.2f'*','~italic(p) == '%.3f'", ct$estimate, ct$p.value)
lims <- range(c(conc$delta_TCGA, conc$delta_ICGC, 0)) * 1.15

txt <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = gene), size = 2, show.legend = FALSE, min.segment.length = 0.2,
                           segment.size = 0.2, box.padding = 0.15, point.padding = 0.1, seed = 2026)
} else geom_text(aes(label = gene), size = 2, vjust = -0.8, show.legend = FALSE)

pC <- ggplot(conc, aes(delta_TCGA, delta_ICGC, colour = cl)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#f2f1ec") +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#f2f1ec") +
  geom_hline(yintercept = 0, colour = "#8a8983", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "#8a8983", linewidth = 0.3) +
  geom_point(aes(shape = shape), size = 1.6, stroke = 0.5) + txt +
  annotate("text", x = Inf, y = -Inf, label = lab1, hjust = 1.04, vjust = -2.2, size = 2,
           colour = "#52514e") +                                # canto inferior direito (vazio)
  annotate("text", x = Inf, y = -Inf, label = lab2, parse = TRUE, hjust = 1.04, vjust = -0.4, size = 2,
           colour = "#52514e") +
  scale_colour_manual(values = setNames(unname(dat$CLUSTER_COLORS[names(CL_LAB)]), CL_LAB), name = NULL) +
  scale_shape_manual(values = c(Concordant = 16, Discordant = 1), name = NULL) +
  coord_cartesian(xlim = lims, ylim = lims) +
  labs(x = expression("Cliff's" ~ delta * ", M1 vs M0 (TCGA-KIRC)"),
       y = expression("Cliff's" ~ delta * ", M1 vs M0 (ICGC-RECA)")) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 5.5, colour = "black"),
        legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"),
        legend.key.size = unit(2.5, "mm"), legend.text = element_text(size = 6.5),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(2, 2, 2, 4))

## ---------------------------------------------------------------- montagem
bottom <- pB + pC + plot_layout(widths = c(5, 2.1))
fig <- (pA / bottom) + plot_layout(heights = c(2, 3)) +
  plot_annotation(tag_levels = "A") &                       # tags: A (pA), B (pB), C (pC)
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave(file.path(FIG_DIR, "Figure5.pdf"),  fig, width = W, height = H, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure5.tiff"), fig, width = W, height = H, dpi = 600,
       compression = "lzw", bg = "white")
ggsave(file.path(FIG_DIR, "Figure5.png"),  fig, width = W, height = H, dpi = 300, bg = "white")
cat(sprintf("Figura 5: %.1f x %.1f in (2 colunas) -> %s\n", W, H, FIG_DIR))
