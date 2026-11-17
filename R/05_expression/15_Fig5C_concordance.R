### Concordância TCGA x ICGC do efeito M1 vs M0 dos 14 TFs (painel da Fig 5) ###
# Revisão CSBJ-D-26-01039 (R1.M4, R2.M2, R3.3)
# Lê a tabela suplementar gerada pelo figure5_expression.R (não refaz testes).
# Saídas: fig5_concordancia_ICGC.(png|tiff|pdf) e concordancia_TCGA_ICGC.csv

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

FIG5_DIR <- "~/Kidney/Figures/Figure5"
MM_FILE  <- file.path(FIG5_DIR, "SuppTable_expression_M1_vs_M0.csv")
TN_FILE  <- file.path(FIG5_DIR, "SuppTable_expression_tumor_vs_normal.csv")   # só para o cluster de cada TF
OUT_DIR  <- FIG5_DIR
CL_COL   <- c("Differentiation-associated cluster" = "#1a9850", "Proliferative cluster" = "#542788")

pe <- path.expand; strip <- function(x) sub("\\..*", "", x)
mm <- read.csv(pe(MM_FILE)); tn <- read.csv(pe(TN_FILE))
cl <- tn |> transmute(ensembl = strip(ensembl), cluster) |> distinct(ensembl, .keep_all = TRUE)
cl$cluster[cl$cluster %in% c("Immune cluster", "Immune-related cluster")] <- names(CL_COL)[1]

d <- mm |> mutate(ensembl = strip(ensembl)) |>
  select(ensembl, gene, cohort, cliffs_delta, p_adj) |>
  tidyr::pivot_wider(names_from = cohort, values_from = c(cliffs_delta, p_adj)) |>
  rename(delta_TCGA = `cliffs_delta_TCGA-KIRC`, delta_ICGC = `cliffs_delta_ICGC-RECA`,
         padj_TCGA = `p_adj_TCGA-KIRC`, padj_ICGC = `p_adj_ICGC-RECA`) |>
  left_join(cl, by = "ensembl") |>
  mutate(concordant = sign(delta_TCGA) == sign(delta_ICGC))

ct <- cor.test(d$delta_TCGA, d$delta_ICGC, method = "spearman", exact = FALSE)
by_cl <- d |> group_by(cluster) |> summarise(n = n(), concordant = sum(concordant), .groups = "drop")
write.csv(d, file.path(pe(OUT_DIR), "concordancia_TCGA_ICGC.csv"), row.names = FALSE)

lab <- sprintf("Same direction: %d/%d TFs\nSpearman ρ = %.2f, p = %.3f",
               sum(d$concordant), nrow(d), ct$estimate, ct$p.value)
g <- ggplot(d, aes(delta_TCGA, delta_ICGC, color = cluster)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#f2f1ec", alpha = 0.7) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#f2f1ec", alpha = 0.7) +
  geom_hline(yintercept = 0, color = "#8a8983", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#8a8983", linewidth = 0.4) +
  geom_point(aes(shape = concordant), size = 2.6) +
  geom_text(aes(label = gene), size = 2.6, vjust = -0.9, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.3, size = 2.8, color = "#52514e") +
  scale_color_manual(values = CL_COL, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Concordant", `FALSE` = "Discordant"),
                     name = NULL) +
  labs(x = "Cliff's δ, M1 vs M0 (TCGA-KIRC)", y = "Cliff's δ, M1 vs M0 (ICGC-RECA)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(pe(OUT_DIR), "fig5_concordancia_ICGC.pdf"), g, width = 4.6, height = 4.6)
ggsave(file.path(pe(OUT_DIR), "fig5_concordancia_ICGC.png"), g, width = 4.6, height = 4.6, dpi = 300)
ggsave(file.path(pe(OUT_DIR), "fig5_concordancia_ICGC.tiff"), g, width = 4.6, height = 4.6, dpi = 300,
       compression = "lzw")

cat(sprintf("Mesma direção: %d/%d | Spearman rho = %.2f, p = %.3f\n", sum(d$concordant), nrow(d),
            ct$estimate, ct$p.value))
print(by_cl)
print(d |> transmute(gene, cluster, delta_TCGA = round(delta_TCGA, 2), delta_ICGC = round(delta_ICGC, 2),
                     concordant))
