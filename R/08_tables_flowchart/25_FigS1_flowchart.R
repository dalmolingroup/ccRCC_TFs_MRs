### Figura suplementar S1: fluxograma dos pacientes (TCGA-KIRC e ICGC-RECA) ###
# Revisão CSBJ-D-26-01039 (R2 comentário adicional 1, R3.2)
# Só ggplot2. Edite as contagens no bloco N; o que estiver como NA aparece como "[n]" na figura.
# Saídas: FigS1_flowchart.(tiff|png|pdf)

suppressPackageStartupMessages(library(ggplot2))

OUT_DIR <- "results_supp_figures"

## ---------------------------------------------------------------- contagens
N <- list(
  tcga_samples = 602, tcga_patients = 537,
  normal = 72, primary = 530,
  mx = 25, network = 505, duplicates = 4,          # confira: 505 - 501 = 4 alíquotas duplicadas
  analysed = 501, m0 = 422, m1 = 79,
  train = 350, train_m1 = 55, test = 151, test_m1 = 24,
  clinical = 494, surv_min = 481, surv_max = 492,
  vhl = 340, vhl_mut = 163, vhl_wt = 177, vhl_na = 161,
  deconv = 501,                                      # 'with_deconvolution' do results_purity_deconv/fluxo.csv
  icgc = 91, icgc_m_known = 90, icgc_m0 = 81, icgc_m1 = 9, icgc_normal = 45, icgc_deaths = 30
)
f <- function(x) if (is.na(x)) "[n]" else format(x, big.mark = ",")

## ---------------------------------------------------------------- caixas
# x, y = centro; w, h = largura e altura (unidades arbitrárias); type = main | excl | split | use | icgc
B <- function(id, x, y, w, h, label, type) data.frame(id, x, y, w, h, label, type)
boxes <- rbind(
  B("t1", 3.0, 13.0, 3.6, 1.0, sprintf("TCGA-KIRC RNA-seq (GDC)\n%s samples from %s patients", f(N$tcga_samples), f(N$tcga_patients)), "main"),
  B("e1", 7.0, 13.0, 3.0, 1.0, sprintf("Solid-tissue normal samples: %s\n(tumour-normal comparison only)", f(N$normal)), "excl"),
  B("t2", 3.0, 11.2, 3.6, 0.8, sprintf("Primary tumour samples: %s", f(N$primary)), "main"),
  B("e2", 7.0, 11.2, 3.0, 0.8, sprintf("Unknown metastatic status\n(MX or missing): %s", f(N$mx)), "excl"),
  B("t3", 3.0, 9.4, 3.6, 1.0, sprintf("Primary tumour samples: %s\nRegulatory network inference (RTN)", f(N$network)), "main"),
  B("e3", 7.0, 9.4, 3.0, 0.8, sprintf("Duplicate aliquots collapsed: %s\n(one sample per patient)", f(N$duplicates)), "excl"),
  B("t4", 3.0, 7.6, 3.6, 1.0, sprintf("Patients analysed: %s\n(%s M0, %s M1)", f(N$analysed), f(N$m0), f(N$m1)), "main"),
  B("s1", 1.6, 5.8, 2.4, 1.0, sprintf("Training set (70%%)\nn = %s (%s M1)", f(N$train), f(N$train_m1)), "split"),
  B("s2", 4.3, 5.8, 2.4, 1.0, sprintf("Test set (30%%)\nn = %s (%s M1)", f(N$test), f(N$test_m1)), "split"),
  B("s3", 7.0, 5.8, 2.4, 1.0, sprintf("Nested cross-validation\n5 outer folds, n = %s", f(N$analysed)), "split"),
  B("u1", 1.15, 3.7, 2.1, 1.3, sprintf("Clinical comparison\n(grade, pT, age)\nn = %s", f(N$clinical)), "use"),
  B("u2", 3.45, 3.7, 2.1, 1.3, sprintf("Survival (TCGA-CDR)\nOS, PFI, DSS\nn = %s-%s", f(N$surv_min), f(N$surv_max)), "use"),
  B("u3", 5.75, 3.7, 2.1, 1.3, sprintf("VHL subgroups\nn = %s (%s mutated,\n%s wild-type)", f(N$vhl), f(N$vhl_mut), f(N$vhl_wt)), "use"),
  B("u4", 8.05, 3.7, 2.1, 1.3, sprintf("Purity and immune\ndeconvolution\nn = %s", f(N$deconv)), "use"),
  B("i1", 11.6, 13.0, 3.4, 1.3, sprintf("ICGC RECA-EU (ccRCC)\n%s primary tumours\n%s matched normal samples", f(N$icgc), f(N$icgc_normal)), "icgc"),
  B("i2", 11.6, 10.7, 3.4, 1.3, sprintf("Tumour-normal differential\nexpression (%s pairs):\nsignature for MRA", f(N$icgc_normal)), "icgc"),
  B("i3", 11.6, 8.2, 3.4, 1.5, sprintf("Locked signature, regulon\nactivity and expression\n%s tumours with known M stage\n(%s M0, %s M1)", f(N$icgc_m_known), f(N$icgc_m0), f(N$icgc_m1)), "icgc"),
  B("i4", 11.6, 5.8, 3.4, 1.0, sprintf("Overall survival (exploratory)\n%s patients, %s deaths", f(N$icgc), f(N$icgc_deaths)), "icgc")
)
FILL <- c(main = "#e8eef7", excl = "#f6f5f2", split = "#ece6f4", use = "#e7f2ea", icgc = "#fdf1e4")
LINE <- c(main = "#2a5c9a", excl = "#8a8983", split = "#542788", use = "#1a9850", icgc = "#c46a1b")

## ---------------------------------------------------------------- setas
b <- function(id) boxes[boxes$id == id, ]
v_arrow <- function(a, z) { A <- b(a); Z <- b(z); data.frame(x = A$x, y = A$y - A$h / 2, xend = Z$x, yend = Z$y + Z$h / 2) }
h_arrow <- function(a, z) { A <- b(a); Z <- b(z); data.frame(x = A$x + A$w / 2, y = A$y, xend = Z$x - Z$w / 2, yend = Z$y) }
arrows <- rbind(v_arrow("t1", "t2"), v_arrow("t2", "t3"), v_arrow("t3", "t4"),
                h_arrow("t1", "e1"), h_arrow("t2", "e2"), h_arrow("t3", "e3"),
                v_arrow("i1", "i2"), v_arrow("i2", "i3"), v_arrow("i3", "i4"))
# t4 -> conjuntos e conjuntos -> análises (barramento)
t4 <- b("t4"); bus1 <- 6.75
lines <- rbind(data.frame(x = t4$x, y = t4$y - t4$h / 2, xend = t4$x, yend = bus1),
               data.frame(x = 1.6, y = bus1, xend = 7.0, yend = bus1))
for (s in c("s1", "s2", "s3")) arrows <- rbind(arrows, data.frame(x = b(s)$x, y = bus1, xend = b(s)$x, yend = b(s)$y + b(s)$h / 2))
bus2 <- 4.85
lines <- rbind(lines, data.frame(x = 0.1, y = 5.8, xend = 0.1, yend = bus2),            # da coorte de 501 para as análises
               data.frame(x = 0.1, y = bus2, xend = 8.05, yend = bus2))
lines <- rbind(lines, data.frame(x = t4$x - t4$w / 2, y = t4$y, xend = 0.1, yend = t4$y),
               data.frame(x = 0.1, y = t4$y, xend = 0.1, yend = 5.8))
for (u in c("u1", "u2", "u3", "u4")) arrows <- rbind(arrows, data.frame(x = b(u)$x, y = bus2, xend = b(u)$x, yend = b(u)$y + b(u)$h / 2))
# ICGC -> MRA na rede do TCGA (tracejado, em cotovelo entre as caixas de exclusão)
i2 <- b("i2"); t3 <- b("t3"); ymra <- 10.3; xdrop <- t3$x + 1.2
mra_l <- data.frame(x = i2$x - i2$w / 2, y = ymra, xend = xdrop, yend = ymra)
mra_a <- data.frame(x = xdrop, y = ymra, xend = xdrop, yend = t3$y + t3$h / 2)
## ---------------------------------------------------------------- figura
headers <- data.frame(x = c(4.4, 11.6), y = 14.15, label = c("TCGA-KIRC (discovery)", "ICGC-RECA"))
g <- ggplot() +
  geom_segment(data = lines, aes(x, y, xend = xend, yend = yend), linewidth = 0.35, color = "#52514e") +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend), linewidth = 0.35, color = "#52514e",
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_segment(data = mra_l, aes(x, y, xend = xend, yend = yend), linewidth = 0.35, color = "#c46a1b", linetype = "dashed") +
  geom_segment(data = mra_a, aes(x, y, xend = xend, yend = yend), linewidth = 0.35, color = "#c46a1b",
               linetype = "dashed", arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("text", x = 7.0, y = ymra + 0.22, label = "Master regulator analysis: 356 MRs",
           size = 2.2, color = "#c46a1b") +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                              fill = type, color = type), linewidth = 0.4) +
  geom_text(data = boxes, aes(x, y, label = label), size = 2.35, lineheight = 0.95, color = "#1f1e1c") +
  geom_text(data = headers, aes(x, y, label = label), size = 3, fontface = "bold", color = "#1f1e1c") +
  scale_fill_manual(values = FILL, guide = "none") + scale_color_manual(values = LINE, guide = "none") +
  coord_cartesian(xlim = c(-0.1, 13.4), ylim = c(2.9, 14.5), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
W <- 190; H <- 165
ggsave(file.path(OUT_DIR, "FigS1_flowchart.tiff"), g, width = W, height = H, units = "mm", dpi = 500,
       compression = "lzw", bg = "white")
ggsave(file.path(OUT_DIR, "FigS1_flowchart.png"), g, width = W, height = H, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "FigS1_flowchart.pdf"), g, width = W, height = H, units = "mm", device = cairo_pdf)
if (any(vapply(N, is.na, TRUE))) message("[aviso] contagens em aberto: ", paste(names(N)[vapply(N, is.na, TRUE)], collapse = ", "))
cat("Salvo em", OUT_DIR, "\n")