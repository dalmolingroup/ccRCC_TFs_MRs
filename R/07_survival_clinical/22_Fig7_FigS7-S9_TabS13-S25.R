### Seção 3.5 — Figura 7, figuras e tabelas suplementares ###
# Revisão CSBJ-D-26-01039 (R1.M2, R1.M5, R1.M6, R2.M4, R3.6, R3.7)
#
# NÃO refaz as análises: lê as saídas dos quatro scripts da seção (+ activity_vs_clinical.R) e
# monta as figuras e tabelas. Rode antes, cada um com a linha extra indicada (salva os dados por paciente):
#
#   activity_vs_clinical.R    -> results_activity_clinical/
#       após  df <- inner_join(act, clin, by = "patient"):
#       write.csv(df, file.path(OUT_DIR, "dados_por_paciente.csv"), row.names = FALSE)
#   survival_pfi_dss_act.R    -> results_survival/
#       após  df <- act |> inner_join(clin, ...) |> left_join(meta, ...):
#       write.csv(df, file.path(OUT_DIR, "dados_sobrevivencia.csv"), row.names = FALSE)
#   VHL_influence.R           -> results_vhl/
#       após  wt <- filter(dv, VHL == "Wild-type")   (linha ~153):
#       write.csv(dv, file.path(OUT_DIR, "dados_VHL_por_paciente.csv"), row.names = FALSE)
#   tumor_purity_and_deconvolution.R -> results_purity_deconv/   (sem mudança)
#   survival_ICGC.R (versão com os dois clusters) -> results_survival_ICGC/
#
# Os rótulos antigos ("Immune cluster", immune_cluster, immune) são reconhecidos e trocados.
#
# Saídas (OUT_DIR):
#   Figure7.(tiff|png|pdf)          A atividade x grau e estágio; B forest do Cox; C KM do PFI nos M0
#   FigS_KM_all.(tiff|png|pdf)      KMs de OS, PFI e DSS, todos e M0, dos dois clusters
#   FigS_VHL.(tiff|png|pdf)         atividade por VHL; M1 x M0 nos selvagens; Cox por subgrupo
#   FigS_purity.(tiff|png|pdf)      correlações com pureza/células; Cox e logística com ajuste
#   TableS_*.csv + SupplementaryTables_section3.5.xlsx
#   + no console: números para o texto

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(survival); library(cowplot); library(patchwork)
})
# evita conflito com AnnotationDbi/MASS/biomaRt que estejam carregados na sessão
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename; count <- dplyr::count


## ---------------------------------------------------------------- CONFIG
ACT_DIR  <- "results_activity_clinical"
SURV_DIR <- "results_survival"
VHL_DIR  <- "results_vhl"
PUR_DIR  <- "results_purity_deconv"
ICGC_DIR <- "~/Kidney/results_survival_ICGC"
OUT_DIR  <- "results_section3.5"
ENDPOINTS <- c("OS", "PFI", "DSS")
TIME_DIV  <- 30.44           # tempos do CDR em dias -> meses

CL_LAB <- c(Cluster_1 = "Proliferative cluster", Cluster_2 = "Differentiation-associated cluster")
CL_COL <- c(Cluster_1 = "#542788", Cluster_2 = "#1a9850")
GREY   <- "#9a9893"; INK2 <- "#52514e"

pe <- path.expand
dir.create(pe(OUT_DIR), showWarnings = FALSE, recursive = TRUE)
rd <- function(dir, f, required = FALSE) {
  p <- file.path(pe(dir), f)
  if (file.exists(p)) return(read.csv(p, check.names = FALSE))
  if (required) stop("Arquivo não encontrado: ", p, "\n(veja no cabeçalho a linha a acrescentar no script de origem)")
  message("[aviso] não encontrado: ", p); NULL
}
norm_cl <- function(x) {
  y <- tolower(as.character(x))
  ifelse(grepl("immune|differentiation|cluster[_ ]?2", y), "Cluster_2",
         ifelse(grepl("prolif|cluster[_ ]?1", y), "Cluster_1", NA))
}
score_cols <- function(d) {            # colunas dos escores dos clusters numa tabela por paciente
  cand <- names(d)[grepl("cluster$|cluster[_ ]?[12]$|^immune$|^proliferative$", names(d), ignore.case = TRUE)]
  cand <- cand[!is.na(norm_cl(cand))]
  setNames(cand, norm_cl(cand))[c("Cluster_1", "Cluster_2")]
}
long_scores <- function(d) {
  sc <- score_cols(d); if (anyNA(sc)) stop("Escores dos dois clusters não encontrados: ", paste(names(d), collapse = ", "))
  bind_rows(lapply(names(sc), function(k) d |> mutate(cl = k, score = .data[[sc[[k]]]])))
}
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
fmt_hr <- function(h, l, u) sprintf("%.2f (%.2f–%.2f)", h, l, u)
th <- theme_minimal(base_size = 7) +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(color = INK2),
        strip.text = element_text(face = "bold", size = 7), legend.key.size = unit(3, "mm"))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(pe(OUT_DIR), paste0(name, ".tiff")), p, width = w, height = h, units = "mm", dpi = 500,
         compression = "lzw", bg = "white")
  ggsave(file.path(pe(OUT_DIR), paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(pe(OUT_DIR), paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
}
tabs <- list()
put <- function(d, name) { if (is.null(d) || !nrow(d)) return(invisible()); tabs[[name]] <<- d
write.csv(d, file.path(pe(OUT_DIR), paste0("TableS_", name, ".csv")), row.names = FALSE) }

## ---------------------------------------------------------------- Kaplan-Meier (ggplot, sem survminer)
km_fit <- function(d, endpoint, k) {
  tc <- paste0(endpoint, ".time")
  d <- d[!is.na(d[[endpoint]]) & !is.na(d[[tc]]) & d[[tc]] > 0 & !is.na(d$score), ]
  d$t <- d[[tc]] / TIME_DIV; d$e <- d[[endpoint]]
  d$g <- factor(ifelse(d$score > median(d$score), "High", "Low"), levels = c("High", "Low"))
  fit <- survfit(Surv(t, e) ~ g, data = d)
  st <- rep(sub("g=", "", names(fit$strata)), fit$strata)
  cur <- data.frame(g = st, time = fit$time, surv = fit$surv, cens = fit$n.censor)
  cur <- bind_rows(data.frame(g = c("High", "Low"), time = 0, surv = 1, cens = 0), cur)
  lr <- survdiff(Surv(t, e) ~ g, data = d)
  p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  nn <- table(d$g); ev <- tapply(d$e, d$g, sum)
  labs_g <- setNames(sprintf("%s (n = %d; %d events)", names(nn), as.integer(nn), as.integer(ev)), names(nn))
  list(d = d, fit = fit, cur = cur, p = p, labs = labs_g, col = c(High = unname(CL_COL[k]), Low = GREY))
}
km_plot <- function(km, title = NULL, ylab = "Probability", risk = FALSE, brk = 24) {
  xmax <- max(km$d$t)
  g <- ggplot(km$cur, aes(time, surv, color = g)) +
    geom_step(linewidth = 0.5) +
    geom_point(data = filter(km$cur, cens > 0), shape = 3, size = 0.7, stroke = 0.3) +
    annotate("text", x = 0, y = 0.05, hjust = 0, size = 2.2, color = INK2,
             label = paste0("Log-rank p = ", fmt_p(km$p))) +
    scale_color_manual(values = km$col, labels = km$labs, name = NULL) +
    scale_x_continuous(breaks = seq(0, xmax, brk), limits = c(0, xmax)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Months", y = ylab, title = title) + th +
    theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.text = element_text(size = 5.5), plot.title = element_text(size = 7, face = "bold"))
  if (!risk) return(g)
  tt <- seq(0, xmax, brk)
  s <- summary(km$fit, times = tt, extend = TRUE)
  rt <- data.frame(g = sub("g=", "", as.character(s$strata)), time = s$time, n = s$n.risk)
  r <- ggplot(rt, aes(time, factor(g, c("Low", "High")), label = n, color = g)) + geom_text(size = 2.1) +
    scale_color_manual(values = km$col, guide = "none") +
    scale_x_continuous(breaks = tt, limits = c(0, xmax)) +
    labs(x = NULL, y = NULL, title = "Number at risk") + th +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), plot.title = element_text(size = 6))
  g / r + plot_layout(heights = c(1, 0.22))
}

## ================================================================= dados
act_res <- rd(ACT_DIR, "atividade_vs_clinico.csv")
act_pt  <- rd(ACT_DIR, "dados_por_paciente.csv")
cox     <- rd(SURV_DIR, "cox_resultados.csv", required = TRUE)
cox_full <- rd(SURV_DIR, "cox_cluster_modelo_completo.csv")
flow    <- rd(SURV_DIR, "fluxo_sobrevivencia.csv")
surv_pt <- rd(SURV_DIR, "dados_sobrevivencia.csv", required = TRUE)

cox <- cox |> mutate(cl = norm_cl(variable), cohort_lab = recode(cohort, all = "All patients", M0_only = "M0 patients"))
cl_cox <- cox |> filter(!is.na(cl))

## ================================================================= FIGURA 7
# A: atividade x grau e estágio (todos os pacientes)
if (!is.null(act_pt)) {
  la <- long_scores(act_pt)
  la <- bind_rows(
    la |> filter(!is.na(grade_ord)) |> transmute(cl, score, feature = "Histological grade",
                                                 group = factor(paste0("G", grade_ord), paste0("G", 1:4))),
    if ("stage_ord" %in% names(la)) la |> filter(!is.na(stage_ord)) |>
      transmute(cl, score, feature = "AJCC stage", group = factor(c("I", "II", "III", "IV")[stage_ord], c("I", "II", "III", "IV"))))
  la$group <- factor(as.character(la$group), c(paste0("G", 1:4), "I", "II", "III", "IV"))
  la$feature <- factor(la$feature, c("Histological grade", "AJCC stage"))
  la$cl_lab <- factor(CL_LAB[la$cl], CL_LAB)
  ann <- act_res |> filter(!is.na(norm_cl(variable)), cohort == "all", feature %in% c("Grade", "Stage")) |>
    transmute(cl_lab = factor(CL_LAB[norm_cl(variable)], CL_LAB),
              feature = factor(ifelse(feature == "Grade", "Histological grade", "AJCC stage"), levels(la$feature)),
              txt = sprintf("ρ = %.2f, p = %s", rho, fmt_p(p_spearman)))
  ns <- la |> count(cl_lab, feature, group)
  pA <- ggplot(la, aes(group, score)) +
    geom_jitter(width = 0.18, size = 0.35, alpha = 0.3, color = "#8a8983") +
    geom_boxplot(aes(fill = cl_lab), width = 0.6, outlier.shape = NA, alpha = 0.75, linewidth = 0.3) +
    geom_text(data = ns, aes(group, -Inf, label = n), vjust = -0.5, size = 1.9, color = INK2) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.4, size = 2.2,
              color = INK2, inherit.aes = FALSE) +
    scale_fill_manual(values = setNames(CL_COL, CL_LAB), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.18))) +
    facet_grid(cl_lab ~ feature, scales = "free", space = "free_x") +
    labs(x = NULL, y = "Cluster activity (mean dES)") + th + theme(panel.grid.major.x = element_blank())
} else pA <- ggplot() + theme_void() + labs(title = "dados_por_paciente.csv ausente")

# B: forest do Cox multivariável
fb <- cl_cox |> filter(model == "multivariable") |>
  mutate(cl_lab = factor(CL_LAB[cl], rev(CL_LAB)), endpoint = factor(endpoint, rev(ENDPOINTS)),
         txt = sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high))
pB <- ggplot(fb, aes(HR, endpoint, color = cl_lab)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.5,
                 position = position_dodge(0.6)) +
  geom_point(size = 1.6, position = position_dodge(0.6)) +
  geom_text(aes(x = Inf, label = txt), hjust = 1.02, size = 1.9, position = position_dodge(0.6),
            show.legend = FALSE) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.6))) +
  scale_color_manual(values = setNames(CL_COL, CL_LAB), breaks = CL_LAB, name = NULL) +
  facet_wrap(~ cohort_lab, ncol = 1) +
  labs(x = "HR per SD of activity (log scale)", y = NULL) + th +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

# C: KM do PFI nos M0, dois clusters
m0 <- surv_pt |> filter(M == "M0")
sc <- score_cols(surv_pt)
kmC <- lapply(names(sc), function(k) {
  km <- km_fit(m0 |> mutate(score = .data[[sc[[k]]]]), "PFI", k)
  km_plot(km, title = paste0(CL_LAB[[k]], ", M0"), ylab = "Progression-free interval", risk = TRUE)
})
pC <- wrap_plots(kmC, ncol = 1)

fig7 <- plot_grid(
  pA,
  plot_grid(pB, patchworkGrob(pC), ncol = 2, rel_widths = c(0.9, 1), labels = c("B", "C"), label_size = 10),
  ncol = 1, rel_heights = c(1, 1.15), labels = c("A", ""), label_size = 10)
save_fig(fig7, "Figure7", 190, 215)

## ================================================================= Fig S: todos os KMs
kms <- list()
for (k in names(sc)) for (coh in c("All patients", "M0 patients")) for (e in ENDPOINTS) {
  d <- if (coh == "All patients") surv_pt else m0
  km <- km_fit(d |> mutate(score = .data[[sc[[k]]]]), e, k)
  kms[[length(kms) + 1]] <- km_plot(km, title = sprintf("%s — %s, %s", CL_LAB[[k]], e, coh), ylab = e)
}
save_fig(wrap_plots(kms, ncol = 3) + plot_annotation(tag_levels = "A"), "FigS_KM_all", 190, 250)

## ================================================================= Tabelas: Cox, fluxo, clínico
bh_cl <- cl_cox |> group_by(model) |> mutate(p_BH = p.adjust(p, "BH")) |> ungroup()
# pacientes com o desfecho disponível em cada coorte -> excluídos por covariável ausente = disponíveis - n do modelo
avail <- bind_rows(lapply(c("all", "M0_only"), function(coh) {
  d <- if (coh == "all") surv_pt else dplyr::filter(surv_pt, M == "M0")
  data.frame(cohort = coh, endpoint = ENDPOINTS, available = sapply(ENDPOINTS, function(e) {
    tc <- paste0(e, ".time"); sum(!is.na(d[[e]]) & !is.na(d[[tc]]) & d[[tc]] > 0) }))
}))
bh_cl <- bh_cl |> left_join(avail, by = c("cohort", "endpoint")) |> mutate(excluded = available - n)
put(bh_cl |> arrange(cl, model, cohort, factor(endpoint, ENDPOINTS)) |>
      transmute(Cluster = CL_LAB[cl], Cohort = cohort_lab, Endpoint = endpoint,
                Model = ifelse(model == "multivariable", "Multivariable (age, sex, stage, grade)", "Univariable"),
                `Patients with endpoint` = available, `Excluded (missing covariates)` = excluded,
                n, Events = events, `HR per SD (95% CI)` = fmt_hr(HR, CI_low, CI_high), p = fmt_p(p),
                `BH-adjusted p (12 cluster tests per model)` = fmt_p(p_BH),
                `Schoenfeld p (score)` = fmt_p(PH_p_variable), `Schoenfeld p (global)` = fmt_p(PH_p_global)),
    "cox_clusters")
if (!is.null(cox_full)) put(cox_full |> mutate(Cluster = CL_LAB[norm_cl(cluster)]) |>
                              transmute(Cluster, Cohort = recode(cohort, all = "All patients", M0_only = "M0 patients"), Endpoint = endpoint,
                                        Term = sub("score_SD", "Cluster activity (per SD)", term), n, Events = events,
                                        `HR (95% CI)` = fmt_hr(HR, CI_low, CI_high), p = fmt_p(p)), "cox_full_models")
if (!is.null(flow)) put(flow, "survival_flow")
tf_cox <- cox |> filter(is.na(cl), model == "multivariable")
if (nrow(tf_cox)) put(tf_cox |> transmute(Regulator = label, Cohort = cohort_lab, Endpoint = endpoint, n, Events = events,
                                          `HR per SD (95% CI)` = fmt_hr(HR, CI_low, CI_high), p = fmt_p(p),
                                          `BH-adjusted p (14 TFs)` = fmt_p(FDR)), "cox_individual_TFs")
if (!is.null(act_res)) put(act_res |>
                             mutate(Variable = ifelse(!is.na(norm_cl(variable)), CL_LAB[norm_cl(variable)], label)) |>
                             transmute(Variable, Cohort = recode(cohort, all = "All patients", M0_only = "M0 patients (stage IV excluded)"),
                                       Feature = feature, n, `Spearman ρ` = round(rho, 2), `p (Spearman)` = fmt_p(p_spearman),
                                       `p (Kruskal–Wallis)` = fmt_p(p_kruskal), `BH-adjusted p (TFs)` = fmt_p(FDR_spearman)),
                           "activity_clinical")

## ================================================================= Fig S + Tabela: VHL
vt <- list(act = rd(VHL_DIR, "act_tab_dois_clusters.csv"), cor = rd(VHL_DIR, "cor_tab_dois_clusters.csv"),
           cox = rd(VHL_DIR, "cox_tab_dois_clusters.csv"), int = rd(VHL_DIR, "inter_tab_dois_clusters.csv"))
vpt <- rd(VHL_DIR, "dados_VHL_por_paciente.csv")
if (!is.null(vt$act)) {
  vt <- lapply(vt, function(d) if (!is.null(d)) mutate(d, cl = norm_cl(cluster), Cluster = CL_LAB[cl]) else NULL)
  put(vt$act |> transmute(Cluster, Comparison = comparison, `Group 1` = sprintf("%s (n = %d)", group1, n1),
                          `Group 2` = sprintf("%s (n = %d)", group2, n2), `Cliff's δ` = round(cliffs_delta, 2),
                          p = fmt_p(p)), "VHL_activity")
  if (!is.null(vt$cor)) put(vt$cor |> transmute(Cluster, VHL, Feature = feature, n, `Spearman ρ` = round(rho, 2),
                                                p = fmt_p(p)), "VHL_grade_stage")
  if (!is.null(vt$cox)) put(vt$cox |> left_join(vt$int[, c("cl", "cohort", "endpoint", "p_interaction")],
                                                by = c("cl", "cohort", "endpoint")) |>
                              transmute(Cluster, Cohort = recode(cohort, all = "All patients", M0_only = "M0 patients"),
                                        `VHL status` = subgroup, Endpoint = endpoint, n, Events = events,
                                        `HR per SD (95% CI)` = fmt_hr(HR, CI_low, CI_high), p = fmt_p(p),
                                        `p interaction (cluster × VHL)` = fmt_p(p_interaction)), "VHL_cox")
  panels <- list()
  if (!is.null(vpt)) {
    lv <- long_scores(vpt) |> filter(!is.na(VHL)) |> mutate(cl_lab = factor(CL_LAB[cl], CL_LAB))
    lab_d <- function(cmp) vt$act |> filter(comparison == cmp) |>
      transmute(cl_lab = factor(Cluster, CL_LAB), txt = sprintf("δ = %.2f, p = %s", cliffs_delta, fmt_p(p)))
    box <- function(d, x, ann, fills) ggplot(d, aes(.data[[x]], score)) +
      geom_jitter(width = 0.18, size = 0.35, alpha = 0.3, color = "#8a8983") +
      geom_boxplot(aes(fill = .data[[x]]), width = 0.55, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
      geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.4, size = 2.2,
                color = INK2, inherit.aes = FALSE) +
      scale_fill_manual(values = fills, guide = "none") +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
      facet_wrap(~ cl_lab, scales = "free_y", labeller = label_wrap_gen(22)) + labs(x = NULL, y = "Cluster activity (mean dES)") + th
    nV <- table(vpt$VHL)
    panels$A <- box(lv |> mutate(VHL = factor(VHL, c("Wild-type", "Mutated"))), "VHL", lab_d("VHL (all patients)"),
                    c("Wild-type" = "#86b6ef", Mutated = "#1c5cab")) +
      scale_x_discrete(labels = function(x) sprintf("%s\n(n = %d)", x, as.integer(nV[x])))
    lw <- lv |> filter(VHL == "Wild-type", M %in% c("M0", "M1"))
    nM <- table(lw$M[lw$cl == "Cluster_1"])
    panels$B <- box(lw, "M", lab_d("M1 vs M0 (VHL Wild-type)"), c(M0 = "#2166ac", M1 = "#d73027")) +
      scale_x_discrete(labels = function(x) sprintf("%s\n(n = %d)", x, as.integer(nM[x]))) +
      labs(subtitle = "VHL wild-type tumours")
  }
  if (!is.null(vt$cox)) {
    fv <- vt$cox |> left_join(vt$int[, c("cl", "cohort", "endpoint", "p_interaction")], by = c("cl", "cohort", "endpoint")) |>
      mutate(cl_lab = factor(Cluster, CL_LAB), endpoint = factor(endpoint, ENDPOINTS),
             cohort = recode(cohort, all = "All", M0_only = "M0"),
             row = factor(paste(endpoint, cohort), rev(paste(rep(ENDPOINTS, each = 2), c("All", "M0")))),
             subgroup = factor(subgroup, c("Wild-type", "Mutated")))
    pint <- fv |> distinct(cl_lab, row, p_interaction) |> mutate(txt = paste0("p int = ", sprintf("%.2f", p_interaction)))
    panels$C <- ggplot(fv, aes(HR, row, color = subgroup)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.45, position = position_dodge(0.6)) +
      geom_point(size = 1.3, position = position_dodge(0.6)) +
      geom_text(data = pint, aes(x = Inf, y = row, label = txt), hjust = 1.05, size = 1.9, color = INK2,
                inherit.aes = FALSE) +
      scale_x_log10(expand = expansion(mult = c(0.05, 0.45))) +
      scale_color_manual(values = c("Wild-type" = "#2a78d6", Mutated = "#eb6834"), name = "VHL") +
      facet_wrap(~ cl_lab) + labs(x = "HR per SD (log scale), adjusted for age, sex, stage and grade", y = NULL) +
      th + theme(legend.position = "bottom", panel.grid.major.y = element_blank())
  }
  if (length(panels)) {
    lay <- if (length(panels) == 3) (panels$A | panels$B) / panels$C + plot_layout(heights = c(1, 1.2)) else wrap_plots(panels, ncol = 1)
    save_fig(lay + plot_annotation(tag_levels = "A"), "FigS_VHL", 190, 170)
  }
}

## ================================================================= Fig S + Tabela: pureza e deconvolução
cors <- rd(PUR_DIR, "correlacoes_atividade_vs_TME.csv")
logit <- rd(PUR_DIR, "M1_vs_M0_logistica_ajustada.csv")
coxp <- rd(PUR_DIR, "cox_ajustado_pureza.csv")
ADJ_LAB <- c(none = "None", purity = "Purity", immune_score = "ImmuneScore", stromal_score = "StromalScore",
             cpe = "CPE purity", purity_immune = "Purity + ImmuneScore")
if (!is.null(cors)) {
  cc <- cors |> filter(!is.na(norm_cl(variable)), grepl("cluster$|cluster[_ ]?[12]$|^immune$|^proliferative$", variable)) |>
    mutate(cl = norm_cl(variable), Cluster = CL_LAB[cl])
  put(cc |> transmute(Cluster, Method = method, Feature = cell_label, `Spearman ρ` = round(rho, 2), p = fmt_p(p),
                      `BH-adjusted p` = fmt_p(FDR)) |> arrange(Method, Feature, Cluster), "purity_correlations")
  # concordância entre os métodos: mesma categoria celular no quanTIseq e no CIBERSORT
  cat_of <- function(x) dplyr::case_when(
    grepl("CD8", x) ~ "CD8 T cells", grepl("regulatory|Treg", x, ignore.case = TRUE) ~ "Tregs",
    grepl("CD4", x) ~ "CD4 T cells", grepl("Dendritic|DC", x) ~ "Dendritic cells", grepl("NK", x) ~ "NK cells",
    grepl("Neutrophil", x) ~ "Neutrophils", grepl("Monocyte", x) ~ "Monocytes",
    grepl("Macrophage", x) ~ "Macrophages", grepl("B cell|Plasma", x) ~ "B cells", TRUE ~ NA_character_)
  conc <- cc |> filter(method %in% c("quanTIseq", "CIBERSORT")) |> mutate(category = cat_of(cell_label)) |>
    filter(!is.na(category), !is.na(rho)) |> group_by(Cluster, category, method) |>
    summarise(cells = paste(sprintf("%s (%.2f%s)", cell_label, rho, ifelse(!is.na(FDR) & FDR < 0.05, "*", "")), collapse = "; "),
              sign = { s <- sign(rho[!is.na(FDR) & FDR < 0.05]); if (length(s) && all(s == s[1])) s[1] else 0 }, .groups = "drop") |>
    tidyr::pivot_wider(names_from = method, values_from = c(cells, sign)) |>
    mutate(Concordant = !is.na(sign_quanTIseq) & !is.na(sign_CIBERSORT) & sign_quanTIseq != 0 &
             sign_quanTIseq == sign_CIBERSORT) |>
    transmute(Cluster, `Cell type` = category, `quanTIseq (ρ)` = cells_quanTIseq, `CIBERSORT (ρ)` = cells_CIBERSORT,
              `Concordant between methods` = ifelse(is.na(cells_quanTIseq) | is.na(cells_CIBERSORT), "Only one method",
                                                    ifelse(Concordant, "Yes", "No"))) |> arrange(Cluster, `Cell type`)
  put(conc, "purity_method_concordance")
  hm <- cc |> filter(!is.na(rho)) |>
    mutate(method = factor(sub("ESTIMATE/CPE", "ESTIMATE", method), c("ESTIMATE", "quanTIseq", "CIBERSORT")),
           Cluster = factor(sub(" cluster", "", Cluster), sub(" cluster", "", CL_LAB)),
           star = ifelse(!is.na(FDR) & FDR < 0.05, "*", ""))
  pP <- ggplot(hm, aes(Cluster, cell_label, fill = rho)) + geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f%s", rho, star)), size = 1.9) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = expression(rho)) +
    facet_grid(method ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL, caption = "* BH-adjusted p < 0.05") + th +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0))
}
if (!is.null(logit)) {
  lg <- logit |> mutate(cl = norm_cl(score), Cluster = CL_LAB[cl], adj = ADJ_LAB[adjusted_for])
  put(lg |> transmute(Cluster, `Additional adjustment` = adj, n, M1, `OR per SD (95% CI)` = fmt_hr(OR_per_SD, CI_low, CI_high),
                      p = fmt_p(p)), "purity_logistic_M1")
  pL <- ggplot(lg, aes(OR_per_SD, factor(adj, rev(ADJ_LAB)), color = Cluster)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.45) + geom_point(size = 1.3) +
    scale_x_log10() + scale_color_manual(values = setNames(CL_COL, CL_LAB), guide = "none") +
    facet_wrap(~ factor(Cluster, CL_LAB), scales = "free_x", labeller = label_wrap_gen(22)) +
    labs(x = "OR per SD for M1 (log scale)", y = "Adjustment") + th
}
if (!is.null(coxp)) {
  cp <- coxp |> mutate(cl = norm_cl(score), Cluster = CL_LAB[cl], adj = ADJ_LAB[adjusted_for],
                       cohort = recode(cohort, all = "All patients", M0_only = "M0 patients"))
  put(cp |> transmute(Cluster, Cohort = cohort, Endpoint = endpoint, `Additional adjustment` = adj, n, Events = events,
                      `HR per SD (95% CI)` = fmt_hr(HR, CI_low, CI_high), p = fmt_p(p)), "purity_cox")
  pX <- ggplot(cp |> mutate(endpoint = factor(endpoint, ENDPOINTS)),
               aes(HR, factor(adj, rev(ADJ_LAB)), color = factor(Cluster, CL_LAB))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.45, position = position_dodge(0.6)) +
    geom_point(size = 1.2, position = position_dodge(0.6)) + scale_x_log10() +
    scale_color_manual(values = setNames(CL_COL, CL_LAB), name = NULL) +
    facet_grid(cohort ~ endpoint) +
    labs(x = "HR per SD (log scale), adjusted for age, sex, stage, grade + listed variable", y = "Adjustment") +
    th + theme(legend.position = "bottom", panel.border = element_rect(fill = NA, color = "grey85"))
}
if (exists("pP") && exists("pX")) {
  right <- if (exists("pL")) (pL / pX) + plot_layout(heights = c(0.55, 1)) else pX
  save_fig(wrap_elements(full = pP) + wrap_elements(full = right) + plot_layout(widths = c(0.75, 1)) +
             plot_annotation(tag_levels = "A"), "FigS_purity", 190, 200)
}

## ================================================================= Tabela: ICGC
ic <- rd(ICGC_DIR, "cox_ICGC_clusters.csv")
if (is.null(ic)) { ic <- rd(ICGC_DIR, "cox_ICGC_cluster.csv"); if (!is.null(ic) && !"cluster" %in% names(ic)) {
  ic$cluster <- "Differentiation-associated cluster"; message("[aviso] ICGC só com o cluster de diferenciação") } }
if (!is.null(ic)) put(ic |> mutate(Cluster = CL_LAB[norm_cl(cluster)]) |>
                        transmute(Cluster, Cohort = cohort, Model = model, n, Deaths = events,
                                  `HR per SD (95% CI)` = fmt_hr(HR, CI_low, CI_high), p = fmt_p(p),
                                  `Schoenfeld p` = fmt_p(PH_p), `Events per variable` = round(events_per_variable, 1)), "ICGC_OS")

if (requireNamespace("openxlsx", quietly = TRUE) && length(tabs))
  openxlsx::write.xlsx(setNames(tabs, substr(names(tabs), 1, 31)), file.path(pe(OUT_DIR), "SupplementaryTables_section3.5.xlsx"), firstRow = TRUE, colWidths = "auto")

## ================================================================= números para o texto
cat("\n================ NÚMEROS PARA A SEÇÃO 3.5 ================\n")
if (!is.null(act_res)) {
  cat("\nAtividade x clínico (Spearman):\n")
  print(act_res |> filter(!is.na(norm_cl(variable))) |>
          transmute(cluster = CL_LAB[norm_cl(variable)], cohort, feature, n, rho = round(rho, 2), p = fmt_p(p_spearman)),
        row.names = FALSE)
}
cat("\nCox multivariável dos clusters (HR por DP; BH entre os 12 testes):\n")
print(as.data.frame(bh_cl |> filter(model == "multivariable") |>
                      transmute(cluster = CL_LAB[cl], cohort, endpoint, n, events, HR = fmt_hr(HR, CI_low, CI_high),
                                p = fmt_p(p), p_BH = fmt_p(p_BH), PH = fmt_p(PH_p_variable))), row.names = FALSE)
cat(sprintf("\nMenor p de Schoenfeld (escores): %s\n", fmt_p(min(bh_cl$PH_p_variable, na.rm = TRUE))))
fu <- sapply(ENDPOINTS, function(e) {           # seguimento mediano (KM reverso)
  tc <- paste0(e, ".time"); d <- surv_pt[!is.na(surv_pt[[e]]) & surv_pt[[tc]] > 0, ]
  s <- summary(survfit(Surv(d[[tc]] / TIME_DIV, 1 - d[[e]]) ~ 1))$table["median"]; round(unname(s), 1)
})
cat("Seguimento mediano (meses, KM reverso):", paste(ENDPOINTS, fu, collapse = "; "), "\n")
if (!is.null(ic)) { cat("\nICGC (OS):\n"); print(as.data.frame(tabs$ICGC_OS), row.names = FALSE) }
cat("\nSaídas em", pe(OUT_DIR), "\n")