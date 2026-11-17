"""
Figuras do pipeline revisado (classificação M1 vs M0, TCGA-KIRC)

Lê os arquivos que o ml_pipeline_revision.py salva em OUT_DIR, então pode ser
rodado sozinho sobre resultados já existentes, sem refazer o pipeline:

    python plots_ml_revision.py                     # usa BALANCING do pipeline
    python plots_ml_revision.py undersample         # outro modo
    python plots_ml_revision.py compare             # class_weight vs undersample

Figuras (PNG 300 dpi + PDF vetorial em <OUT_DIR>/figuras/):
  fig_roc_pr_holdout      ROC e PR no teste intocado: selecionados vs 356 TFs
                          (nova Fig 3; R2 Fig 3, R2.M1, R3.1)
  fig_matriz_confusao     matrizes de confusão no limiar escolhido no treino
  fig_nested_cv           AUC por fold + curva ROC das predições out-of-fold
  fig_estabilidade        frequência de seleção do Boruta (R2.adic, R3.11)
  fig_importancia         importância das features do modelo final
  fig_balanceamento       class weighting vs undersampling (R2.adic)
  Figure3                 painel da Figura 3 do artigo (A ROC/PR, B confusão, C nested CV, D estabilidade);
                          TIFF (LZW), PNG e PDF a 300 dpi

Símbolos dos genes: opcional, um CSV com colunas `ensembl,symbol`
(ex.: gerado com biomaRt). Sem ele, os eixos mostram os IDs Ensembl.
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (average_precision_score, confusion_matrix,
                             precision_recall_curve, roc_auc_score, roc_curve)

# ---------------------------------------------------------------- config
OUT_DIR = "results_ml_revision"
SYMBOL_FILE = "0_data/ensembl_to_symbol.csv"     # opcional: colunas ensembl,symbol
SIGNATURE_FILE = "0_data/signature_class_weight.txt"   # assinatura única do artigo (14 TFs, class weighting)
STABLE_CUTOFF = 0.8
TOP_N_STABILITY = 25

# paleta categórica validada (slot 1 azul, slot 2 laranja) + tinta de texto
C_SEL = "#2a78d6"        # TFs selecionados / class weighting
C_ALL = "#eb6834"        # 356 TFs / undersampling
C_MUTED = "#b8b6af"      # barras sem destaque
C_REF = "#8a8983"        # linhas de referência (acaso, prevalência, corte)
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e6e5e0"
SEQ = ["#ffffff", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # azul sequencial

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,          # texto editável no Illustrator/Inkscape
    "savefig.dpi": 300, "savefig.bbox": "tight",
    "figure.constrained_layout.use": True,
})


# ---------------------------------------------------------------- utilitários
def _path(out_dir, base, balancing, ext):
    """Arquivo com sufixo do modo (novo) ou sem sufixo (execuções antigas)."""
    for name in (f"{base}_{balancing}.{ext}", f"{base}.{ext}"):
        p = os.path.join(out_dir, name)
        if os.path.exists(p):
            return p
    return None


def load_symbols():
    if not os.path.exists(SYMBOL_FILE):
        return {}
    s = pd.read_csv(SYMBOL_FILE)
    s.columns = [c.lower() for c in s.columns]
    s = s.dropna(subset=["symbol"])
    return dict(zip(s["ensembl"].str.split(".").str[0], s["symbol"]))


def label(gene, symbols):
    return symbols.get(str(gene).split(".")[0], gene)


def _norm(g):
    """ENSG sem versão e em maiúsculas: compara IDs de fontes diferentes."""
    return str(g).strip().strip('"').split(".")[0].upper()


def load_signature(out_dir=OUT_DIR, balancing="class_weight"):
    """Assinatura = TFs confirmados pelo Boruta no treino (JSON do modo principal).
    O arquivo SIGNATURE_FILE só é usado se o JSON não existir."""
    p = _path(out_dir, "resultados", balancing, "json")
    if p and balancing in os.path.basename(p):
        conf = json.load(open(p)).get("holdout", {}).get("boruta_confirmed", [])
        if conf:
            return list(conf)
    if not os.path.exists(SIGNATURE_FILE):
        return []
    lines = [l.strip().strip('"') for l in open(SIGNATURE_FILE) if l.strip()]
    return [l for l in lines if l.upper().startswith("ENSG") or l.isalnum()]


def save(fig, fig_dir, name):
    os.makedirs(fig_dir, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(fig_dir, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  salvo: {name}.png / .pdf")


def p_txt(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def ci_txt(m):
    lo, hi = m["AUC_ci95"]
    return f"AUC {m['AUC']:.2f} [{lo:.2f}–{hi:.2f}]"


# ---------------------------------------------------------------- figuras
def fig_roc_pr(res, pred, fig_dir, tag):
    h = res["holdout"]
    y = pred["y"].values
    models = [("p_selecionado", "modelo_selecionado", C_SEL,
               f"Selected TFs (n = {len(h['boruta_confirmed'])})"),
              ("p_356", "modelo_356_TFs", C_ALL, "All 356 TFs")]
    prev = y.mean()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for col, key, color, name in models:
        m = h[key]
        fpr, tpr, _ = roc_curve(y, pred[col])
        a1.plot(fpr, tpr, color=color, lw=2, label=f"{name}\n{ci_txt(m)}")
        # ponto de operação no limiar escolhido no treino
        a1.plot(1 - m["specificity"], m["sensitivity"], "o", ms=8, color=color,
                mec="white", mew=2, zorder=5)
        pr, rc, _ = precision_recall_curve(y, pred[col])
        a2.plot(rc, pr, color=color, lw=2, label=f"{name}\nPR-AUC {m['PR_AUC']:.2f}")
    a1.plot([0, 1], [0, 1], ls="--", lw=1, color=C_REF)
    a2.axhline(prev, ls="--", lw=1, color=C_REF)
    a2.text(0.98, prev - 0.02, f"prevalence = {prev:.2f}", ha="right", va="top",
            color=INK_2, fontsize=7.5)

    d = h["diferenca_AUC_selecionado_vs_356"]
    a1.set(xlabel="1 − specificity", ylabel="Sensitivity", xlim=(-0.01, 1.01), ylim=(-0.01, 1.01),
           title="ROC, untouched test set")
    a1.text(0.97, 0.04, f"ΔAUC {d['diff']:+.2f} [{d['ci95'][0]:+.2f}; {d['ci95'][1]:+.2f}]\n"
            f"paired bootstrap {p_txt(d['p_bootstrap'])}",
            ha="right", va="bottom", color=INK_2, fontsize=7.5)
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1 if len(models) > 2 else 2)
    a2.set(xlabel="Recall (sensitivity)", ylabel="Precision", xlim=(-0.01, 1.01), ylim=(0, 1.01),
           title="Precision–recall, untouched test set")
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1 if len(models) > 2 else 2)
    for a in (a1, a2):
        a.set_aspect("equal")
    fig.text(0.5, -0.12 if len(models) == 2 else -0.24, f"n = {len(y)} ({int(y.sum())} M1). Dots: operating point at the "
             "threshold chosen on training data (Youden, out-of-fold).",
             ha="center", color=INK_2, fontsize=7.5)
    save(fig, fig_dir, f"fig_roc_pr_holdout_{tag}")


def fig_confusion(res, fig_dir, tag):
    h = res["holdout"]
    panels = [("modelo_selecionado", f"Selected TFs (n = {len(h['boruta_confirmed'])})"),
              ("modelo_356_TFs", "All 356 TFs")]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0))
    for ax, (key, title) in zip(axes, panels):
        m = h[key]
        cm = m["confusion_matrix"]
        mat = np.array([[cm["TN"], cm["FP"]], [cm["FN"], cm["TP"]]])
        rows = mat / mat.sum(axis=1, keepdims=True)     # % por classe verdadeira
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ)
        ax.imshow(rows, cmap=cmap, vmin=0, vmax=1)
        for i in range(2):
            for j in range(2):
                dark = rows[i, j] > 0.55
                ax.text(j, i, f"{mat[i, j]}\n({rows[i, j]:.0%})", ha="center", va="center",
                        color="white" if dark else INK, fontsize=9)
        ax.set_xticks([0, 1], ["M0", "M1"]); ax.set_yticks([0, 1], ["M0", "M1"])
        ax.set(xlabel="Predicted", ylabel="Observed")
        ax.set_title(f"{title}\nsens {m['sensitivity']:.2f} · spec {m['specificity']:.2f} · "
                     f"F1 {m['F1']:.2f}", fontsize=8.5)
        ax.grid(False)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.text(0.5, -0.07, "Threshold chosen on training data; percentages are per observed class.",
             ha="center", color=INK_2, fontsize=7.5)
    save(fig, fig_dir, f"fig_matriz_confusao_{tag}")


def fig_nested(res, oof, fig_dir, tag):
    n = res["nested_cv"]
    folds = n["folds"]
    aucs = np.array([f["AUC"] for f in folds])
    ncols = 2 if oof is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(3.6 * ncols, 3.2), squeeze=False)
    a1 = axes[0, 0]
    x = np.arange(1, len(aucs) + 1)
    a1.axhspan(aucs.mean() - aucs.std(), aucs.mean() + aucs.std(), color=C_SEL, alpha=0.12, lw=0)
    a1.axhline(aucs.mean(), color=C_SEL, lw=2)
    a1.plot(x, aucs, "o", ms=8, color=C_SEL, mec="white", mew=2, zorder=5)
    for xi, f in zip(x, folds):
        a1.text(xi, f["AUC"] + 0.015, f"{f['n_features']} TFs", ha="center", va="bottom",
                fontsize=7, color=INK_2)
    a1.axhline(0.5, ls="--", lw=1, color=C_REF)
    a1.set(xticks=x, xlabel="Outer fold", ylabel="AUC", ylim=(0.45, 1.0), xlim=(0.4, len(x) + 0.6),
           title=f"Nested CV: AUC {aucs.mean():.2f} ± {aucs.std():.2f}")
    a1.grid(axis="x", visible=False)

    if oof is not None:
        a2 = axes[0, 1]
        y, p = oof["y"].values, oof["p"].values
        fpr, tpr, _ = roc_curve(y, p)
        a2.plot(fpr, tpr, color=C_SEL, lw=2, label=f"Pooled out-of-fold\nAUC {roc_auc_score(y, p):.2f}")
        a2.plot([0, 1], [0, 1], ls="--", lw=1, color=C_REF)
        a2.set(xlabel="1 − specificity", ylabel="Sensitivity", title="Nested CV, pooled ROC",
               xlim=(-0.01, 1.01), ylim=(-0.01, 1.01))
        a2.set_aspect("equal")
        a2.legend(loc="lower right")
    fig.text(0.5, -0.06, "Boruta and model fitting repeated inside each outer training fold; "
             "band = mean ± SD.", ha="center", color=INK_2, fontsize=7.5)
    save(fig, fig_dir, f"fig_nested_cv_{tag}")


def stability_series(res, out_dir, balancing):
    p = _path(out_dir, "estabilidade_boruta_frequencia", balancing, "csv")
    if p:
        s = pd.read_csv(p, index_col=0).iloc[:, 0]
    else:  # execuções antigas: só top20 + os 11 originais no JSON
        e = res["estabilidade"]
        s = pd.Series({**e["top20"], **e.get("frequencia_11_TFs_originais", {})})
    return s.sort_values(ascending=False)


def fig_stability(res, out_dir, balancing, symbols, signature, fig_dir, tag):
    if "estabilidade" not in res:
        print("  sem estabilidade no JSON; pulando")
        return
    s = stability_series(res, out_dir, balancing)
    orig = {_norm(g) for g in signature}
    keep = list(s.index[:TOP_N_STABILITY]) + [g for g in s.index if _norm(g) in orig]
    s = s.loc[list(dict.fromkeys(keep))].sort_values()
    n_in = sum(_norm(g) in orig for g in s.index)
    if n_in < len(orig):
        print(f"  [aviso] só {n_in} de {len(orig)} TFs da assinatura encontrados na estabilidade")

    fig, ax = plt.subplots(figsize=(4.6, 0.19 * len(s) + 1.0))
    colors = [C_SEL if _norm(g) in orig else C_MUTED for g in s.index]
    ax.barh(range(len(s)), s.values, height=0.72, color=colors, edgecolor="white", linewidth=1)
    ax.axvline(STABLE_CUTOFF, ls="--", lw=1, color=C_REF)
    ax.text(STABLE_CUTOFF, len(s) - 0.2, f" {STABLE_CUTOFF:.0%}", color=INK_2, fontsize=7.5,
            va="bottom")
    for i, v in enumerate(s.values):
        ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=7, color=INK_2)
    ax.set_yticks(range(len(s)), [label(g, symbols) for g in s.index])
    ax.set(xlim=(0, 1.1), xlabel="Selection frequency",
           title=f"Boruta stability, {res['estabilidade']['n_reamostragens']} resamplings")
    ax.grid(axis="y", visible=False)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=C_SEL),
                       plt.Rectangle((0, 0), 1, 1, color=C_MUTED)],
              labels=[f"In the signature ({len(orig)} TFs)", "Not in the signature"],
              loc="upper center", bbox_to_anchor=(0.4, -0.12 * 30 / (len(s) + 8)), ncol=1)
    save(fig, fig_dir, f"fig_estabilidade_{tag}")


def fig_importance(res, out_dir, balancing, symbols, signature, fig_dir, tag):
    p = _path(out_dir, "importancia_features_modelo_final", balancing, "csv")
    if not p:
        print("  sem importancia_features_modelo_final; pulando")
        return
    s = pd.read_csv(p, index_col=0).iloc[:, 0].sort_values()
    orig = {_norm(g) for g in signature}
    fig, ax = plt.subplots(figsize=(4.2, 0.22 * len(s) + 1.0))
    ax.barh(range(len(s)), s.values, height=0.72, edgecolor="white", linewidth=1,
            color=[C_SEL if _norm(g) in orig else C_MUTED for g in s.index])
    ax.set_yticks(range(len(s)), [label(g, symbols) for g in s.index])
    ax.set(xlabel="Mean decrease in impurity", title="Final model, feature importance")
    ax.grid(axis="y", visible=False)
    if orig and not all(_norm(g) in orig for g in s.index):
        ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=C_SEL),
                           plt.Rectangle((0, 0), 1, 1, color=C_MUTED)],
                  labels=["In the signature", "Other"], loc="upper center",
                  bbox_to_anchor=(0.4, -0.25 * 10 / (len(s) + 4)), ncol=2)
    save(fig, fig_dir, f"fig_importancia_{tag}")


def fig_balancing(out_dir, symbols, signature, fig_dir):
    res = {}
    for b in ("class_weight", "undersample"):
        p = _path(out_dir, "resultados", b, "json")
        if p and b in os.path.basename(p):
            res[b] = json.load(open(p))
    if len(res) < 2:
        print("  precisa de resultados_class_weight.json e resultados_undersample.json; pulando")
        return
    names = {"class_weight": ("Class weighting", C_SEL), "undersample": ("Undersampling", C_ALL)}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.8, 3.8), gridspec_kw={"width_ratios": [1, 1.6]})

    # painel A: AUC por fold do nested CV
    for k, (b, r) in enumerate(res.items()):
        aucs = np.array([f["AUC"] for f in r["nested_cv"]["folds"]])
        name, c = names[b]
        xs = np.full(len(aucs), k) + np.linspace(-0.12, 0.12, len(aucs))
        a1.plot(xs, aucs, "o", ms=8, color=c, mec="white", mew=2, zorder=5)
        a1.plot([k - 0.25, k + 0.25], [aucs.mean()] * 2, color=c, lw=2)
        a1.text(k, 0.47, f"{aucs.mean():.2f} ± {aucs.std():.2f}", ha="center", fontsize=7.5,
                color=INK_2)
    a1.axhline(0.5, ls="--", lw=1, color=C_REF)
    a1.set_xticks(range(len(res)), [names[b][0] for b in res])
    a1.set(ylabel="Nested CV AUC (per fold)", ylim=(0.45, 1.0), xlim=(-0.6, len(res) - 0.4),
           title="Discrimination")
    a1.grid(axis="x", visible=False)

    # painel B: estabilidade dos TFs da assinatura nos dois modos
    freq = {b: stability_series(r, out_dir, b) for b, r in res.items() if "estabilidade" in r}
    orig = {_norm(g) for g in signature}
    common = set.intersection(*[set(f.index) for f in freq.values()]) if freq else set()
    genes = [g for g in common if _norm(g) in orig]
    if len(genes) < len(orig):
        print(f"  [aviso] só {len(genes)} de {len(orig)} TFs da assinatura nas duas estabilidades")
    genes = sorted(genes, key=lambda g: freq.get("class_weight", list(freq.values())[0])[g])
    y = np.arange(len(genes))
    for k, (b, f) in enumerate(freq.items()):
        name, c = names[b]
        a2.barh(y + (0.2 if k == 0 else -0.2), [f[g] for g in genes], height=0.38, color=c,
                edgecolor="white", linewidth=1, label=name)
    a2.axvline(STABLE_CUTOFF, ls="--", lw=1, color=C_REF)
    a2.set_yticks(y, [label(g, symbols) for g in genes])
    a2.set(xlim=(0, 1.05), xlabel="Selection frequency", title="Stability of the signature TFs")
    a2.grid(axis="y", visible=False)
    for k, f in enumerate(freq.values()):
        for yi, g in zip(y, genes):
            a2.text(f[g] + 0.01, yi + (0.2 if k == 0 else -0.2), f"{f[g]:.2f}", va="center",
                    fontsize=6.5, color=INK_2)
    a2.set_xlim(0, 1.12)
    a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    save(fig, fig_dir, "fig_balanceamento_class_weight_vs_undersample")


# ---------------------------------------------------------------- Figura 3 (painel do artigo)
def fig3_panel(res, pred, oof, out_dir, balancing, symbols, signature, fig_dir):
    """Figura 3 do manuscrito: (A) ROC e PR no teste; (B) matrizes de confusão;
    (C) nested CV (AUC por fold e ROC das predições out-of-fold); (D) estabilidade.
    Salva TIFF (LZW), PNG e PDF a 300 dpi."""
    h = res["holdout"]
    y = pred["y"].values
    n_sig = len(h["boruta_confirmed"])
    models = [("p_selecionado", "modelo_selecionado", C_SEL, f"Signature ({n_sig} TFs)"),
              ("p_356", "modelo_356_TFs", C_ALL, "All 356 MRs")]

    fig = plt.figure(figsize=(11.5, 8.2))
    outer = fig.add_gridspec(1, 2, width_ratios=[2.35, 1])
    left = outer[0].subgridspec(2, 1, height_ratios=[1, 1])
    top = left[0].subgridspec(1, 3, width_ratios=[1, 1, 1.05])
    bot = left[1].subgridspec(1, 2, width_ratios=[1.25, 1])
    aR, aP = fig.add_subplot(top[0]), fig.add_subplot(top[1])
    cm_gs = top[2].subgridspec(2, 1)
    aC1, aC2 = fig.add_subplot(cm_gs[0]), fig.add_subplot(cm_gs[1])
    aN, aO = fig.add_subplot(bot[0]), fig.add_subplot(bot[1])
    aS = fig.add_subplot(outer[1])

    # (A) ROC e PR
    prev = y.mean()
    for col, key, color, name in models:
        m = h[key]
        fpr, tpr, _ = roc_curve(y, pred[col])
        aR.plot(fpr, tpr, color=color, lw=2, label=f"{name}, {ci_txt(m)}")
        aR.plot(1 - m["specificity"], m["sensitivity"], "o", ms=7, color=color, mec="white", mew=1.5, zorder=5)
        pr, rc, _ = precision_recall_curve(y, pred[col])
        aP.plot(rc, pr, color=color, lw=2, label=f"{name}, PR-AUC {m['PR_AUC']:.2f}")
    aR.plot([0, 1], [0, 1], ls="--", lw=1, color=C_REF)
    d = h["diferenca_AUC_selecionado_vs_356"]
    aR.text(0.97, 0.04, f"ΔAUC {d['diff']:+.2f} [{d['ci95'][0]:+.2f}; {d['ci95'][1]:+.2f}]\n"
            f"paired bootstrap {p_txt(d['p_bootstrap'])}", ha="right", va="bottom", color=INK_2, fontsize=7)
    aR.set(xlabel="1 − specificity", ylabel="Sensitivity", xlim=(-0.01, 1.01), ylim=(-0.01, 1.01), title="ROC, test set")
    aR.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=7)
    aP.axhline(prev, ls="--", lw=1, color=C_REF)
    aP.text(0.98, prev + 0.02, f"prevalence = {prev:.2f}", ha="right", va="bottom", color=INK_2, fontsize=7)
    aP.set(xlabel="Recall", ylabel="Precision", xlim=(-0.01, 1.01), ylim=(0, 1.01), title="Precision–recall, test set")
    aP.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=7)
    for a in (aR, aP):
        a.set_aspect("equal")

    # (B) matrizes de confusão
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("seq", SEQ)
    for ax, (key, title) in zip((aC1, aC2), [("modelo_selecionado", f"Signature ({n_sig} TFs)"),
                                             ("modelo_356_TFs", "All 356 MRs")]):
        m = h[key]; cm = m["confusion_matrix"]
        mat = np.array([[cm["TN"], cm["FP"]], [cm["FN"], cm["TP"]]])
        rows = mat / mat.sum(axis=1, keepdims=True)
        ax.imshow(rows, cmap=cmap, vmin=0, vmax=1)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{mat[i, j]} ({rows[i, j]:.0%})", ha="center", va="center",
                        color="white" if rows[i, j] > 0.55 else INK, fontsize=7.5)
        ax.set_xticks([0, 1], ["M0", "M1"]); ax.set_yticks([0, 1], ["M0", "M1"])
        ax.set_ylabel("Observed", fontsize=7.5)
        ax.set_title(f"{title}\nsens {m['sensitivity']:.2f} · spec {m['specificity']:.2f}", fontsize=7.5)
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_visible(False)
    aC2.set_xlabel("Predicted", fontsize=7.5)

    # (C) nested CV
    folds = res["nested_cv"]["folds"]
    aucs = np.array([f["AUC"] for f in folds]); x = np.arange(1, len(aucs) + 1)
    aN.axhspan(aucs.mean() - aucs.std(), aucs.mean() + aucs.std(), color=C_SEL, alpha=0.12, lw=0)
    aN.axhline(aucs.mean(), color=C_SEL, lw=2)
    aN.plot(x, aucs, "o", ms=7, color=C_SEL, mec="white", mew=1.5, zorder=5)
    for xi, f in zip(x, folds):
        aN.text(xi, f["AUC"] + 0.015, f"{f['n_features']} TFs", ha="center", va="bottom", fontsize=7, color=INK_2)
    aN.axhline(0.5, ls="--", lw=1, color=C_REF)
    aN.set(xticks=x, xlabel="Outer fold", ylabel="AUC", ylim=(0.45, 1.0), xlim=(0.4, len(x) + 0.6),
           title=f"Nested CV: AUC {aucs.mean():.2f} ± {aucs.std():.2f}")
    aN.grid(axis="x", visible=False)
    if oof is not None:
        yo, po = oof["y"].values, oof["p"].values
        fpr, tpr, _ = roc_curve(yo, po)
        aO.plot(fpr, tpr, color=C_SEL, lw=2, label=f"Pooled out-of-fold, AUC {roc_auc_score(yo, po):.2f}")
        aO.plot([0, 1], [0, 1], ls="--", lw=1, color=C_REF)
        aO.legend(loc="lower right", fontsize=7)
    aO.set(xlabel="1 − specificity", ylabel="Sensitivity", title="Nested CV, pooled ROC",
           xlim=(-0.01, 1.01), ylim=(-0.01, 1.01))
    aO.set_aspect("equal")

    # (D) estabilidade
    if "estabilidade" in res:
        st = stability_series(res, out_dir, balancing)
        orig = {_norm(g) for g in signature}
        keep = list(st.index[:TOP_N_STABILITY]) + [g for g in st.index if _norm(g) in orig]
        st = st.loc[list(dict.fromkeys(keep))].sort_values()
        aS.barh(range(len(st)), st.values, height=0.72, edgecolor="white", linewidth=1,
                color=[C_SEL if _norm(g) in orig else C_MUTED for g in st.index])
        aS.axvline(STABLE_CUTOFF, ls="--", lw=1, color=C_REF)
        for i, v in enumerate(st.values):
            aS.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=6.5, color=INK_2)
        aS.set_yticks(range(len(st)), [label(g, symbols) for g in st.index], fontsize=7.5)
        aS.set(xlim=(0, 1.12), xlabel="Selection frequency",
               title=f"Boruta stability, {res['estabilidade']['n_reamostragens']} resamplings")
        aS.grid(axis="y", visible=False)
        aS.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=C_SEL), plt.Rectangle((0, 0), 1, 1, color=C_MUTED)],
                  labels=["Signature TF", "Not selected"], loc="lower right", fontsize=7)

    for ax, lab, (dx, dy) in [(aR, "A", (-0.28, 1.1)), (aC1, "B", (-0.55, 1.32)),
                              (aN, "C", (-0.12, 1.06)), (aS, "D", (-0.42, 1.03))]:
        ax.text(dx, dy, lab, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="left")

    os.makedirs(fig_dir, exist_ok=True)
    name = os.path.join(fig_dir, "Figure3")
    fig.savefig(name + ".png", dpi=300)
    fig.savefig(name + ".pdf")
    fig.savefig(name + ".tiff", dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("  salvo: Figure3.tiff / .png / .pdf (300 dpi)")


# ---------------------------------------------------------------- execução
def make_all_plots(out_dir=OUT_DIR, balancing="class_weight"):
    fig_dir = os.path.join(out_dir, "figuras")
    p = _path(out_dir, "resultados", balancing, "json")
    if not p or balancing not in os.path.basename(p):
        print(f"resultados_{balancing}.json não encontrado em {out_dir}")
        return
    res = json.load(open(p))
    symbols, signature = load_symbols(), load_signature()
    tag = balancing
    print(f"Figuras ({balancing}) -> {fig_dir}")

    pp = _path(out_dir, "predicoes_teste", balancing, "csv")
    if pp:
        if not os.path.basename(pp).endswith(f"_{balancing}.csv"):
            print("  [aviso] predicoes_teste.csv sem sufixo: confira se é deste modo de balanceamento")
        fig_roc_pr(res, pd.read_csv(pp, index_col=0), fig_dir, tag)
    else:
        print("  sem predicoes_teste; pulando ROC/PR")
    fig_confusion(res, fig_dir, tag)
    po = _path(out_dir, "predicoes_nested_oof", balancing, "csv")
    fig_nested(res, pd.read_csv(po, index_col=0) if po else None, fig_dir, tag)
    fig_stability(res, out_dir, balancing, symbols, signature, fig_dir, tag)
    fig_importance(res, out_dir, balancing, symbols, signature, fig_dir, tag)
    if pp:
        try:
            fig3_panel(res, pd.read_csv(pp, index_col=0), pd.read_csv(po, index_col=0) if po else None,
                       out_dir, balancing, symbols, signature, fig_dir)
        except Exception as e:
            print(f"  [aviso] Figura 3 não gerada: {e}")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    if arg == "compare":
        fig_balancing(OUT_DIR, load_symbols(), load_signature(), os.path.join(OUT_DIR, "figuras"))
    else:
        try:
            import ml_pipeline_revision as P
            default = P.BALANCING
        except Exception:
            default = "class_weight"
        make_all_plots(OUT_DIR, arg or default)