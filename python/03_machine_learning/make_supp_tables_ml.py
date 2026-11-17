"""
Tabelas suplementares do machine learning, a partir dos resultados já salvos (não refaz nada).

  Tabela S (split)      distribuição dos estratos no treino e no teste + n e M1 de cada fold
                        do nested CV  (R1.m2, R3.1)
  Tabela S (baselines)  desempenho da assinatura e dos modelos de comparação no teste
                        (R3.12, R2.M3, R3.8)

Uso:  python make_supp_tables_ml.py
Saídas (CSV, prontas para colar no Word/Excel) em results_ml_revision/tabelas_suplementares/
"""
import json
import os

import pandas as pd

OUT_DIR = "results_ml_revision"
BALANCING = "class_weight"
MAIN_JSON = f"{OUT_DIR}/resultados_{BALANCING}.json"
BASE_JSON = f"{OUT_DIR}/resultados_baselines_{BALANCING}.json"
CLIN_JSON = f"{OUT_DIR}/clinico/resultados_clinico_{BALANCING}.json"
SPLIT_CSV = f"{OUT_DIR}/split_distribuicao_estratos_{BALANCING}.csv"
TAB_DIR = f"{OUT_DIR}/tabelas_suplementares"
os.makedirs(TAB_DIR, exist_ok=True)


def load(path):
    if os.path.exists(path):
        return json.load(open(path))
    print(f"[aviso] {path} não encontrado; a parte correspondente fica de fora")
    return {}


def f2(x):
    return "" if x is None else f"{x:.2f}"


def auc_ci(m):
    if not m or "AUC" not in m:
        return ""
    ci = m.get("AUC_ci95")
    return f"{m['AUC']:.2f} ({ci[0]:.2f}–{ci[1]:.2f})" if ci else f"{m['AUC']:.2f}"


def acc_ci(m):
    ci = m.get("accuracy_ci95")
    return f"{m['accuracy']:.2f} ({ci[0]:.2f}–{ci[1]:.2f})" if ci else f2(m.get("accuracy"))


def p_txt(p):
    return "< 0.001" if p < 0.001 else f"{p:.3f}" if p < 0.01 else f"{p:.2f}"


def diff_txt(d, label):
    if not d:
        return ""
    return f"{label}: {d['diff']:+.3f} ({d['ci95'][0]:+.3f} to {d['ci95'][1]:+.3f}); p = {p_txt(d['p_bootstrap'])}"


main, base, clin = load(MAIN_JSON), load(BASE_JSON), load(CLIN_JSON)
h = main.get("holdout", {})

# ======================================================================== Tabela S: split e folds
rows = []
if h:
    rows += [{"Set": "Training", "n": h["n_treino"], "M1": h["M1_treino"], "M0": h["n_treino"] - h["M1_treino"]},
             {"Set": "Test", "n": h["n_teste"], "M1": h["M1_teste"], "M0": h["n_teste"] - h["M1_teste"]}]
split = pd.DataFrame(rows)
if os.path.exists(SPLIT_CSV):                        # contagem por estrato (M x grau)
    st = pd.read_csv(SPLIT_CSV, index_col=0)
    st.index = [{"treino": "Training", "teste": "Test"}.get(i, i) for i in st.index]
    split = split.merge(st, left_on="Set", right_index=True, how="left")
if "split" in h:
    split["χ² p (strata, training vs test)"] = [p_txt(h["split"]["p"]), ""]
n_total = (h.get("n_treino", 0) + h.get("n_teste", 0)) or None
m1_total = (h.get("M1_treino", 0) + h.get("M1_teste", 0)) or None

folds = pd.DataFrame(main.get("nested_cv", {}).get("folds", []))
if len(folds):
    folds = pd.DataFrame({
        "Outer fold": folds["fold"] + 1,
        "Training n": n_total - folds["n_teste"] if n_total else "",
        "Training M1": m1_total - folds["M1_teste"] if m1_total else "",
        "Test n": folds["n_teste"], "Test M1": folds["M1_teste"],
        "Selected TFs": folds["n_features"],
        "AUC": folds["AUC"].round(2), "Sensitivity": folds["sensitivity"].round(2),
        "Specificity": folds["specificity"].round(2)})

with open(f"{TAB_DIR}/TableS_split_e_folds.csv", "w", encoding="utf-8") as fh:
    fh.write("Holdout split (stratified by metastatic status and grade)\n")
    split.to_csv(fh, index=False)
    fh.write("\nNested cross-validation (5 stratified outer folds)\n")
    folds.to_csv(fh, index=False)
print("=== Split ===\n", split.to_string(index=False))
print("\n=== Folds do nested CV ===\n", folds.to_string(index=False))

# ======================================================================== Tabela S: baselines
def row(model, features, m, n_test, comparison="", nested=""):
    if not m:
        return None
    return {"Model": model, "Features": features, "Test set": n_test,
            "AUC (95% CI)": auc_ci(m), "PR-AUC": f2(m.get("PR_AUC")),
            "Accuracy (95% CI)": acc_ci(m) if "accuracy" in m else "",
            "Sensitivity": f2(m.get("sensitivity")), "Specificity": f2(m.get("specificity")),
            "Precision": f2(m.get("precision")), "F1": f2(m.get("F1")),
            "Nested CV AUC (mean ± SD)": nested, "Comparison (ΔAUC, paired bootstrap)": comparison}


n_test = f"n = {h.get('n_teste', '')} ({h.get('M1_teste', '')} M1)"
nested_sig = main.get("nested_cv", {})
nested_sig = f"{nested_sig['AUC_media']:.2f} ± {nested_sig['AUC_dp']:.2f}" if "AUC_media" in nested_sig else ""
n_sig = len(h.get("boruta_confirmed", []))
sp = base.get("sem_priorizacao", {})

tab = [
    row("Majority class", "None (always M0)", base.get("classe_majoritaria"), n_test),
    row("Signature", f"{n_sig} TFs selected by Boruta from the 356 MRs", h.get("modelo_selecionado"), n_test,
        nested=nested_sig),
    row("All master regulators", "356 MRs (no selection)", h.get("modelo_356_TFs"), n_test,
        diff_txt(h.get("diferenca_AUC_selecionado_vs_356"), "signature − 356 MRs")),
    row("Without network prioritisation",
        f"{sp.get('n_features', '')} TFs selected by Boruta from all {base.get('n_reguladores', '')} regulators"
        + (f" ({sp['n_features_que_sao_MRs']} are MRs)" if "n_features_que_sao_MRs" in sp else ""),
        sp if "AUC" in sp else None, n_test,
        diff_txt(base.get("diferenca_AUC_assinatura_vs_sem_priorizacao"), "signature − without prioritisation")),
]

ch = clin.get("holdout", {})
if ch:
    cm, cc, cn = ch.get("metricas", {}), ch.get("comparacoes", {}), clin.get("nested_cv", {}).get("resumo", {})
    n_c = f"n = {ch.get('n_teste', '')} ({ch.get('M1_teste', '')} M1); complete clinical data"
    nst = lambda k: f"{cn[k]['AUC_media']:.2f} ± {cn[k]['AUC_dp']:.2f}" if k in cn else ""
    tab += [
        row("Clinical", "Grade, pT, age (logistic regression)", cm.get("clinico"), n_c, nested=nst("clinico")),
        row("Signature (complete-case subset)", "TFs selected by Boruta", cm.get("TFs"), n_c,
            diff_txt(cc.get("TFs_vs_clinico"), "signature − clinical"), nst("TFs")),
        row("Clinical + signature", "Grade, pT, age + signature score", cm.get("clinico_TFs"), n_c,
            diff_txt(cc.get("clinico_TFs_vs_clinico"), "combined − clinical"), nst("clinico_TFs")),
    ]
tab = pd.DataFrame([r for r in tab if r])
tab.to_csv(f"{TAB_DIR}/TableS_baselines.csv", index=False, encoding="utf-8")

# regressão multivariável (escore da assinatura ajustado por grau, pT e idade)
mv = clin.get("multivariavel")
if mv:
    t = pd.DataFrame(mv["tabela"]).T
    t = t.rename(index={"grade": "Grade", "pT": "Pathological T stage", "age": "Age",
                        "tf_score": "Signature score", "const": "Intercept"})
    t = t.drop(index="Intercept", errors="ignore")
    out = pd.DataFrame({"Variable": t.index,
                        "OR per SD (95% CI)": [f"{r['OR_por_DP']:.2f} ({r['IC95_inf']:.2f}–{r['IC95_sup']:.2f})"
                                               for _, r in t.iterrows()],
                        "p": [p_txt(float(r["p"])) for _, r in t.iterrows()]})
    lrt = mv.get("LRT_tf_score", {})
    with open(f"{TAB_DIR}/TableS_baselines.csv", "a", encoding="utf-8") as fh:
        fh.write("\nMultivariable logistic regression for M1 vs M0 (complete-case cohort)\n")
        out.to_csv(fh, index=False)
        if lrt:
            fh.write(f"Likelihood-ratio test for the signature score: chi2 = {lrt['chi2']:.2f}, "
                     f"df = {lrt['gl']}, p = {p_txt(lrt['p'])}\n")

print("\n=== Baselines ===\n", tab.to_string(index=False))
print(f"\nTabelas salvas em {TAB_DIR}/")