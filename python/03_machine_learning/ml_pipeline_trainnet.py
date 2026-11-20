"""ANÁLISE DE SENSIBILIDADE - Revisor 2, Comentário 2 (Table S26).
Cópia do ml_pipeline_revision.py com 3 mudanças: EXPR_FILE (candidatos da rede inferida só no treino,
R/02_network_MRA/03b_rtn_network_train_only.R), OUT_DIR e N_STABILITY = 0.
Mesmo split treino/teste do artigo (determinístico). Use APENAS o resultado do holdout:
o nested CV não é independente (a rede só-treino contém os pacientes dos folds) e a validação no ICGC
não é válida (o arquivo do ICGC só tem os 356 MRs originais).
"""
"""
Pipeline revisado - classificação M1 vs M0 (TCGA-KIRC)
Responde a: R1.m2, R2.M1, R2.adic (métricas, undersampling), R3.1, R3.2, R3.11, R3.12, R3.13

Diferenças em relação à versão submetida:
  1. Split estratificado (M x grau) feito ANTES de qualquer balanceamento
  2. Balanceamento apenas no treino (class_weight ou undersampling)
  3. Boruta apenas no treino; modelo final avaliado UMA vez no teste intocado
  4. Limiar de decisão (Youden) escolhido no treino, aplicado uma vez no teste
  5. Modelos de 11 TFs, 356 TFs e baseline avaliados no MESMO conjunto de teste
  6. Nested CV + estabilidade de seleção
  7. Seeds fixas e tudo salvo em OUT_DIR (arquivos com sufixo do modo de balanceamento)
  8. Figuras geradas ao final (plots_ml_revision.py)
"""
import json
import os

import numpy as np
import pandas as pd

# BorutaPy antigo usa aliases removidos do numpy >= 1.24
np.int = int
np.float = float
np.bool = bool

from boruta import BorutaPy
from imblearn.under_sampling import RandomUnderSampler
from scipy.stats import chi2_contingency
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, average_precision_score,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score, roc_curve)
from sklearn.model_selection import StratifiedKFold, train_test_split

# ---------------------------------------------------------------- config
SEED = 2026
EXPR_FILE = "0_data/trainnet/expr_MR_KIRC_trainnet.csv"   # amostras x 365 MRs da rede SÓ-TREINO (03b)
META_FILE = "0_data/metadata_KIRC_full.csv"     # índice = mesmo barcode do EXPR_FILE
ORIGINAL_FEATURES = ""   # assinatura original (11 TFs) fora do artigo; "" pula as etapas dela
OUT_DIR = "results_ml_trainnet"

M_COL = "ajcc_pathologic_m"
STRAT_COL = "neoplasm.histologic.grade"       # ajuste ao nome real; None = estratificar só por M
STRAT_LOW = ["G1", "G2", "GX"]
STRAT_HIGH = ["G3", "G4"]

BALANCING = "class_weight"      # "class_weight" (principal) ou "undersample" (comparação)
TEST_SIZE = 0.30
RF_PARAMS = dict(n_estimators=1000, max_depth=5, n_jobs=-1)
BORUTA_MAX_ITER = 100
N_BOOT = 2000                   # bootstrap para IC 95%
N_OUTER = 5                     # folds externos do nested CV
MAKE_PLOTS = True               # gera as figuras em OUT_DIR/figuras (plots_ml_revision.py)
N_STABILITY = 0                 # sem estabilidade nesta análise de sensibilidade

os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------- dados
def load_data():
    X = pd.read_csv(EXPR_FILE, index_col=0)
    meta = pd.read_csv(META_FILE, index_col=0)
    flow = {"meta_inicial": len(meta)}

    meta = meta[meta[M_COL].isin(["M0", "M1"])]
    flow["apos_filtro_M0_M1"] = len(meta)

    meta = meta[meta.index.str[13:15] == "01"]           # tumor primário
    flow["apos_filtro_tumor_primario"] = len(meta)

    meta = meta[~meta.index.str[:12].duplicated(keep="first")]  # 1 amostra por paciente
    flow["apos_1_amostra_por_paciente"] = len(meta)

    common = meta.index.intersection(X.index)
    X, meta = X.loc[common], meta.loc[common]
    flow["com_expressao"] = len(meta)

    y = (meta[M_COL] == "M1").astype(int)
    flow["M0"], flow["M1"] = int((y == 0).sum()), int((y == 1).sum())
    return X, y, meta, flow


def make_strata(meta, y, min_size=10):
    if STRAT_COL is None or STRAT_COL not in meta.columns:
        print(f"[aviso] coluna de estratificação '{STRAT_COL}' ausente; estratificando só por M")
        return y.astype(str)
    g = np.where(meta[STRAT_COL].isin(STRAT_LOW), "low",
                 np.where(meta[STRAT_COL].isin(STRAT_HIGH), "high", "unk"))
    strata = pd.Series(y.astype(str).values + "_" + g, index=y.index)

    # estratos pequenos (ex.: grau GX ou ausente) vão para o maior estrato do mesmo M,
    # para que todo estrato sobreviva a split 70/30 + CV 5-fold + reamostragem 80%
    counts = strata.value_counts()
    for s in counts[counts < min_size].index:
        m = s.split("_")[0]
        same_m = counts[[c for c in counts.index if c.startswith(m + "_") and counts[c] >= min_size]]
        strata[strata == s] = same_m.idxmax() if len(same_m) else m
    n_unk = int((g == "unk").sum())
    if n_unk:
        print(f"[info] {n_unk} paciente(s) sem grau G1-G4 alocados ao maior estrato do mesmo M")
    return strata


# ---------------------------------------------------------------- modelos
def make_rf(seed):
    params = dict(RF_PARAMS, random_state=seed)
    if BALANCING == "class_weight":
        params["class_weight"] = "balanced_subsample"
    return RandomForestClassifier(**params)


def balance(X, y, seed):
    if BALANCING == "undersample":
        rus = RandomUnderSampler(random_state=seed)
        Xb, yb = rus.fit_resample(X, y)
        return Xb, yb
    return X, y


def boruta_select(Xtr, ytr, seed):
    Xb, yb = balance(Xtr, ytr, seed)
    sel = BorutaPy(make_rf(seed), n_estimators="auto", max_iter=BORUTA_MAX_ITER,
                   random_state=seed, verbose=0)
    sel.fit(np.asarray(Xb), np.asarray(yb))
    confirmed = Xtr.columns[sel.support_].tolist()
    tentative = Xtr.columns[sel.support_weak_].tolist()
    return confirmed, tentative


def fit_predict(Xtr, ytr, Xte, feats, seed):
    Xb, yb = balance(Xtr[feats], ytr, seed)
    rf = make_rf(seed).fit(Xb, yb)
    return rf.predict_proba(Xte[feats])[:, 1], rf


def train_threshold(Xtr, ytr, feats, seed, n_splits=5):
    """Limiar de Youden escolhido SÓ no treino (predições out-of-fold)."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.zeros(len(ytr))
    for tr, va in skf.split(Xtr, ytr):
        p, _ = fit_predict(Xtr.iloc[tr], ytr.iloc[tr], Xtr.iloc[va], feats, seed)
        oof[va] = p
    fpr, tpr, thr = roc_curve(ytr, oof)
    return float(thr[np.argmax(tpr - fpr)])


# ---------------------------------------------------------------- métricas
def boot_ci(y, p, fn, n=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    y, p = np.asarray(y), np.asarray(p)
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) < 2:
            continue
        vals.append(fn(y[i], p[i]))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return [float(lo), float(hi)]


def paired_auc_diff(y, p1, p2, n=N_BOOT, seed=SEED):
    """Bootstrap pareado da diferença de AUC (p1 - p2)."""
    rng = np.random.default_rng(seed)
    y, p1, p2 = map(np.asarray, (y, p1, p2))
    diffs = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) < 2:
            continue
        diffs.append(roc_auc_score(y[i], p1[i]) - roc_auc_score(y[i], p2[i]))
    diffs = np.array(diffs)
    p_two_sided = 2 * min((diffs <= 0).mean(), (diffs >= 0).mean())
    return {"diff": float(roc_auc_score(y, p1) - roc_auc_score(y, p2)),
            "ci95": [float(x) for x in np.percentile(diffs, [2.5, 97.5])],
            "p_bootstrap": float(min(p_two_sided, 1.0))}


def all_metrics(y, p, thr=0.5):
    y = np.asarray(y)
    pred = (np.asarray(p) >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    out = {
        "n": int(len(y)), "n_M1": int(y.sum()),
        "AUC": float(roc_auc_score(y, p)),
        "AUC_ci95": boot_ci(y, p, roc_auc_score),
        "PR_AUC": float(average_precision_score(y, p)),
        "accuracy": float(accuracy_score(y, pred)),
        "accuracy_ci95": boot_ci(y, p, lambda a, b: accuracy_score(a, b >= thr)),
        "sensitivity": float(recall_score(y, pred, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else float("nan"),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "F1": float(f1_score(y, pred, zero_division=0)),
        "confusion_matrix": {"TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp)},
    }
    return out


# ---------------------------------------------------------------- etapas
def split_report(y, strata, idx_tr, idx_te):
    tab = pd.crosstab(pd.Series(np.where(y.index.isin(idx_tr), "treino", "teste"), index=y.index),
                      strata)
    chi2, p, _, _ = chi2_contingency(tab)
    tab.to_csv(f"{OUT_DIR}/split_distribuicao_estratos_{BALANCING}.csv")
    return {"tabela": tab.to_dict(), "chi2": float(chi2), "p": float(p)}


def holdout(X, y, strata):
    # split por posições: o sklearn não indexa índices de string do pyarrow.
    # Mesmo n, estratos e seed => exatamente o mesmo split das execuções anteriores
    pos_tr, pos_te = train_test_split(np.arange(len(y)), test_size=TEST_SIZE,
                                      stratify=np.asarray(strata), random_state=SEED)
    idx_tr, idx_te = y.index[pos_tr], y.index[pos_te]
    Xtr, Xte, ytr, yte = X.iloc[pos_tr], X.iloc[pos_te], y.iloc[pos_tr], y.iloc[pos_te]

    confirmed, tentative = boruta_select(Xtr, ytr, SEED)
    p_sel, rf_sel = fit_predict(Xtr, ytr, Xte, confirmed, SEED)
    p_all, _ = fit_predict(Xtr, ytr, Xte, X.columns.tolist(), SEED)
    p_major = np.full(len(yte), ytr.mean())  # baseline: prevalência do treino
    thr_sel = train_threshold(Xtr, ytr, confirmed, SEED)
    thr_all = train_threshold(Xtr, ytr, X.columns.tolist(), SEED)

    res = {
        "split": split_report(y, strata, idx_tr, idx_te),
        "n_treino": int(len(ytr)), "n_teste": int(len(yte)),
        "M1_treino": int(ytr.sum()), "M1_teste": int(yte.sum()),
        "boruta_confirmed": confirmed, "boruta_tentative": tentative,
        "modelo_selecionado": {**all_metrics(yte, p_sel, thr_sel), "limiar_treino": thr_sel},
        "modelo_356_TFs": {**all_metrics(yte, p_all, thr_all), "limiar_treino": thr_all},
        "baseline_majoritario": {"AUC": 0.5, "accuracy": float(accuracy_score(yte, (p_major >= 0.5).astype(int)))},
        "diferenca_AUC_selecionado_vs_356": paired_auc_diff(yte, p_sel, p_all),
    }
    pd.DataFrame({"y": yte, "p_selecionado": p_sel, "p_356": p_all}).to_csv(
        f"{OUT_DIR}/predicoes_teste_{BALANCING}.csv")
    pd.Series(rf_sel.feature_importances_, index=confirmed).sort_values(ascending=False).to_csv(
        f"{OUT_DIR}/importancia_features_modelo_final_{BALANCING}.csv")
    return res, Xtr, ytr, strata.iloc[pos_tr]


def nested_cv(X, y, strata):
    outer = StratifiedKFold(n_splits=N_OUTER, shuffle=True, random_state=SEED)
    oof = pd.Series(np.nan, index=y.index)
    folds, counts = [], pd.Series(0, index=X.columns)
    for k, (tr, te) in enumerate(outer.split(np.zeros(len(y)), np.asarray(strata))):
        Xtr, Xte, ytr, yte = X.iloc[tr], X.iloc[te], y.iloc[tr], y.iloc[te]
        feats, _ = boruta_select(Xtr, ytr, SEED + k)
        if not feats:
            print(f"[aviso] fold {k}: Boruta não confirmou features; usando todas")
            feats = X.columns.tolist()
        counts[feats] += 1
        p, _ = fit_predict(Xtr, ytr, Xte, feats, SEED + k)
        oof.iloc[te] = p
        thr = train_threshold(Xtr, ytr, feats, SEED + k)
        m = all_metrics(yte, p, thr)
        folds.append({"fold": k, "n_teste": int(len(yte)), "M1_teste": int(yte.sum()),
                      "n_features": len(feats), "AUC": m["AUC"],
                      "sensitivity": m["sensitivity"], "specificity": m["specificity"],
                      "limiar_treino": thr})
    aucs = [f["AUC"] for f in folds]
    counts.sort_values(ascending=False).to_csv(f"{OUT_DIR}/nested_cv_frequencia_selecao_{BALANCING}.csv")
    pd.DataFrame({"y": y, "p": oof}).to_csv(f"{OUT_DIR}/predicoes_nested_oof_{BALANCING}.csv")
    return {"folds": folds, "AUC_media": float(np.mean(aucs)), "AUC_dp": float(np.std(aucs)),
            "AUC_pooled_oof": float(roc_auc_score(y, oof.values))}


def stability(Xtr, ytr, strata_tr, original):
    counts = pd.Series(0, index=Xtr.columns)
    for r in range(N_STABILITY):
        pos, _ = train_test_split(np.arange(len(ytr)), train_size=0.8,
                                  stratify=np.asarray(strata_tr), random_state=SEED + 100 + r)
        feats, _ = boruta_select(Xtr.iloc[pos], ytr.iloc[pos], SEED + 100 + r)
        counts[feats] += 1
    freq = (counts / N_STABILITY).sort_values(ascending=False)
    freq.to_csv(f"{OUT_DIR}/estabilidade_boruta_frequencia_{BALANCING}.csv")
    out = {"n_reamostragens": N_STABILITY, "top20": freq.head(20).round(3).to_dict()}
    if original:
        out["frequencia_11_TFs_originais"] = freq.reindex(original).round(3).to_dict()
    return out


def fixed_signature(X, y, strata, feats, label="original_11"):
    """Avalia uma assinatura FIXA (sem Boruta) no mesmo split e nos mesmos folds do pipeline.

    Atenção: os 11 TFs originais foram escolhidos na versão submetida, com dados que incluem
    pacientes que agora estão no teste. A estimativa é, portanto, otimista e serve para
    reportar o desempenho da assinatura original com o protocolo correto, não como validação.
    """
    feats = [f for f in feats if f in X.columns]
    pos_tr, pos_te = train_test_split(np.arange(len(y)), test_size=TEST_SIZE,
                                      stratify=np.asarray(strata), random_state=SEED)
    Xtr, Xte, ytr, yte = X.iloc[pos_tr], X.iloc[pos_te], y.iloc[pos_tr], y.iloc[pos_te]

    p_fix, _ = fit_predict(Xtr, ytr, Xte, feats, SEED)
    thr = train_threshold(Xtr, ytr, feats, SEED)
    p_all, _ = fit_predict(Xtr, ytr, Xte, X.columns.tolist(), SEED)
    out = {"features": feats, "n_features": len(feats),
           "nota": "assinatura definida na versão submetida com pacientes que agora estão no teste: "
                   "estimativa otimista",
           "holdout": {**all_metrics(yte, p_fix, thr), "limiar_treino": thr},
           "diferenca_AUC_vs_356": paired_auc_diff(yte, p_fix, p_all)}

    # junta às predições do holdout para as figuras e compara com a assinatura nova, se existir
    pfile = f"{OUT_DIR}/predicoes_teste_{BALANCING}.csv"
    if os.path.exists(pfile):
        pred = pd.read_csv(pfile, index_col=0)
        if len(pred) == len(yte) and list(pred.index.astype(str)) == list(yte.index.astype(str)):
            pred[f"p_{label}"] = p_fix
            pred.to_csv(pfile)
            out["diferenca_AUC_vs_selecionado"] = paired_auc_diff(yte, p_fix, pred["p_selecionado"].values)
        else:
            print(f"[aviso] {pfile} não corresponde ao split atual; comparação com a assinatura nova pulada")

    # nested CV com features fixas: mesmos folds externos, sem seleção
    outer = StratifiedKFold(n_splits=N_OUTER, shuffle=True, random_state=SEED)
    aucs = []
    for k, (tr, te) in enumerate(outer.split(np.zeros(len(y)), np.asarray(strata))):
        p, _ = fit_predict(X.iloc[tr], y.iloc[tr], X.iloc[te], feats, SEED + k)
        aucs.append(float(roc_auc_score(y.iloc[te], p)))
    out["nested_cv"] = {"AUC_folds": aucs, "AUC_media": float(np.mean(aucs)), "AUC_dp": float(np.std(aucs))}
    return out


# ---------------------------------------------------------------- validação externa (ICGC)
EXT_EXPR = "0_data/expr_MR_ICGC_RT91.csv"    # 91 tumores RT x genes (mesmos IDs do TCGA); R/01_data/06b_make_ICGC_external_RT91.R
EXT_META = "0_data/metadata_ICGC_RT91.csv"   # gerado pelo mesmo script 06b
EXT_M_COL = "pathologic_m"
N_PERM = 5000


def zscore(df):
    return (df - df.mean()) / df.std(ddof=0).replace(0, 1)


def external_validation(X_tcga, y_tcga, feats):
    """Modelo travado: treinado em TODO o TCGA com as features dadas, aplicado ao ICGC sem
    reseleção nem reajuste. z-score por gene dentro de cada coorte (não usa rótulos)."""
    Xe = pd.read_csv(EXT_EXPR, index_col=0)
    me = pd.read_csv(EXT_META, index_col=0)
    if len(set(feats) & set(Xe.index)) > len(set(feats) & set(Xe.columns)):   # veio genes x amostras
        Xe = Xe.T
    if not set(feats) & set(Xe.columns):
        raise ValueError(f"Nenhuma feature encontrada em {EXT_EXPR}: confira se os IDs (ENSG) batem com os do TCGA")
    # somente tumores primários (RT), 1 por doador: amostras RA herdariam o M do doador (revisão R2, C1/C3)
    if "tissue_type" in me.columns:
        me = me[me["tissue_type"] == "RT"]
    if "submitted_donor_id" in me.columns:
        me = me[~me["submitted_donor_id"].duplicated()]
    me = me[me[EXT_M_COL].isin(["M0", "M1"])]
    common = me.index.astype(str).intersection(Xe.index.astype(str))
    Xe.index, me.index = Xe.index.astype(str), me.index.astype(str)
    Xe, ye = Xe.loc[common], (me.loc[common, EXT_M_COL] == "M1").astype(int)

    avail = [f for f in feats if f in Xe.columns]
    missing = sorted(set(feats) - set(avail))
    Xt, Xez = zscore(X_tcga[avail]), zscore(Xe[avail])
    thr = train_threshold(Xt, y_tcga, avail, SEED)
    p, _ = fit_predict(Xt, y_tcga, Xez, avail, SEED)

    auc = roc_auc_score(ye, p)
    rng = np.random.default_rng(SEED)
    perm = np.array([roc_auc_score(rng.permutation(ye.values), p) for _ in range(N_PERM)])
    return {"n": int(len(ye)), "n_M1": int(ye.sum()), "features_usadas": avail,
            "features_ausentes": missing, **all_metrics(ye, p, thr),
            "p_permutacao": float((np.sum(perm >= auc) + 1) / (N_PERM + 1)),
            "predicoes": dict(zip(ye.index, map(float, p)))}


def run_fixed_and_external(X, y, strata, original, results):
    """Etapas rápidas (sem Boruta): assinatura original fixa e validação no ICGC."""
    if original:
        print("[extra] Assinatura original (fixa) no mesmo split...")
        results["assinatura_original_11"] = fixed_signature(X, y, strata, original)
    if os.path.exists(EXT_EXPR) and os.path.exists(EXT_META):
        print("[extra] Validação externa no ICGC...")
        ext = {}
        if original:
            ext["assinatura_original_11"] = external_validation(X, y, original)
        if "holdout" in results:
            ext["assinatura_nova"] = external_validation(X, y, results["holdout"]["boruta_confirmed"])
        for k, v in ext.items():
            pd.Series(v.pop("predicoes"), name="p").to_csv(f"{OUT_DIR}/predicoes_ICGC_{k}_{BALANCING}.csv")
        results["externa_ICGC"] = ext
    else:
        print(f"[extra] {EXT_EXPR} ou {EXT_META} ausente: validação no ICGC pulada")
    return results


def print_extra(results):
    if "assinatura_original_11" in results:
        o = results["assinatura_original_11"]
        m = o["holdout"]
        print(f"Original 11 TFs (fixa, {o['n_features']} presentes): AUC {m['AUC']:.3f} "
              f"[{m['AUC_ci95'][0]:.3f}-{m['AUC_ci95'][1]:.3f}] | nested CV "
              f"{o['nested_cv']['AUC_media']:.3f} ± {o['nested_cv']['AUC_dp']:.3f}")
        print("  vs 356 TFs:", o["diferenca_AUC_vs_356"])
        if "diferenca_AUC_vs_selecionado" in o:
            print("  vs assinatura nova:", o["diferenca_AUC_vs_selecionado"])
    for k, v in results.get("externa_ICGC", {}).items():
        print(f"ICGC, {k}: n = {v['n']} ({v['n_M1']} M1), {len(v['features_usadas'])} features, "
              f"AUC {v['AUC']:.3f} [{v['AUC_ci95'][0]:.3f}-{v['AUC_ci95'][1]:.3f}], "
              f"permutação p = {v['p_permutacao']:.3f}; ausentes: {v['features_ausentes']}")


# ---------------------------------------------------------------- baselines (R3.12)
# O modelo do artigo usa a EXPRESSÃO dos 356 MRs (EXPR_FILE). O baseline "sem priorização pela rede"
# repete o mesmo procedimento (Boruta só no treino + RF) na expressão de TODOS os reguladores da rede.
BASE_EXPR_FILE = "0_data/expr_all_regulators_voom_KIRC.csv"   # amostras x todos os reguladores (voom)
CLINICAL_JSON = "results_clinical_baseline/resultados_clinical_baseline.json"   # opcional, do clinical_baseline.py


def _key(c):
    return str(c).split(".")[0].upper()


def run_baselines():
    """`python ml_pipeline_revision.py baselines`: no MESMO split, compara a assinatura com
    (i) classe majoritária e (ii) o mesmo procedimento sem a priorização da MRA (todos os reguladores)."""
    X, y, meta, flow = load_data()
    strata = make_strata(meta, y)
    jfile = f"{OUT_DIR}/resultados_{BALANCING}.json"
    if not os.path.exists(jfile):
        raise SystemExit(f"{jfile} não encontrado: rode o pipeline principal antes")
    sig = json.load(open(jfile))["holdout"]["boruta_confirmed"]

    E = pd.read_csv(BASE_EXPR_FILE, index_col=0)
    E.columns = [_key(c) for c in E.columns]
    E = E.loc[:, ~E.columns.duplicated()]
    if len(y.index.intersection(E.index)) < len(y):
        raise SystemExit("A matriz de todos os reguladores não cobre todas as amostras do pipeline: confira os IDs")
    E = E.loc[y.index]
    n_mr_in = sum(_key(c) in set(E.columns) for c in X.columns)
    print(f"Todos os reguladores: {E.shape[1]} TFs (contém {n_mr_in} dos {X.shape[1]} MRs); {len(y)} amostras")

    pos_tr, pos_te = train_test_split(np.arange(len(y)), test_size=TEST_SIZE,
                                      stratify=np.asarray(strata), random_state=SEED)
    ytr, yte = y.iloc[pos_tr], y.iloc[pos_te]
    Xtr, Xte, Etr, Ete = X.iloc[pos_tr], X.iloc[pos_te], E.iloc[pos_tr], E.iloc[pos_te]

    out = {"n_teste": int(len(yte)), "M1_teste": int(yte.sum()), "n_reguladores": int(E.shape[1]),
           "n_MRs": int(X.shape[1]), "assinatura": sig}
    p_sig, _ = fit_predict(Xtr, ytr, Xte, sig, SEED)
    thr_sig = train_threshold(Xtr, ytr, sig, SEED)
    out["assinatura_MRA"] = {**all_metrics(yte, p_sig, thr_sig), "limiar_treino": thr_sig, "n_features": len(sig)}
    out["classe_majoritaria"] = {"AUC": 0.5, "PR_AUC": float(yte.mean()),
                                 "accuracy": float(accuracy_score(yte, np.zeros(len(yte), int))),
                                 "sensitivity": 0.0, "specificity": 1.0}
    print("Boruta na expressão de todos os reguladores (treino)...")
    conf_e, _ = boruta_select(Etr, ytr, SEED)
    preds = pd.DataFrame({"y": yte, "p_assinatura": p_sig})
    if conf_e:
        p_all, _ = fit_predict(Etr, ytr, Ete, conf_e, SEED)
        thr_all = train_threshold(Etr, ytr, conf_e, SEED)
        out["sem_priorizacao"] = {**all_metrics(yte, p_all, thr_all), "limiar_treino": thr_all,
                                  "features": conf_e, "n_features": len(conf_e),
                                  "n_features_que_sao_MRs": int(sum(c in {_key(m) for m in X.columns} for c in conf_e))}
        out["diferenca_AUC_assinatura_vs_sem_priorizacao"] = paired_auc_diff(yte, p_sig, p_all)
        preds["p_sem_priorizacao"] = p_all
    else:
        out["sem_priorizacao"] = {"nota": "Boruta não confirmou nenhuma feature"}
    preds.to_csv(f"{OUT_DIR}/predicoes_teste_baselines_{BALANCING}.csv")
    if os.path.exists(CLINICAL_JSON):
        out["clinico"] = json.load(open(CLINICAL_JSON))
    with open(f"{OUT_DIR}/resultados_baselines_{BALANCING}.json", "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

    print("\n=== Baselines (teste intocado, n = %d; %d M1) ===" % (len(yte), yte.sum()))
    rows = [("Classe majoritária", out["classe_majoritaria"]),
            (f"Assinatura MRA ({len(sig)} TFs)", out["assinatura_MRA"])]
    if "AUC_ci95" in out["sem_priorizacao"]:
        sp = out["sem_priorizacao"]
        rows.append((f"Sem priorização ({sp['n_features']} TFs; {sp['n_features_que_sao_MRs']} MRs)", sp))
    for name, m in rows:
        ci = f" [{m['AUC_ci95'][0]:.2f}-{m['AUC_ci95'][1]:.2f}]" if "AUC_ci95" in m else ""
        print(f"{name:40s} AUC {m['AUC']:.2f}{ci} | PR-AUC {m['PR_AUC']:.2f} | acc {m['accuracy']:.2f} "
              f"| sens {m['sensitivity']:.2f} | esp {m['specificity']:.2f}")
    if "diferenca_AUC_assinatura_vs_sem_priorizacao" in out:
        d = out["diferenca_AUC_assinatura_vs_sem_priorizacao"]
        print(f"ΔAUC assinatura − sem priorização: {d['diff']:.3f} [{d['ci95'][0]:.3f}, {d['ci95'][1]:.3f}], "
              f"p = {d['p_bootstrap']:.3f}")
        print("Features sem priorização:", out["sem_priorizacao"]["features"])
    return out


def main_extra():
    """`python ml_pipeline_revision.py extra`: roda só as etapas rápidas e atualiza o JSON
    existente, sem refazer Boruta, nested CV e estabilidade."""
    X, y, meta, flow = load_data()
    strata = make_strata(meta, y)
    original = [l.strip() for l in open(ORIGINAL_FEATURES) if l.strip()] if os.path.exists(ORIGINAL_FEATURES) else []
    jfile = f"{OUT_DIR}/resultados_{BALANCING}.json"
    results = json.load(open(jfile)) if os.path.exists(jfile) else {}
    results = run_fixed_and_external(X, y, strata, original, results)
    with open(jfile, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print_extra(results)
    if MAKE_PLOTS and "holdout" in results:
        try:
            import plots_ml_revision as PL
            PL.make_all_plots(OUT_DIR, BALANCING)
        except Exception as e:
            print(f"[aviso] figuras não geradas: {e}")


def main():
    X, y, meta, flow = load_data()
    print("Fluxo de pacientes:", flow)
    strata = make_strata(meta, y)
    print("Estratos:\n", strata.value_counts())

    original = []
    if os.path.exists(ORIGINAL_FEATURES):
        original = [l.strip() for l in open(ORIGINAL_FEATURES) if l.strip()]

    results = {"config": {"seed": SEED, "balancing": BALANCING, "strat_col": STRAT_COL,
                          "test_size": TEST_SIZE, "rf": RF_PARAMS},
               "fluxo_pacientes": flow}

    print("\n[1/3] Holdout com teste intocado...")
    results["holdout"], Xtr, ytr, strata_tr = holdout(X, y, strata)
    if original:
        sel = set(results["holdout"]["boruta_confirmed"])
        results["holdout"]["sobreposicao_com_11_originais"] = sorted(sel & set(original))

    print("[2/3] Nested CV...")
    results["nested_cv"] = nested_cv(X, y, strata)

    if N_STABILITY > 0:
        print("[3/3] Estabilidade do Boruta...")
        results["estabilidade"] = stability(Xtr, ytr, strata_tr, original)

    results = run_fixed_and_external(X, y, strata, original, results)

    with open(f"{OUT_DIR}/resultados_{BALANCING}.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    if MAKE_PLOTS:
        try:
            import plots_ml_revision as PL
            PL.make_all_plots(OUT_DIR, BALANCING)
        except Exception as e:  # figura nunca derruba o pipeline
            print(f"[aviso] figuras não geradas: {e}")

    h = results["holdout"]
    print("\n=== Resumo ===")
    print(f"Features confirmadas ({len(h['boruta_confirmed'])}):", h["boruta_confirmed"])
    for name in ["modelo_selecionado", "modelo_356_TFs"]:
        m = h[name]
        print(f"{name}: AUC {m['AUC']:.3f} [{m['AUC_ci95'][0]:.3f}-{m['AUC_ci95'][1]:.3f}] "
              f"| sens {m['sensitivity']:.2f} | esp {m['specificity']:.2f} | PR-AUC {m['PR_AUC']:.3f}")
    print("Diferença AUC (selecionado - 356):", h["diferenca_AUC_selecionado_vs_356"])
    n = results["nested_cv"]
    print(f"Nested CV: AUC {n['AUC_media']:.3f} ± {n['AUC_dp']:.3f}")
    print_extra(results)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "extra":
        main_extra()
    elif len(sys.argv) > 1 and sys.argv[1] == "baselines":
        run_baselines()
    else:
        main()