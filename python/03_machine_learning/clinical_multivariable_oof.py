"""Comentário 5 (R2): associação ajustada do escore da assinatura com metástase, usando predições
OUT-OF-FOLD para TODOS os pacientes.

Revisão CSBJ-D-26-01039 - Revisor 2, Comentário 5:
  "clarify whether the regulator-derived score used in the clinical logistic regression was obtained
   through out-of-fold predictions for every relevant patient, and whether its adjusted association
   was evaluated independently of model fitting."

Escore: logito da probabilidade predita no nested CV (predicoes_nested_oof_<BALANCING>.csv).
Cada paciente recebe a predição de um modelo (Boruta + RF) que NÃO o viu no treino, então o escore
de todos os pacientes é out-of-fold e a regressão logística é ajustada depois, sobre predições fixas.
Modelo: M1 ~ grau + pT + idade + escore (por DP), casos completos; LRT para o escore.

Rodar na pasta do ml_pipeline_revision.py:  python clinical_multivariable_oof.py
"""
import json
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2

import ml_pipeline_revision as P

AGE_COL = "age.at.index"
PT_COL = "ajcc_pathologic_T"
GRADE_COL = P.STRAT_COL
OUT = os.path.join(P.OUT_DIR, "clinico")
os.makedirs(OUT, exist_ok=True)


def find_col(meta, name):
    if name in meta.columns:
        return name
    norm = lambda s: s.lower().replace("_", ".")
    hits = [c for c in meta.columns if norm(c) == norm(name)]
    if not hits:
        raise KeyError(f"coluna '{name}' não encontrada")
    return hits[0]


X, y, meta, flow = P.load_data()
gc, tc, ac = find_col(meta, GRADE_COL), find_col(meta, PT_COL), find_col(meta, AGE_COL)
C = pd.DataFrame({
    "grade": meta[gc].astype(str).str.extract(r"G([1-4])", expand=False).astype(float),
    "pT": meta[tc].astype(str).str.extract(r"T([1-4])", expand=False).astype(float),
    "age": pd.to_numeric(meta[ac], errors="coerce")}, index=meta.index)
if C["age"].median() < 0:
    C["age"] = C["age"].abs() / 365.25

# predições out-of-fold do nested CV (todas as 501 amostras)
oof = pd.read_csv(f"{P.OUT_DIR}/predicoes_nested_oof_{P.BALANCING}.csv", index_col=0)
assert oof["p"].notna().all() and len(oof) == len(y), "predições OOF incompletas"
p = oof["p"].clip(1e-4, 1 - 1e-4)
C["tf_score"] = np.log(p / (1 - p)).reindex(C.index)

d = pd.concat([C, y.rename("m1")], axis=1).dropna()
print(f"Casos completos: {len(d)} ({int(d.m1.sum())} M1) de {len(y)}")

Z = d[["grade", "pT", "age", "tf_score"]].apply(lambda s: (s - s.mean()) / s.std())   # OR por DP
full = sm.Logit(d.m1, sm.add_constant(Z)).fit(disp=0)
base = sm.Logit(d.m1, sm.add_constant(Z[["grade", "pT", "age"]])).fit(disp=0)
lr = 2 * (full.llf - base.llf)

ci = full.conf_int()
tab = pd.DataFrame({"OR_per_SD": np.exp(full.params), "CI_low": np.exp(ci[0]), "CI_high": np.exp(ci[1]),
                    "p": full.pvalues}).drop(index="const")
tab.to_csv(f"{OUT}/regressao_multivariavel_OOF.csv")
print("\nRegressão logística multivariável (escore = logito da predição out-of-fold do nested CV):")
print(tab.round(4).to_string())
print(f"\nLRT (escore adicionado ao modelo clínico): chi2 = {lr:.2f}, gl = 1, p = {chi2.sf(lr, 1):.2e}")

json.dump({"n": len(d), "M1": int(d.m1.sum()), "tabela": tab.to_dict(orient="index"),
           "LRT_tf_score": {"chi2": float(lr), "gl": 1, "p": float(chi2.sf(lr, 1))}},
          open(f"{OUT}/regressao_multivariavel_OOF.json", "w"), indent=2)
