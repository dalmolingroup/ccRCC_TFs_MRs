"""Exporta os IDs de treino e teste do split principal (70/30) do ml_pipeline_revision.py.

Revisão CSBJ-D-26-01039 - Revisor 2, Comentário 2.
O split é determinístico (mesmas amostras, mesmos estratos, SEED = 2026), então este script
reproduz exatamente o treino/teste usado no artigo. Rode na mesma pasta do ml_pipeline_revision.py:

    python export_split_ids.py
"""
import json
import os

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

import ml_pipeline_revision as P   # o main() não roda no import

X, y, meta, flow = P.load_data()
strata = P.make_strata(meta, y)
pos_tr, pos_te = train_test_split(np.arange(len(y)), test_size=P.TEST_SIZE,
                                  stratify=np.asarray(strata), random_state=P.SEED)

train = pd.DataFrame({"sample": y.index[pos_tr], "M": y.iloc[pos_tr].map({0: "M0", 1: "M1"}).values})
test = pd.DataFrame({"sample": y.index[pos_te], "M": y.iloc[pos_te].map({0: "M0", 1: "M1"}).values})
train.to_csv("0_data/split_train_ids.csv", index=False)
test.to_csv("0_data/split_test_ids.csv", index=False)

print(f"Treino: {len(train)} ({(train.M == 'M1').sum()} M1) | Teste: {len(test)} ({(test.M == 'M1').sum()} M1)")

# checagem contra o JSON do artigo
jfile = f"{P.OUT_DIR}/resultados_{P.BALANCING}.json"
if os.path.exists(jfile):
    h = json.load(open(jfile)).get("holdout", {})
    ok = h.get("n_treino") == len(train) and h.get("n_teste") == len(test)
    print(f"JSON do artigo: treino = {h.get('n_treino')}, teste = {h.get('n_teste')} -> {'OK' if ok else 'DIFERENTE: confira'}")
print("Salvos: 0_data/split_train_ids.csv, 0_data/split_test_ids.csv")
