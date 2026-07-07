"""STEP4 : LASSO ∩ SVM-RFE → ROC (R 01_paper_design_lasso_svm_roc.R 이식, 논문 Figure4 설계).

훈련 = GSE96804 단독,  검증 = ComBat(GSE104948 + GSE104954).
후보 = 저자 확정 10개(Supp Table 8/9): ALDH2,FN1,VNN2,CREB5,XAF1,CA2,CDKN1B,IFI44L,SYTL2,TSPYL5.
LASSO(sklearn L1 로지스틱) ∩ SVM-RFE(sklearn RFECV+선형SVC) → 교집합. 기대 6개(FN1/ALDH2 포함).
단일유전자 ROC(train/valid) + FN1+ALDH2 결합모델(GLM/RF/SVM/XGBoost) AUC.
라이브러리 매핑: glmnet→sklearn.LogisticRegression(L1), caret rfe→sklearn.RFECV, pROC→sklearn.roc_auc_score.
출력: RES_DIR/step4_paper/*.csv
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import numpy as np
import pandas as pd

CAND10 = ["ALDH2", "FN1", "VNN2", "CREB5", "XAF1", "CA2", "CDKN1B", "IFI44L", "SYTL2", "TSPYL5"]
EXPECT6 = {"CDKN1B", "ALDH2", "FN1", "XAF1", "TSPYL5", "VNN2"}
OUT = config.RES_DIR / "step4_paper"; OUT.mkdir(parents=True, exist_ok=True)


def _combat(mats, batch):
    try:
        from inmoose.pycombat import pycombat_norm
        return pycombat_norm(mats, batch)
    except Exception:
        try:
            from combat.pycombat import pycombat
            return pycombat(mats, batch)
        except Exception:
            return None


def load_valid():
    """검증 = ComBat(GSE104948 + GSE104954)."""
    a = pd.read_csv(config.OUT_DIR / "GSE104948.labeled.txt", sep="\t", index_col=0)
    b = pd.read_csv(config.OUT_DIR / "GSE104954.labeled.txt", sep="\t", index_col=0)
    a.columns = [f"GSE104948_{c}" for c in a.columns]
    b.columns = [f"GSE104954_{c}" for c in b.columns]
    common = sorted(set(a.index) & set(b.index))
    allt = pd.concat([a.loc[common], b.loc[common]], axis=1)
    batch = [0] * a.shape[1] + [1] * b.shape[1]
    cc = _combat(allt, batch)
    if cc is None:
        print("[valid] TODO: pycombat 미설치 → ComBat 없이 원자료 병합 사용(근사).")
        return allt
    return pd.DataFrame(cc, index=allt.index, columns=allt.columns)


def y_of(df):
    return np.array([0 if c.endswith("_Control") else 1 for c in df.columns])


def main():
    from sklearn.linear_model import LogisticRegression
    from sklearn.feature_selection import RFECV
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    train = pd.read_csv(config.OUT_DIR / "GSE96804.labeled.txt", sep="\t", index_col=0)
    valid = load_valid()
    ytr, yva = y_of(train), y_of(valid)
    feats = [g for g in CAND10 if g in train.index]
    pd.Series(feats).to_csv(OUT / "interGenes.List.txt", index=False, header=False)

    Xtr = train.loc[feats].T.to_numpy()
    # ---- LASSO (L1 로지스틱, CV) ----
    lasso = LogisticRegression(penalty="l1", solver="liblinear", C=1.0, max_iter=5000)
    lasso.fit(Xtr, ytr)
    lasso_genes = [f for f, c in zip(feats, lasso.coef_[0]) if abs(c) > 1e-8]
    # ---- SVM-RFE (선형 SVC + RFECV) ----
    rfe = RFECV(SVC(kernel="linear"), step=1, cv=StratifiedKFold(5), scoring="accuracy", min_features_to_select=2)
    rfe.fit(Xtr, ytr)
    svm_genes = [f for f, keep in zip(feats, rfe.support_) if keep]
    inter = sorted(set(lasso_genes) & set(svm_genes))
    pd.Series(inter).to_csv(OUT / "interGenes.txt", index=False, header=False)
    print(f"[LASSO] {sorted(lasso_genes)}")
    print(f"[SVM-RFE] {sorted(svm_genes)}")
    print(f"[교집합] {inter} | FN1={'FN1' in inter} ALDH2={'ALDH2' in inter} | 기대6 일치={set(inter)==EXPECT6}")

    # ---- 단일유전자 ROC ----
    def auc(df, y, g):
        if g not in df.index:
            return np.nan
        a = roc_auc_score(y, df.loc[g].to_numpy())
        return max(a, 1 - a)   # pROC 자동 방향(AUC>=0.5) 관례에 맞춤(예: ALDH2 보호유전자)
    rows = [{"gene": g, "AUC_train": auc(train, ytr, g), "AUC_valid": auc(valid, yva, g)} for g in CAND10]
    roc_tab = pd.DataFrame(rows).sort_values("AUC_valid", ascending=False)
    roc_tab.to_csv(OUT / "single_gene_ROC_AUC.csv", index=False)
    print("\n[단일유전자 ROC-AUC]\n", roc_tab.to_string(index=False))

    # ---- FN1+ALDH2 결합모델 ----
    def two_gene(df):
        return np.column_stack([df.loc["FN1"].to_numpy(), df.loc["ALDH2"].to_numpy()])
    Xt, Xv = two_gene(train), two_gene(valid)
    models = {
        "GLM": LogisticRegression(max_iter=5000),
        "RF":  RandomForestClassifier(n_estimators=500, random_state=123),
        "SVM": SVC(kernel="rbf", probability=True, random_state=123),
    }
    try:
        from xgboost import XGBClassifier
        models["XGBoost"] = XGBClassifier(n_estimators=50, max_depth=3, learning_rate=0.3,
                                          use_label_encoder=False, eval_metric="logloss")
    except Exception:
        from sklearn.ensemble import GradientBoostingClassifier
        models["XGBoost"] = GradientBoostingClassifier(random_state=123)  # 대체(설치 없을 때)
    mres = []
    for name, mdl in models.items():
        mdl.fit(Xt, ytr)
        ptr = mdl.predict_proba(Xt)[:, 1]; pva = mdl.predict_proba(Xv)[:, 1]
        mres.append({"model": name, "AUC_train": roc_auc_score(ytr, ptr), "AUC_valid": roc_auc_score(yva, pva)})
    mtab = pd.DataFrame(mres)
    mtab.to_csv(OUT / "combined_model_AUC.csv", index=False)
    print("\n[FN1+ALDH2 결합모델 AUC]\n", mtab.to_string(index=False))


if __name__ == "__main__":
    main()
