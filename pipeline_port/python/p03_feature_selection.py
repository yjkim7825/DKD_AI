"""
p03_feature_selection.py
========================
원본: machine learning modeling 1.R  (LASSO + SVM-RFE + Venn 교집합)

입력 : '{샘플}_{그룹}' 라벨 학습 발현행렬(data.train.txt) + 후보 유전자 리스트(interGenes)
출력 : LASSO.gene.txt, SVM-RFE.gene.txt, interGenes.txt (교집합)

주의(충실도):
  - LASSO  : glmnet(binomial, alpha=1, cv.glmnet lambda.min)
             -> sklearn LogisticRegressionCV(penalty='l1', solver='saga',
                                              cv=10, scoring='neg_log_loss')
  - SVM-RFE: R caret::rfe(svmRadial) -> sklearn RFECV(SVC(kernel='linear')) 로
             선형 SVM-RFE(Guyon et al.) 구현. (radial 커널은 특징 랭킹이 불가하여
             표준 SVM-RFE 관례에 따라 선형 커널 사용 — 결과 유전자 집합은 대개 유사)

CLI:
  python p03_feature_selection.py \
      --train data/data.train.txt --genes data/interGenes.List.txt \
      --outdir results/ml
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import config
from io_utils import read_expr_labeled, labels_from_columns


def _load_xy(train_file, gene_file):
    df = read_expr_labeled(train_file)
    if gene_file:
        genes = [g.strip() for g in open(gene_file) if g.strip()]
        df = df.loc[df.index.intersection(genes)]
    X = df.T.values                      # samples x genes
    genes = df.index.tolist()
    labels = labels_from_columns(df.columns)
    y = (labels == config.CASE_LABEL).astype(int)   # Control=0, DKD=1
    return X, y, genes


def run_lasso(X, y, genes, outdir):
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.preprocessing import StandardScaler
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegressionCV(
        Cs=50, cv=config.CV_FOLDS, penalty="l1", solver="saga",
        scoring="neg_log_loss", max_iter=5000, random_state=config.RANDOM_SEED,
    ).fit(Xs, y)
    coef = clf.coef_.ravel()
    sel = [g for g, c in zip(genes, coef) if abs(c) > 1e-8]
    pd.Series(sel).to_csv(Path(outdir) / "LASSO.gene.txt",
                          sep="\t", index=False, header=False)
    print(f"[LASSO] {len(sel)} genes -> {sel}")
    return sel


def run_svm_rfe(X, y, genes, outdir):
    from sklearn.svm import SVC
    from sklearn.feature_selection import RFECV
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import StratifiedKFold
    Xs = StandardScaler().fit_transform(X)
    min_size = min(config.SVM_RFE_SIZES)
    selector = RFECV(
        estimator=SVC(kernel="linear"),
        step=1, cv=StratifiedKFold(config.CV_FOLDS),
        scoring="accuracy", min_features_to_select=min_size,
    ).fit(Xs, y)
    sel = [g for g, m in zip(genes, selector.support_) if m]
    pd.Series(sel).to_csv(Path(outdir) / "SVM-RFE.gene.txt",
                          sep="\t", index=False, header=False)
    print(f"[SVM-RFE] {len(sel)} genes -> {sel}")
    return sel


def main(train_file, gene_file, outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    np.random.seed(config.RANDOM_SEED)
    X, y, genes = _load_xy(train_file, gene_file)
    print(f"[data] {X.shape[0]} samples x {X.shape[1]} genes "
          f"(Control={int((y==0).sum())}, DKD={int((y==1).sum())})")

    lasso = set(run_lasso(X, y, genes, outdir))
    svm = set(run_svm_rfe(X, y, genes, outdir))
    inter = sorted(lasso & svm)
    pd.Series(inter).to_csv(outdir / "interGenes.txt",
                            sep="\t", index=False, header=False)
    print(f"[Venn] LASSO∩SVM-RFE = {len(inter)} genes -> {inter}")
    return inter


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", required=True)
    ap.add_argument("--genes", default=None, help="후보 유전자 리스트(interGenes.List.txt)")
    ap.add_argument("--outdir", default=str(config.RESULT_DIR / "ml"))
    a = ap.parse_args()
    main(a.train, a.genes, a.outdir)
