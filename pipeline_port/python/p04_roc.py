"""
p04_roc.py
==========
원본: machine learning modeling 2.R  (pROC 로 유전자별 ROC/AUC + bootstrap 95% CI)

입력 : 학습/검증 발현행렬('{샘플}_{그룹}' 라벨) + 유전자 리스트(LASSO.gene.txt)
출력 : ROC 곡선 PDF (유전자별) + roc_auc_summary.csv

CLI:
  python p04_roc.py --expr data/data.train.txt --genes results/ml/LASSO.gene.txt \
      --title Training --outdir results/roc
  python p04_roc.py --expr data/data.test.txt  --genes results/ml/LASSO.gene.txt \
      --title Validation --outdir results/roc
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import config
from io_utils import read_expr_labeled, labels_from_columns


def bootstrap_auc_ci(y, score, n=config.ROC_BOOTSTRAP, seed=config.RANDOM_SEED):
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    y = np.asarray(y); score = np.asarray(score)
    idx = np.arange(len(y))
    aucs = []
    for _ in range(n):
        b = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[b])) < 2:
            continue
        aucs.append(roc_auc_score(y[b], score[b]))
    lo, hi = np.percentile(aucs, [2.5, 97.5])
    return lo, hi


def run_roc(expr_file, gene_file, title, outdir):
    from sklearn.metrics import roc_auc_score, roc_curve
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    df = read_expr_labeled(expr_file)
    labels = labels_from_columns(df.columns)
    y = (labels == config.CASE_LABEL).astype(int)
    genes = [g.strip() for g in open(gene_file) if g.strip()]

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        have_mpl = True
    except Exception:
        have_mpl = False

    rows = []
    for g in genes:
        if g not in df.index:
            continue
        score = df.loc[g].values.astype(float)
        # AUC 방향 자동 정렬 (pROC 기본 동작과 동일하게 >0.5 유지)
        auc = roc_auc_score(y, score)
        if auc < 0.5:
            score = -score
            auc = roc_auc_score(y, score)
        lo, hi = bootstrap_auc_ci(y, score)
        rows.append({"gene": g, "AUC": round(auc, 3),
                     "CI_low": round(lo, 3), "CI_high": round(hi, 3)})
        if have_mpl:
            fpr, tpr, _ = roc_curve(y, score)
            plt.figure(figsize=(3.5, 3.5))
            plt.plot(fpr, tpr, color="red")
            plt.plot([0, 1], [0, 1], ls="--", c="grey", lw=0.6)
            plt.title(f"{g} of {title}")
            plt.text(0.45, 0.15, f"AUC={auc:.3f}\n95%CI {lo:.3f}-{hi:.3f}",
                     color="red", fontsize=8)
            plt.xlabel("1 - Specificity"); plt.ylabel("Sensitivity")
            plt.tight_layout()
            plt.savefig(outdir / f"{title}_ROC.{g}.pdf")
            plt.close()

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / f"roc_auc_summary_{title}.csv", index=False)
    print(f"[ROC] {title}:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--title", default="Training")
    ap.add_argument("--outdir", default=str(config.RESULT_DIR / "roc"))
    a = ap.parse_args()
    run_roc(a.expr, a.genes, a.title, a.outdir)
