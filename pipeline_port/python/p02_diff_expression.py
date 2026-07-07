"""
p02_diff_expression.py
======================
원본: differential expression analysis.R  (limma lmFit + eBayes + topTable, 2-group)

입력 : '{샘플}_{그룹}' 라벨 발현행렬 (예: merge.normalize.txt / *_twoGroups.normalize.txt)
출력 : all_*.txt (전체), diff_*.txt (유의 DEG), up/down 유전자 리스트, 화산도(선택)

CLI:
  python p02_diff_expression.py \
      --input data/GSE142025_twoGroups.normalize.txt \
      --ref Control --alt DKD --tag DKD_vs_Control
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

import config
from deg_utils import limma_two_group
from io_utils import read_expr_labeled, labels_from_columns


def run_deg(input_file, ref, alt, tag, outdir,
            logfc=config.LOGFC_FILTER, adjp=config.ADJP_FILTER, plot=True):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_expr_labeled(input_file)
    labels = labels_from_columns(df.columns)
    keep = np.isin(labels, [ref, alt])
    df = df.loc[:, keep]
    labels = labels[keep]
    print(f"[DEG] groups: {pd.Series(labels).value_counts().to_dict()}")

    res = limma_two_group(df.values, labels, ref=ref, alt=alt)
    table = pd.DataFrame(res, index=df.index).sort_values("P.Value")
    table.to_csv(outdir / f"all_{tag}.txt", sep="\t")

    sig = table[(table["logFC"].abs() > logfc) & (table["adj.P.Val"] < adjp)]
    sig.to_csv(outdir / f"diff_{tag}.txt", sep="\t")
    up = sig[sig["logFC"] > 0].index.tolist()
    down = sig[sig["logFC"] < 0].index.tolist()
    pd.Series(up).to_csv(outdir / f"up_genes_{tag}.txt", sep="\t", index=False, header=False)
    pd.Series(down).to_csv(outdir / f"down_genes_{tag}.txt", sep="\t", index=False, header=False)
    print(f"[DEG] {len(sig)} DEG (up={len(up)}, down={len(down)}) -> diff_{tag}.txt")

    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            t = table.copy()
            sig_mask = (t["adj.P.Val"] < adjp) & (t["logFC"].abs() > logfc)
            color = np.where(sig_mask & (t["logFC"] > 0), "red",
                     np.where(sig_mask & (t["logFC"] < 0), "blue", "grey"))
            plt.figure(figsize=(5.5, 4.5))
            plt.scatter(t["logFC"], -np.log10(t["adj.P.Val"] + 1e-300),
                        c=color, s=6)
            plt.axvline(logfc, ls="--", c="k", lw=0.5)
            plt.axvline(-logfc, ls="--", c="k", lw=0.5)
            plt.xlabel("log2FC"); plt.ylabel("-log10(adj.P.Val)")
            plt.title(f"{alt} vs {ref}")
            plt.tight_layout()
            plt.savefig(outdir / f"vol_{tag}.pdf")
            plt.close()
        except Exception as e:  # matplotlib 미설치 등
            print(f"[DEG] 화산도 생략: {e}")
    return sig


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--ref", default=config.CONTROL_LABEL)
    ap.add_argument("--alt", default=config.CASE_LABEL)
    ap.add_argument("--tag", default="DKD_vs_Control")
    ap.add_argument("--outdir", default=str(config.RESULT_DIR / "deg"))
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    run_deg(a.input, a.ref, a.alt, a.tag, a.outdir, plot=not a.no_plot)
