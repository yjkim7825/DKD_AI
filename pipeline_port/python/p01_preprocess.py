"""
p01_preprocess.py
=================
원본: data preprocessing 1/2/3.R  (probe매핑 -> 정규화 -> 배치보정 병합)

CLI 예시:
  # 1) GEO series matrix + 플랫폼 주석 -> gene x sample 발현행렬
  python p01_preprocess.py probe2gene \
      --series data/GSE37263_series_matrix.txt.gz \
      --platform data/GPL5175.txt --symbol-col "Gene Symbol" --symbol-index 1 \
      --out data/GSE37263_geneMatrix.txt

  # 2) 정규화 + 3군 라벨링 (Control/Early/Late 샘플 리스트 지정)
  python p01_preprocess.py normalize \
      --matrix data/GSE142025_geneMatrix.txt \
      --control data/s1.txt --early data/s2.txt --late data/s3.txt \
      --geoid GSE142025 --out data/GSE142025_threeGroups.normalize.txt

  # 3) 여러 데이터셋 교집합 유전자 병합 + ComBat 배치보정
  python p01_preprocess.py combat \
      --inputs data/GSE96804_labeled.txt data/GSE104948_labeled.txt \
      --out data/merge.normalize.txt
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd

from deg_utils import quantile_normalize, auto_log2, combat
from io_utils import (load_series_matrix, load_platform_map, avereps,
                      read_expr_labeled, write_expr_labeled)


def cmd_probe2gene(a):
    expr = load_series_matrix(a.series)                    # probe x sample
    mapping = load_platform_map(a.platform, id_col=a.id_col,
                                symbol_col=a.symbol_col,
                                symbol_index=a.symbol_index)
    common = expr.index.intersection(mapping.index)
    expr = expr.loc[common]
    genes = mapping.loc[common].values
    gene_expr = avereps(expr, genes)                       # gene x sample
    gene_expr = gene_expr.dropna(how="all")
    write_expr_labeled(gene_expr, a.out)
    print(f"[probe2gene] {gene_expr.shape[0]} genes x {gene_expr.shape[1]} samples -> {a.out}")


def _read_sample_list(path):
    return [s.strip() for s in open(path) if s.strip()]


def cmd_normalize(a):
    df = pd.read_csv(a.matrix, sep="\t", index_col=0)
    df = df.groupby(level=0).mean()                        # avereps
    mat = auto_log2(df.values)
    mat = quantile_normalize(mat)                          # normalizeBetweenArrays
    norm = pd.DataFrame(mat, index=df.index, columns=df.columns)

    groups = []
    cols = []
    for path, label in [(a.control, "Control"), (a.early, "Early"), (a.late, "Late")]:
        if path is None:
            continue
        for s in _read_sample_list(path):
            if s in norm.columns:
                cols.append(s)
                groups.append(label)
    out = norm[cols].copy()
    out.columns = [f"{c}_{g}" for c, g in zip(cols, groups)]
    write_expr_labeled(out, a.out)
    print(f"[normalize] {out.shape} -> {a.out}  (labels: {pd.Series(groups).value_counts().to_dict()})")


def cmd_combat(a):
    frames = [read_expr_labeled(p) for p in a.inputs]
    # 교집합 유전자
    inter = frames[0].index
    for f in frames[1:]:
        inter = inter.intersection(f.index)
    inter = sorted(inter)

    batch = []
    parts = []
    for i, f in enumerate(frames):
        sub = f.loc[inter]
        sub.columns = [f"B{i+1}_{c}" for c in sub.columns]
        parts.append(sub)
        batch += [i] * sub.shape[1]
    merged = pd.concat(parts, axis=1)

    corrected = combat(merged.values, np.array(batch))
    out = pd.DataFrame(corrected, index=merged.index, columns=merged.columns)
    write_expr_labeled(out, a.out)
    print(f"[combat] {out.shape}  ({len(frames)} batches) -> {a.out}")


def build_parser():
    p = argparse.ArgumentParser(description="P01 preprocessing (probe->gene / normalize / combat)")
    sub = p.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("probe2gene")
    p1.add_argument("--series", required=True)
    p1.add_argument("--platform", required=True)
    p1.add_argument("--id-col", default="ID")
    p1.add_argument("--symbol-col", default="Gene Symbol")
    p1.add_argument("--symbol-index", type=int, default=0,
                    help="'A // B' 형태 심볼에서 사용할 토큰 인덱스 (R의 [,2]는 1)")
    p1.add_argument("--out", required=True)
    p1.set_defaults(func=cmd_probe2gene)

    p2 = sub.add_parser("normalize")
    p2.add_argument("--matrix", required=True)
    p2.add_argument("--control")
    p2.add_argument("--early")
    p2.add_argument("--late")
    p2.add_argument("--geoid", default="dataset")
    p2.add_argument("--out", required=True)
    p2.set_defaults(func=cmd_normalize)

    p3 = sub.add_parser("combat")
    p3.add_argument("--inputs", nargs="+", required=True)
    p3.add_argument("--out", required=True)
    p3.set_defaults(func=cmd_combat)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)
