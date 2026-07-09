"""STEP1(추가) : GSE96804 complementary DEG (R 03_gse96804_deg.R 이식).

논문 Methods(p.4) "complementary limma(v3.64.1) DEG of GSE96804, adjP<0.05, |logFC|>0.585".
MR 통합(DKD 고발현∩risk / 저발현∩protective)의 근거 → FN1/ALDH2 방향 정의.
⚠️ limma eBayes(moderated t) → Welch t + BH 근사(R 이식 공통 방침). 방향·임계는 R 과 동일.
출력: RES_DIR/step1_deg/ (DEG_all/diff_GSE96804*, DKD_high/low_genes, FN1_ALDH2_direction).
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import numpy as np
import pandas as pd
from scipy import stats

# 02_deg 의 정규화/BH 재사용
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import importlib.util
_spec = importlib.util.spec_from_file_location("deg02", str(pathlib.Path(__file__).resolve().parent / "02_deg.py"))
deg02 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(deg02)

S1 = config.RES_DIR / "step1_deg"; S1.mkdir(parents=True, exist_ok=True)


def main():
    mat = pd.read_csv(config.OUT_DIR / "GSE96804.labeled.txt", sep="\t", index_col=0)
    groups = [c.split("_")[-1] for c in mat.columns]        # Control / DKD
    mat = deg02.quantile_normalize(mat)                     # normalizeBetweenArrays 근사

    cols_d = [c for c, g in zip(mat.columns, groups) if g == config.GROUP_DKD]
    cols_c = [c for c, g in zip(mat.columns, groups) if g == config.GROUP_CONTROL]
    D = mat[cols_d].to_numpy(); C = mat[cols_c].to_numpy()
    logFC = D.mean(1) - C.mean(1)                           # DKD - Control
    t, p = stats.ttest_ind(D, C, axis=1, equal_var=False)
    adj = deg02.bh_adjust(np.nan_to_num(p, nan=1.0))
    res = pd.DataFrame({"id": mat.index, "logFC": logFC, "P.Value": p, "adj.P.Val": adj}).sort_values("adj.P.Val")
    res.to_csv(S1 / "DEG_all_GSE96804_DKD_vs_Control.txt", sep="\t", index=False)

    sig = res[(res["logFC"].abs() > config.LOGFC_FILTER) & (res["adj.P.Val"] < config.ADJP_FILTER)]
    sig.to_csv(S1 / "DEG_diff_GSE96804_DKD_vs_Control.txt", sep="\t", index=False)
    up = sig[sig["logFC"] > 0]["id"].tolist(); dn = sig[sig["logFC"] < 0]["id"].tolist()
    pd.Series(up).to_csv(S1 / "GSE96804.DKD_high_genes.txt", index=False, header=False)
    pd.Series(dn).to_csv(S1 / "GSE96804.DKD_low_genes.txt", index=False, header=False)

    # FN1/ALDH2 방향 (논문: FN1=DKD고발현·risk, ALDH2=저발현·protective)
    foc = res[res["id"].isin(["FN1", "ALDH2"])].copy()
    foc["direction"] = np.where((foc["logFC"].abs() > config.LOGFC_FILTER) & (foc["adj.P.Val"] < config.ADJP_FILTER),
                                np.where(foc["logFC"] > 0, "Up", "Down"), "Not")
    foc.to_csv(S1 / "GSE96804.FN1_ALDH2_direction.csv", index=False)
    print(f"[GSE96804 DKD vs Control] 유의 {len(sig)} (up {len(up)} / down {len(dn)})")
    print(foc.to_string(index=False))


if __name__ == "__main__":
    main()
