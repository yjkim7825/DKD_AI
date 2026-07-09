"""STEP1-2 : limma DEG 3종 (R 02_deg.R 이식).

대비: Late_vs_Early / Late_vs_Control / Early_vs_Control.
필터: |logFC| > LOGFC_FILTER(0.585) & adj.P.Val < ADJP_FILTER(0.05).
⚠️ limma 의 eBayes(moderated t)는 파이썬 직접 대응이 없어 per-gene Welch t-검정 + BH 로 근사.
   방향(alt-ref)·임계값은 R 과 동일. 결과 개수는 eBayes 차이로 소폭 달라질 수 있음(PORT_REPORT 참조).
출력: RES_DIR/DEG_all_*.txt, DEG_diff_*.txt
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import numpy as np
import pandas as pd
from scipy import stats

S1 = config.RES_DIR / "step1_deg"; S1.mkdir(parents=True, exist_ok=True)  # STEP 구분 저장(R 과 동일)


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """limma::normalizeBetweenArrays(quantile) 근사 — 열(샘플) 간 분위수 정규화."""
    ranked_mean = df.stack().groupby(df.rank(method="first").stack().astype(int)).mean()
    return df.rank(method="min").stack().astype(int).map(ranked_mean).unstack()


def bh_adjust(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (R p.adjust(method="BH") 와 동일)."""
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    ranked = p[order] * n / (np.arange(1, n + 1))          # 오름차순 랭크별 p*n/rank
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]      # 위(큰 p)에서부터 단조화
    adj = np.empty(n); adj[order] = np.clip(ranked, 0, 1)   # 원래 순서로 복원
    return adj


def run_contrast(mat, groups, alt, ref, tag):
    cols_a = [c for c, g in zip(mat.columns, groups) if g == alt]
    cols_r = [c for c, g in zip(mat.columns, groups) if g == ref]
    A = mat[cols_a].to_numpy(); R = mat[cols_r].to_numpy()
    logFC = A.mean(1) - R.mean(1)
    t, p = stats.ttest_ind(A, R, axis=1, equal_var=False)
    adj = bh_adjust(np.nan_to_num(p, nan=1.0))
    res = pd.DataFrame({"id": mat.index, "logFC": logFC, "P.Value": p, "adj.P.Val": adj})
    res = res.sort_values("adj.P.Val")
    res.to_csv(S1 / f"DEG_all_{tag}.txt", sep="\t", index=False)
    sig = res[(res["logFC"].abs() > config.LOGFC_FILTER) & (res["adj.P.Val"] < config.ADJP_FILTER)]
    sig.to_csv(S1 / f"DEG_diff_{tag}.txt", sep="\t", index=False)
    up = int((sig["logFC"] > 0).sum()); dn = int((sig["logFC"] < 0).sum())
    print(f"[DEG {tag}] 유의 {len(sig)} (up {up} / down {dn})")
    return len(sig), up, dn


def main():
    mat = pd.read_csv(config.OUT_DIR / "GSE142025.labeled.txt", sep="\t", index_col=0)
    groups = [c.split("_")[-1] for c in mat.columns]
    mat = quantile_normalize(mat)
    n_lve, up, dn = run_contrast(mat, groups, config.GROUP_LATE,  config.GROUP_EARLY,   "Late_vs_Early")
    run_contrast(mat, groups, config.GROUP_LATE,  config.GROUP_CONTROL, "Late_vs_Control")
    run_contrast(mat, groups, config.GROUP_EARLY, config.GROUP_CONTROL, "Early_vs_Control")
    # 논문 대조 CSV (Late_vs_Early: 논문 2,833 / up1,557 / down1,276)
    pd.DataFrame({"metric": ["Late_vs_Early total", "up", "down"], "paper": [2833, 1557, 1276],
                  "ours": [n_lve, up, dn]}).to_csv(S1 / "compare_paper_vs_ours.GSE142025.csv", index=False)


if __name__ == "__main__":
    main()
