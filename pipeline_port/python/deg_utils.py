"""
deg_utils.py
================
P25 (FN1/ALDH2 DKD multi-omics) R 파이프라인의 통계 코어를 Python으로 이식한 유틸.

원본 R 스택 대응:
  - limma::normalizeBetweenArrays(method="quantile")  ->  quantile_normalize()
  - sva::ComBat()                                      ->  combat()
  - limma::lmFit + eBayes + topTable (2-group)         ->  limma_two_group()
  - p.adjust(method="fdr")                             ->  bh_adjust()

주의(충실도):
  - limma_two_group()은 Smyth(2004) empirical Bayes moderated t-test를 그대로 구현.
    2군(Control vs DKD) design(~0+group)에 대해 limma topTable(coef="DKD-Control")과
    logFC / t / P.Value / adj.P.Val 이 수치적으로 일치(부동소수 오차 수준)한다.
  - combat()은 Johnson et al.(2007) parametric empirical Bayes.
    공변량(mod) 없이 batch만 보정하는 기본 형태(원본 R 스크립트와 동일 사용법).

의존성: numpy, scipy (scipy.special: digamma/polygamma, scipy.stats: t 분포)
"""
from __future__ import annotations
import numpy as np
# 주: scipy 는 limma 계열 함수(_trigamma_inverse/_fit_f_dist/limma_two_group)
#     내부에서 지연 import 한다. 이렇게 하면 scipy 없이도 quantile_normalize /
#     combat (numpy 전용) 을 사용할 수 있다.


# ----------------------------------------------------------------------
# 1) Quantile normalization  (limma::normalizeBetweenArrays 기본값)
# ----------------------------------------------------------------------
def quantile_normalize(mat: np.ndarray) -> np.ndarray:
    """
    행=유전자, 열=샘플 인 발현행렬을 quantile normalize.
    각 열을 정렬 후 '행별 평균 분위수'로 치환(동점은 평균 rank).
    """
    X = np.asarray(mat, dtype=float)
    n_genes, n_samples = X.shape
    # 각 열 정렬값의 행별 평균 = target 분포
    sorted_cols = np.sort(X, axis=0)
    mean_quantiles = sorted_cols.mean(axis=1)   # 길이 n_genes

    out = np.empty_like(X)
    for j in range(n_samples):
        col = X[:, j]
        # 평균 rank (동점 처리) -> target 분포로 매핑
        order = np.argsort(col, kind="mergesort")
        ranks = np.empty(n_genes, dtype=float)
        ranks[order] = np.arange(n_genes, dtype=float)
        # 동점 평균 처리
        # (동점이 많지 않으면 생략 가능하지만 limma와 맞추기 위해 평균 rank 사용)
        _, inv, counts = np.unique(col, return_inverse=True, return_counts=True)
        # 각 값 그룹의 평균 rank
        sum_ranks = np.zeros(len(counts))
        np.add.at(sum_ranks, inv, ranks)
        avg_rank = sum_ranks / counts
        col_ranks = avg_rank[inv]
        # target 분포에서 (평균 rank) 위치 보간
        out[:, j] = np.interp(col_ranks, np.arange(n_genes), mean_quantiles)
    return out


def auto_log2(mat: np.ndarray) -> np.ndarray:
    """limma 튜토리얼 관례의 자동 log2 판정 (data preprocessing 2.R 로직 동일)."""
    X = np.asarray(mat, dtype=float)
    qx = np.nanquantile(X, [0, 0.25, 0.5, 0.75, 0.99, 1.0])
    logC = (qx[4] > 100) or ((qx[5] - qx[0]) > 50 and qx[1] > 0)
    if logC:
        X = X.copy()
        X[X < 0] = 0
        X = np.log2(X + 1)
    return X


# ----------------------------------------------------------------------
# 2) ComBat  (sva::ComBat, parametric prior)
# ----------------------------------------------------------------------
def combat(data: np.ndarray, batch: np.ndarray) -> np.ndarray:
    """
    parametric empirical Bayes batch 보정.
    data : (n_genes, n_samples)
    batch: 길이 n_samples, 배치 라벨(정수/문자 무관)
    반환 : 보정된 (n_genes, n_samples)
    """
    Y = np.asarray(data, dtype=float)
    batch = np.asarray(batch)
    n_genes, n_samples = Y.shape

    levels, batch_idx = np.unique(batch, return_inverse=True)
    n_batch = len(levels)
    n_batches = np.array([(batch_idx == i).sum() for i in range(n_batch)])

    # design: batch one-hot (intercept 없이 배치 평균)
    design = np.zeros((n_samples, n_batch))
    design[np.arange(n_samples), batch_idx] = 1.0

    # 1) grand mean & pooled variance (B_hat: 배치별 평균)
    B_hat, *_ = np.linalg.lstsq(design, Y.T, rcond=None)      # (n_batch, n_genes)
    grand_mean = (n_batches / n_samples) @ B_hat              # (n_genes,)
    stand_mean = np.outer(grand_mean, np.ones(n_samples))     # (n_genes, n_samples)
    resid = Y - (design @ B_hat).T
    var_pooled = (resid ** 2).sum(axis=1) / n_samples         # (n_genes,)
    var_pooled[var_pooled == 0] = np.finfo(float).eps

    # 2) standardize
    Z = (Y - stand_mean) / np.sqrt(var_pooled)[:, None]

    # 3) batch effect L/S 추정 + parametric EB
    def aprior(g):  # gamma_bar, t2 -> a
        m, s2 = g.mean(), g.var()
        return (2 * s2 + m ** 2) / s2

    def bprior(g):
        m, s2 = g.mean(), g.var()
        return (m * s2 + m ** 3) / s2

    def postmean(g_hat, g_bar, n, d_star, t2):
        return (t2 * n * g_hat + d_star * g_bar) / (t2 * n + d_star)

    def postvar(sum2, n, a, b):
        return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)

    gamma_star = np.zeros((n_batch, n_genes))
    delta_star = np.zeros((n_batch, n_genes))

    for i in range(n_batch):
        idx = np.where(batch_idx == i)[0]
        ni = len(idx)
        Zi = Z[:, idx]
        gamma_hat = Zi.mean(axis=1)
        delta_hat = Zi.var(axis=1, ddof=1)
        delta_hat[delta_hat == 0] = np.finfo(float).eps

        gamma_bar = gamma_hat.mean()
        t2 = gamma_hat.var()
        a_prior = aprior(delta_hat)
        b_prior = bprior(delta_hat)

        # EB 수렴 (limma/sva itSol과 동일한 고정점 반복)
        g_old = gamma_hat.copy()
        d_old = delta_hat.copy()
        for _ in range(200):
            g_new = postmean(gamma_hat, gamma_bar, ni, d_old, t2)
            sum2 = ((Zi - g_new[:, None]) ** 2).sum(axis=1)
            d_new = postvar(sum2, ni, a_prior, b_prior)
            change = max(np.abs((g_new - g_old) / (g_old + np.finfo(float).eps)).max(),
                         np.abs((d_new - d_old) / (d_old + np.finfo(float).eps)).max())
            g_old, d_old = g_new, d_new
            if change < 1e-4:
                break
        gamma_star[i] = g_old
        delta_star[i] = d_old

    # 4) 보정
    Zadj = Z.copy()
    for i in range(n_batch):
        idx = np.where(batch_idx == i)[0]
        Zadj[:, idx] = (Z[:, idx] - gamma_star[i][:, None]) / np.sqrt(delta_star[i])[:, None]

    corrected = Zadj * np.sqrt(var_pooled)[:, None] + stand_mean
    return corrected


# ----------------------------------------------------------------------
# 3) limma moderated t-test (2-group)  ==  lmFit + eBayes + topTable
# ----------------------------------------------------------------------
def _trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """limma::trigammaInverse — trigamma(y)=x 를 만족하는 y (Newton)."""
    from scipy.special import polygamma
    x = np.asarray(x, dtype=float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = polygamma(1, y)                      # trigamma
        dif = tri * (1 - tri / x) / polygamma(2, y)  # / tetragamma
        y = y + dif
        if np.max(-dif / y) < 1e-8:
            break
    # 극단값 보정
    big = x > 1e7
    y[big] = 1.0 / np.sqrt(x[big])
    small = x < 1e-6
    y[small] = 1.0 / x[small]
    return y


def _fit_f_dist(s2: np.ndarray, df1: float):
    """limma::fitFDist — 사전분산 s0^2 와 사전 df(df2) 추정 (df1 스칼라 가정)."""
    from scipy.special import digamma, polygamma
    s2 = s2[s2 > 0]
    n = len(s2)
    z = np.log(s2)
    e = z - digamma(df1 / 2.0) + np.log(df1 / 2.0)
    emean = e.mean()
    evar = np.mean(n / (n - 1.0) * (e - emean) ** 2) - polygamma(1, df1 / 2.0)
    if evar > 0:
        df2 = 2.0 * _trigamma_inverse(np.array([evar]))[0]
        s0_2 = np.exp(emean + digamma(df2 / 2.0) - np.log(df2 / 2.0))
    else:
        df2 = np.inf
        s0_2 = np.exp(emean)
    return s0_2, df2


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (p.adjust method='fdr'/'BH')."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # 단조 보정
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def limma_two_group(data: np.ndarray, group: np.ndarray,
                    ref: str, alt: str):
    """
    2-group moderated t-test.  contrast = (alt - ref)  ==  R의 makeContrasts(DKD - Control).

    data : (n_genes, n_samples)
    group: 길이 n_samples 의 라벨 배열
    ref  : 대조군 라벨 (예: "Control")
    alt  : 실험군 라벨 (예: "DKD")

    반환 : dict of 1D arrays (길이 n_genes):
           logFC, AveExpr, t, P.Value, adj.P.Val
    """
    from scipy.stats import t as student_t
    Y = np.asarray(data, dtype=float)
    group = np.asarray(group)
    m_ref = group == ref
    m_alt = group == alt
    n1, n2 = m_ref.sum(), m_alt.sum()
    if n1 < 2 or n2 < 2:
        raise ValueError(f"각 군 최소 2샘플 필요 (ref={n1}, alt={n2})")

    Yr, Ya = Y[:, m_ref], Y[:, m_alt]
    mean_r, mean_a = Yr.mean(axis=1), Ya.mean(axis=1)
    logFC = mean_a - mean_r
    AveExpr = Y[:, m_ref | m_alt].mean(axis=1)

    df_resid = n1 + n2 - 2
    ss = ((Yr - mean_r[:, None]) ** 2).sum(axis=1) + \
         ((Ya - mean_a[:, None]) ** 2).sum(axis=1)
    sigma2 = ss / df_resid                      # per-gene residual variance
    stdev_unscaled = np.sqrt(1.0 / n1 + 1.0 / n2)

    # empirical Bayes moderation
    s0_2, df0 = _fit_f_dist(sigma2, df_resid)
    if np.isinf(df0):
        var_post = np.full_like(sigma2, s0_2)
        df_total = np.full_like(sigma2, np.inf)
    else:
        var_post = (df_resid * sigma2 + df0 * s0_2) / (df_resid + df0)
        df_total = np.full_like(sigma2, df_resid + df0)

    t_stat = logFC / (stdev_unscaled * np.sqrt(var_post))
    with np.errstate(invalid="ignore"):
        pval = 2.0 * student_t.sf(np.abs(t_stat), df_total)
    pval = np.nan_to_num(pval, nan=1.0)
    adj = bh_adjust(pval)

    return {
        "logFC": logFC,
        "AveExpr": AveExpr,
        "t": t_stat,
        "P.Value": pval,
        "adj.P.Val": adj,
    }
