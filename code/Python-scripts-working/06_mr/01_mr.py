"""STEP6 : 2-표본 MR (R 01_mr.R 이식). IVW 주분석을 numpy 로 직접 구현.

노출 = FN1/ALDH2 cis-eQTL(GWAS-VCF, FORMAT ES:SE:LP:AF:SS:ID).
도구변수 = 저자 Supp Table 8(이미 clump). exposure 통계는 vcf 에서 추출.
결과 = FinnGen R12(3-1, 비압축 tsv) → GWAS Catalog GCST(3-2, gz).
IVW: beta = Σ(bx·by/sey²)/Σ(bx²/sey²), se = √(1/Σ(bx²/sey²)), OR = exp(beta).
라이브러리 매핑: TwoSampleMR::mr(IVW) → numpy 직접. (LD clumping 은 저자 instrument 사용으로 대체.)
출력: RES_DIR/step6_mr/*.csv
"""
import sys, pathlib, gzip
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import numpy as np
import pandas as pd
from math import erf, sqrt

OUT = config.RES_DIR / "step6_mr"; OUT.mkdir(parents=True, exist_ok=True)
FN1_VCF   = config.DATA_ROOT / "2-1. eqtl-a-ENSG00000115414.vcf"
ALDH2_VCF = config.DATA_ROOT / "2-2. eqtl-a-ENSG00000111275.vcf.gz"
FINNGEN   = config.DATA_ROOT / "3-1. finngen_R12_DM_NEPHROPATHY_EXMORE"
GCST      = config.DATA_ROOT / "3-2. GCST90435706" / "GCST90435706.tsv.gz"
SUPP8     = config.DIR_SUPP / "Supplementary Table 8. SNP characteristics and MR instrumental variables for candidate genes.csv"


def norm_sf(z):   # 2-sided p from |z|
    return 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))


def parse_vcf(path, gene):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt = f[8].split(":"); vals = dict(zip(fmt, f[9].split(":")))
            lp = float(vals["LP"])
            rows.append((vals["ID"], f[0], f[1], f[4], f[3],
                         float(vals["ES"]), float(vals["SE"]), 10 ** (-lp),
                         float(vals.get("AF", "nan")), gene))
    return pd.DataFrame(rows, columns=["SNP", "chr", "pos", "effect_allele", "other_allele",
                                       "beta", "se", "pval", "eaf", "exposure"])


def read_outcome(path, colmap, gzipped, snps):
    opener = gzip.open if gzipped else open
    usecols = list(colmap.keys())
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols,
                     compression="gzip" if gzipped else None)
    df = df.rename(columns=colmap)
    df = df[df["SNP"].isin(snps)]
    return df[["SNP", "beta", "se", "effect_allele", "other_allele", "pval"]]


def harmonise(expo, out):
    m = expo.merge(out, on="SNP", suffixes=(".exp", ".out"))
    by = []
    for _, r in m.iterrows():
        ea_e, oa_e = r["effect_allele.exp"].upper(), r["other_allele.exp"].upper()
        ea_o, oa_o = str(r["effect_allele.out"]).upper(), str(r["other_allele.out"]).upper()
        if (ea_o, oa_o) == (ea_e, oa_e):
            by.append(r["beta.out"])
        elif (ea_o, oa_o) == (oa_e, ea_e):
            by.append(-r["beta.out"])       # allele flip → 부호 반전
        else:
            by.append(np.nan)               # 불일치 → 제외
    m["by"] = by
    return m.dropna(subset=["by"])


# ---- MR 방법들 (TwoSampleMR 대응, numpy 직접구현) ----
def ivw(bx, by, sey):
    """IVW: beta = Σ(bx·by/sey²)/Σ(bx²/sey²), se = √(1/Σ(bx²/sey²))."""
    w = 1.0 / sey ** 2
    beta = np.sum(bx * by * w) / np.sum(bx ** 2 * w)
    se = sqrt(1.0 / np.sum(bx ** 2 * w))
    return beta, se, norm_sf(beta / se)


def mr_egger(bx, by, sey):
    """MR-Egger: 가중회귀 by~bx(절편). slope=인과추정, intercept=수평다면발현. weights=1/sey²."""
    s = np.sign(bx); s[s == 0] = 1
    bxo, byo = bx * s, by * s                       # bx>0 로 정렬(TwoSampleMR 관례)
    W = 1.0 / sey ** 2
    X = np.column_stack([np.ones_like(bxo), bxo])
    WX = X * W[:, None]
    XtWX = X.T @ WX
    coef = np.linalg.solve(XtWX, X.T @ (W * byo))   # [intercept, slope]
    resid = byo - X @ coef
    dof = max(len(bxo) - 2, 1)
    sigma2 = np.sum(W * resid ** 2) / dof
    cov = sigma2 * np.linalg.inv(XtWX)
    se_int, se_slope = sqrt(cov[0, 0]), sqrt(cov[1, 1])
    return coef[1], se_slope, norm_sf(coef[1] / se_slope), coef[0], se_int, norm_sf(coef[0] / se_int)


def _wmedian(b_iv, w):
    order = np.argsort(b_iv); b = b_iv[order]; w = w[order]
    cw = (np.cumsum(w) - 0.5 * w) / np.sum(w)
    k = np.max(np.where(cw < 0.5))
    return b[k] + (b[k + 1] - b[k]) * (0.5 - cw[k]) / (cw[k + 1] - cw[k])


def weighted_median(bx, by, sex, sey, nboot=1000):
    """Weighted median (TwoSampleMR): b_iv=by/bx, weights=bx²/sey². SE=부트스트랩."""
    b_iv = by / bx
    w = bx ** 2 / sey ** 2
    est = _wmedian(b_iv, w)
    rng = np.random.default_rng(123)
    boots = []
    for _ in range(nboot):
        bxs = rng.normal(bx, sex); bys = rng.normal(by, sey)
        try:
            boots.append(_wmedian(bys / bxs, bxs ** 2 / sey ** 2))
        except Exception:
            pass
    se = np.std(boots) if boots else np.nan
    return est, se, (norm_sf(est / se) if se and not np.isnan(se) else np.nan)


def cochran_q(bx, by, sey):
    """IVW 이질성 Cochran's Q + p."""
    beta, _, _ = ivw(bx, by, sey)
    w = 1.0 / sey ** 2
    Q = np.sum(w * (by - beta * bx) ** 2)
    dof = len(bx) - 1
    from scipy.stats import chi2
    return Q, dof, float(chi2.sf(Q, dof))


def run(expo_all, outcome_df, tag):
    """harmonise → outcome 필터(p>5e-06) → 3방법 MR + 이질성/다면발현 + 민감도."""
    print(f"\n===== MR: {tag} =====")
    rows, het_rows, ple_rows, ss_rows, loo_rows = [], [], [], [], []
    for gene, sub in expo_all.groupby("exposure"):
        h = harmonise(sub, outcome_df)
        n0 = len(h)
        h = h[h["pval.out"] > 5e-06]                        # 원본과 동일: outcome 유의 SNP 제거
        print(f"  [{gene}] pval.outcome>5e-06 필터: {n0} -> {len(h)} SNP")
        if len(h) < 2:
            print(f"  [{gene}] 도구변수 부족({len(h)}) → 생략"); continue
        bx = h["beta.exp"].to_numpy(); by = h["by"].to_numpy()
        sex = h["se.exp"].to_numpy();  sey = h["se.out"].to_numpy()
        snpid = h["SNP"].to_numpy(); n = len(bx)

        b_iv, se_iv, p_iv = ivw(bx, by, sey)
        b_eg, se_eg, p_eg, ic, se_ic, p_ic = mr_egger(bx, by, sey)
        b_wm, se_wm, p_wm = weighted_median(bx, by, sex, sey)
        for meth, b, se, p in [("Inverse variance weighted", b_iv, se_iv, p_iv),
                               ("MR Egger", b_eg, se_eg, p_eg),
                               ("Weighted median", b_wm, se_wm, p_wm)]:
            rows.append({"exposure": gene, "method": meth, "nsnp": n, "b": b, "se": se,
                         "or": np.exp(b), "or_lci95": np.exp(b - 1.96 * se),
                         "or_uci95": np.exp(b + 1.96 * se), "pval": p})
        Q, dof, pQ = cochran_q(bx, by, sey)
        het_rows.append({"exposure": gene, "method": "IVW", "Q": Q, "Q_df": dof, "Q_pval": pQ})
        ple_rows.append({"exposure": gene, "egger_intercept": ic, "se": se_ic, "pval": p_ic})
        # single-SNP (Wald ratio) + leave-one-out (IVW)
        for i in range(n):
            ss_rows.append({"exposure": gene, "SNP": snpid[i], "b": by[i] / bx[i],
                            "se": sey[i] / abs(bx[i])})
            idx = [j for j in range(n) if j != i]
            bL, seL, pL = ivw(bx[idx], by[idx], sey[idx])
            loo_rows.append({"exposure": gene, "SNP_excluded": snpid[i], "b": bL, "se": seL, "pval": pL})

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / f"table.MRresult.{tag}.csv", index=False)
    pd.DataFrame(het_rows).to_csv(OUT / f"table.heterogeneity.{tag}.csv", index=False)
    pd.DataFrame(ple_rows).to_csv(OUT / f"table.pleiotropy.{tag}.csv", index=False)
    pd.DataFrame(ss_rows).to_csv(OUT / f"table.singleSNP.{tag}.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(OUT / f"table.leaveoneout.{tag}.csv", index=False)
    print(tab.to_string(index=False))
    return tab, pd.DataFrame(ple_rows)


def ivw_filter(mrTab, ple, tag):
    """원본 MR2.R: IVW p<0.05 & 3방법 OR방향 일치 & 다면발현 p>0.05."""
    keep = []
    for gene, gd in mrTab.groupby("exposure"):
        if len(gd) != 3:
            continue
        ivwp = gd.loc[gd["method"] == "Inverse variance weighted", "pval"].iloc[0]
        dir_ok = (gd["or"] > 1).all() or (gd["or"] < 1).all()
        if ivwp < 0.05 and dir_ok:
            keep.append(gene)
    if not ple.empty:
        ok = set(ple.loc[ple["pval"] > 0.05, "exposure"])
        keep = [g for g in keep if g in ok]
    if keep:
        mrTab[mrTab["exposure"].isin(keep)].to_csv(OUT / f"IVW.filter.{tag}.csv", index=False)
    print(f"[IVW필터 {tag}] 통과: {', '.join(keep) if keep else '(없음)'}")
    return keep


def main():
    if not FN1_VCF.exists():
        print("[MR] eQTL vcf 없음 → 데이터 확인 필요."); return
    expo = pd.concat([parse_vcf(FN1_VCF, "FN1"), parse_vcf(ALDH2_VCF, "ALDH2")], ignore_index=True)
    s8 = pd.read_csv(SUPP8, skiprows=1)
    keep = {g: set(s8.loc[s8["exposure"] == g, "SNP"]) for g in ["FN1", "ALDH2"]}
    expo = pd.concat([expo[(expo.exposure == g) & (expo.SNP.isin(keep[g]))] for g in keep])
    expo.to_csv(OUT / "exposure.FN1_ALDH2.csv", index=False)
    snps = set(expo["SNP"])

    finn = read_outcome(FINNGEN, {"rsids": "SNP", "beta": "beta", "sebeta": "se",
                                  "alt": "effect_allele", "ref": "other_allele", "pval": "pval"}, False, snps)
    mf, pf = run(expo, finn, "FinnGen")
    ivw_filter(mf, pf, "FinnGen")
    gcst = read_outcome(GCST, {"variant_id": "SNP", "beta": "beta", "standard_error": "se",
                               "effect_allele": "effect_allele", "other_allele": "other_allele",
                               "p_value": "pval"}, True, snps)
    mg, pg = run(expo, gcst, "GCST")
    ivw_filter(mg, pg, "GCST")

    # 저자 Supp9(IVW, GWAS Catalog) 대조 CSV
    try:
        s9 = pd.read_csv(config.DIR_SUPP / "Supplementary Table 9. Mendelian randomization results of candidate genes IVW, MR-Egger, and weighted median analyses.csv", skiprows=1)
        s9 = s9[(s9["method"] == "Inverse variance weighted") & (s9["exposure"].isin(["FN1", "ALDH2"]))]
        s9 = s9[["exposure", "nsnp", "b", "pval"]].rename(columns={"nsnp": "paper_nsnp", "b": "paper_b", "pval": "paper_pval"})
        our = mg[mg["method"] == "Inverse variance weighted"][["exposure", "nsnp", "b", "or", "pval"]]
        our = our.rename(columns={"nsnp": "ours_nsnp", "b": "ours_b", "pval": "ours_pval"})
        cmp = s9.merge(our, on="exposure", how="outer")
        cmp.to_csv(OUT / "compare_paper_vs_ours.MR_IVW.csv", index=False)
        print("\n[Supp9 대조]\n", cmp.to_string(index=False))
    except Exception as e:
        print(f"[Supp9 대조 생략] {e}")
    print("\n[STEP6] 완료 (기대: GCST FN1 OR~2.78 위험 / ALDH2 OR~0.67 보호, 저자 Supp9 일치)")


if __name__ == "__main__":
    main()
