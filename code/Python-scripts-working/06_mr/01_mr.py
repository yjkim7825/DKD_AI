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


def ivw(bx, by, sey):
    w = 1.0 / sey ** 2
    beta = np.sum(bx * by * w) / np.sum(bx ** 2 * w)
    se = sqrt(1.0 / np.sum(bx ** 2 * w))
    z = beta / se
    return beta, se, norm_sf(z), np.exp(beta), int(len(bx))


def run(expo_all, outcome_df, tag):
    print(f"\n===== MR: {tag} =====")
    res = []
    for gene, sub in expo_all.groupby("exposure"):
        h = harmonise(sub, outcome_df)
        if len(h) < 2:
            print(f"  [{gene}] 도구변수 부족({len(h)}) → 생략"); continue
        beta, se, p, orr, nsnp = ivw(h["beta.exp"].to_numpy(), h["by"].to_numpy(), h["se.out"].to_numpy())
        res.append({"exposure": gene, "nsnp": nsnp, "b": beta, "se": se,
                    "or": orr, "or_lci95": np.exp(beta - 1.96 * se), "or_uci95": np.exp(beta + 1.96 * se), "pval": p})
    tab = pd.DataFrame(res)
    tab.to_csv(OUT / f"table.MRresult.{tag}.csv", index=False)
    print(tab.to_string(index=False))
    return tab


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
    run(expo, finn, "FinnGen")
    gcst = read_outcome(GCST, {"variant_id": "SNP", "beta": "beta", "standard_error": "se",
                               "effect_allele": "effect_allele", "other_allele": "other_allele",
                               "p_value": "pval"}, True, snps)
    run(expo, gcst, "GCST")
    print("\n[STEP6] 완료 (기대: GCST 에서 FN1 OR~2.78 위험 / ALDH2 OR~0.67 보호, 저자 Supp9 일치)")


if __name__ == "__main__":
    main()
