"""STEP6 MR — GCST90435706 outcome 만 실행 (FinnGen 2.1GB 는 무거워 스킵).
01_mr.py 의 함수를 재사용. 원본·R 무수정, 출력은 results/step6_mr/ 로만.
"""
import sys, pathlib, importlib.util
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import pandas as pd

_p = pathlib.Path(__file__).with_name("01_mr.py")
_spec = importlib.util.spec_from_file_location("mr", _p)
mr = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(mr)

expo = pd.concat([mr.parse_vcf(mr.FN1_VCF, "FN1"), mr.parse_vcf(mr.ALDH2_VCF, "ALDH2")], ignore_index=True)
s8 = pd.read_csv(mr.SUPP8, skiprows=1)
keep = {g: set(s8.loc[s8["exposure"] == g, "SNP"]) for g in ["FN1", "ALDH2"]}
expo = pd.concat([expo[(expo.exposure == g) & (expo.SNP.isin(keep[g]))] for g in keep])
expo.to_csv(mr.OUT / "exposure.FN1_ALDH2.csv", index=False)
print(f"[exposure] FN1 {int((expo.exposure=='FN1').sum())} SNP / ALDH2 {int((expo.exposure=='ALDH2').sum())} SNP")
snps = set(expo["SNP"])

gcst = mr.read_outcome(mr.GCST, {"variant_id": "SNP", "beta": "beta", "standard_error": "se",
                                 "effect_allele": "effect_allele", "other_allele": "other_allele",
                                 "p_value": "pval"}, True, snps)
print(f"[GCST] outcome SNP 매칭 {len(gcst)}")
mr.run(expo, gcst, "GCST")
print("\n[STEP6] GCST 완료 (FinnGen 스킵). 기대: FN1 b~1.021/OR~2.78, ALDH2 b~-0.395/OR~0.67")
