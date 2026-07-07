"""STEP5 : Hallmark GSEA + KEGG ssGSEA (R 01_gsea.R 이식).

(1) Hallmark GSEA — GSE142025 DEG 3비교의 logFC 랭킹, gseapy.prerank (gmt 4-1).
(2) KEGG ssGSEA — GSE96804 발현, gseapy.ssgsea (gmt 4-2).
FN1(EMT/산화)·ALDH2(산화/대사) 관련 경로 유의성 확인.
라이브러리 매핑: clusterProfiler::GSEA → gseapy.prerank, GSVA ssGSEA → gseapy.ssgsea.
출력: RES_DIR/step5_gsea/
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config
import pandas as pd

OUT = config.RES_DIR / "step5_gsea"; OUT.mkdir(parents=True, exist_ok=True)
HALLMARK_GMT = config.DATA_ROOT / "4-1. h.all.v2026.1.Hs.symbols.gmt"
KEGG_GMT     = config.DATA_ROOT / "4-2. c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
FOCUS = "EPITHELIAL_MESENCHYMAL|REACTIVE_OXYGEN|OXIDATIVE|INFLAMMATORY|TGF_BETA|APOPTOSIS"


def run_hallmark():
    try:
        import gseapy as gp
    except Exception:
        print("[Hallmark] TODO: gseapy 미설치 → 'pip install gseapy' 후 실행."); return
    for tag in ["Late_vs_Early", "Late_vs_Control", "Early_vs_Control"]:
        f = config.RES_DIR / f"DEG_all_{tag}.txt"
        if not f.exists():
            print(f"[Hallmark {tag}] DEG 없음 → STEP1 02_deg 먼저: {f}"); continue
        deg = pd.read_csv(f, sep="\t")
        rnk = deg[["id", "logFC"]].dropna().sort_values("logFC", ascending=False)
        res = gp.prerank(rnk=rnk, gene_sets=str(HALLMARK_GMT), min_size=5, max_size=1000,
                         permutation_num=1000, seed=123, outdir=None, no_plot=True)
        tab = res.res2d
        tab.to_csv(OUT / f"Hallmark.GSEA.{tag}.txt", sep="\t", index=False)
        term_col = "Term" if "Term" in tab.columns else tab.columns[0]
        foc = tab[tab[term_col].str.contains(FOCUS, regex=True, na=False)]
        print(f"[Hallmark {tag}] 총 {len(tab)} | FN1/ALDH2 관련 {len(foc)}행 저장")


def run_kegg_ssgsea():
    try:
        import gseapy as gp
    except Exception:
        print("[KEGG ssGSEA] TODO: gseapy 미설치."); return
    expr_f = config.OUT_DIR / "GSE96804.labeled.txt"
    if not expr_f.exists():
        print("[KEGG ssGSEA] GSE96804 매트릭스 없음 → R STEP2 필요."); return
    expr = pd.read_csv(expr_f, sep="\t", index_col=0)
    ss = gp.ssgsea(data=expr, gene_sets=str(KEGG_GMT), outdir=None, sample_norm_method="rank", no_plot=True)
    scores = ss.res2d
    scores.to_csv(OUT / "KEGG.ssGSEA.scores.txt", sep="\t", index=False)
    print(f"[KEGG ssGSEA] 저장. (FN1↔ECM_RECEPTOR, ALDH2↔대사 경로 상관은 확장 분석 TODO)")


def main():
    run_hallmark()
    run_kegg_ssgsea()


if __name__ == "__main__":
    main()
