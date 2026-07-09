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
S1  = config.RES_DIR / "step1_deg"                    # 원본과 동일 diff_ 입력 위치
HALLMARK_GMT = config.DATA_ROOT / "4-1. h.all.v2026.1.Hs.symbols.gmt"
KEGG_GMT     = config.DATA_ROOT / "4-2. c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt"
FOCUS = "EPITHELIAL_MESENCHYMAL|REACTIVE_OXYGEN|OXIDATIVE|INFLAMMATORY|TGF_BETA|APOPTOSIS"
KFOCUS = "ECM_RECEPTOR|CITRATE_CYCLE|OXIDATIVE_PHOSPHORYLATION|FOCAL_ADHESION|TRYPTOPHAN"
TAGS = ["Late_vs_Early", "Late_vs_Control", "Early_vs_Control"]


def _prerank(rnk, gmt):
    import gseapy as gp
    return gp.prerank(rnk=rnk, gene_sets=str(gmt), min_size=5, max_size=1000,
                      permutation_num=1000, seed=123, outdir=None, no_plot=True).res2d


def _run_gsea(gmt, prefix, focus):
    """원본 GSEA.R/hallmark.gsea.R 와 동일: diff_(유의 DEG) logFC 랭킹으로 GSEA."""
    try:
        import gseapy as gp  # noqa
    except Exception:
        print(f"[{prefix}] TODO: gseapy 미설치 → 'pip install gseapy'."); return
    for tag in TAGS:
        f = S1 / f"DEG_diff_{tag}.txt"                # 원본과 동일 diff_ 입력
        if not f.exists():
            print(f"[{prefix} {tag}] diff DEG 없음 → STEP1 02_deg 먼저: {f}"); continue
        deg = pd.read_csv(f, sep="\t")
        rnk = deg[["id", "logFC"]].dropna().drop_duplicates("id").sort_values("logFC", ascending=False)
        tab = _prerank(rnk, gmt)
        tab.to_csv(OUT / f"{prefix}.GSEA.{tag}.txt", sep="\t", index=False)
        term_col = "Term" if "Term" in tab.columns else tab.columns[0]
        foc = tab[tab[term_col].str.contains(focus, regex=True, na=False)]
        print(f"[{prefix} {tag}] ranked {len(rnk)} | 총 {len(tab)} | 관련 {len(foc)}행")


def run_hallmark():
    _run_gsea(HALLMARK_GMT, "Hallmark", FOCUS)


def run_kegg_gsea():
    """원본 GSEA.R 대응(신규): 동일 diff_ DEG 를 KEGG gmt 로 GSEA."""
    _run_gsea(KEGG_GMT, "KEGG", KFOCUS)


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
    print("[KEGG ssGSEA] 저장. FN1/ALDH2 ↔ 경로 상관 계산 중...")

    # --- FN1/ALDH2 ↔ 경로 Spearman 상관 (R 대조: FN1↔ECM 0.919, ALDH2↔TCA 0.921) ---
    from scipy.stats import spearmanr
    mat = scores.pivot(index="Term", columns="Name", values="NES")   # 경로 × 샘플
    samples = [s for s in mat.columns if s in expr.columns]; mat = mat[samples]
    rows = []
    for g in ["FN1", "ALDH2"]:
        if g not in expr.index:
            continue
        gv = expr.loc[g, samples].astype(float).values
        cs = pd.Series({p: spearmanr(gv, mat.loc[p, samples].astype(float).values).correlation
                        for p in mat.index}).sort_values(ascending=False)
        for p, r in cs.head(5).items():
            rows.append({"gene": g, "pathway": p, "rho": round(float(r), 3), "kind": "top5"})
        # R 대조용 특정 경로 명시(top5 밖이어도 기록)
        for p in ["KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_CITRATE_CYCLE_TCA_CYCLE"]:
            if p in cs.index:
                rows.append({"gene": g, "pathway": p, "rho": round(float(cs[p]), 3), "kind": "R대조"})
        print(f"  [{g}] top: {cs.index[0]} rho={cs.iloc[0]:.3f}")
    pd.DataFrame(rows).to_csv(OUT / "FN1_ALDH2_pathway_corr.csv", index=False)


def main():
    run_hallmark()
    run_kegg_gsea()       # 원본 GSEA.R 대응(신규)
    run_kegg_ssgsea()     # 논문 Fig2H,I


if __name__ == "__main__":
    main()
