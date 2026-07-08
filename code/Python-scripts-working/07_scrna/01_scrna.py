"""STEP7 : 단일세포 RNA-seq (R 01_scrna.R 이식). Seurat → scanpy.

데이터 = GSE209781 (10x, NM01-03=Control / DKD01-03=DKD).
흐름 = Read10x → QC → normalize → HVG → PCA → Harmony 통합 → Leiden 클러스터 →
       마커기반 세포주석 → FN1/ALDH2 세포유형별 발현.
라이브러리 매핑: Seurat→scanpy, RunHarmony→scanpy.external.pp.harmony_integrate,
   FindClusters→sc.tl.leiden, DotPlot/Violin→sc.pl.
⚠️ 무거움: 전체 재실행 금지 대상. py_compile 문법검증만; 실행은 사용자 판단.
출력: RES_DIR/step7_scrna/
"""
import sys, pathlib, tarfile
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import config

OUT = config.RES_DIR / "step7_scrna"; OUT.mkdir(parents=True, exist_ok=True)
RAW = config.DIR_GSE209781

MARKERS = {
    "Podocyte": ["NPHS1", "NPHS2", "PODXL", "PTPRO"],
    "PCT": ["LRP2", "CUBN", "SLC34A1", "SLC22A6", "ALDOB", "GPX3"],
    "LOH": ["UMOD", "SLC12A1", "CLDN16", "KCNJ1"],
    "DCT": ["SLC12A3", "WNK1", "CALB1"],
    "CD": ["AQP2", "AQP3", "SCNN1G", "SCNN1B"],
    "Endothelial": ["PECAM1", "FLT1", "EMCN", "VWF", "PLVAP"],
    "Mesangial": ["PDGFRB", "ACTA2", "RGS5", "NOTCH3", "TAGLN"],
    "Fibroblast": ["COL1A1", "COL3A1", "DCN"],
    "T_cell": ["CD3E", "TRAC", "CD3D", "IL7R"],
    "NK": ["GNLY", "NKG7", "GZMB", "KLRD1"],
    "Mono_Mac": ["CD14", "CD163", "LYZ", "C1QC", "CSF1R"],
    "B_cell": ["CD79A", "CD79B", "MS4A1", "BANK1"],
    "Plasma": ["JCHAIN", "CD38", "IGHG1", "MZB1"],
    "Mast": ["TPSAB1", "CPA3", "KIT"],
    "Neutrophil": ["S100A8", "S100A9", "FCGR3B", "CSF3R"],
}


def main():
    try:
        import scanpy as sc
        import anndata as ad
        import numpy as np
        import pandas as pd
    except Exception:
        print("[scRNA] TODO: scanpy 미설치 → 'pip install scanpy harmonypy leidenalg' 후 실행.")
        return

    scratch = config.SCRATCH_DIR / "sc"; scratch.mkdir(parents=True, exist_ok=True)
    adatas = []
    for tar in sorted(RAW.glob("*.tar.gz")):
        samp = tar.name.split("_", 1)[1].replace(".tar.gz", "")
        d = scratch / samp
        if not (d / samp).exists():
            with tarfile.open(tar) as t:
                t.extractall(d)
        a = sc.read_10x_mtx(d / samp)
        a.var_names_make_unique()
        a.obs["orig.ident"] = samp
        a.obs["group"] = config.GROUP_CONTROL if samp.startswith("NM") else config.GROUP_DKD
        adatas.append(a)
    # anndata 0.13: AnnData.concatenate 제거됨 → anndata.concat 사용(배치별 바코드 유일화)
    if len(adatas) > 1:
        adata = ad.concat(adatas, join="inner", label="sample",
                          keys=[a.obs["orig.ident"].iloc[0] for a in adatas], index_unique="-")
    else:
        adata = adatas[0]

    # QC (R: nFeature 300~5000, mt<10)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    adata = adata[(adata.obs.n_genes_by_counts > 300) & (adata.obs.n_genes_by_counts < 5000)
                  & (adata.obs.pct_counts_mt < 10)].copy()
    print(f"[QC] 세포 {adata.n_obs}")

    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10); sc.tl.pca(adata, n_comps=30)
    # Harmony 통합 — scanpy wrapper 가 버전차로 shape 오류 → harmonypy 직접 호출(Z_corr.T)
    try:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["orig.ident"])
        adata.obsm["X_pca_harmony"] = ho.Z_corr.T
        rep = "X_pca_harmony"
    except Exception as e:
        print(f"[harmony] 통합 실패 → PCA 로 진행(부분실행): {e}")
        rep = "X_pca"
    sc.pp.neighbors(adata, use_rep=rep, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)

    # 마커기반 자동주석 (클러스터별 마커세트 평균발현 최대)
    adraw = adata.raw.to_adata()   # 정규화·전체유전자(스케일 전) — 마커/발현 조회용
    scores = {}
    for ct, gs in MARKERS.items():
        gs = [g for g in gs if g in adraw.var_names]
        if not gs:
            continue
        sub = adraw[:, gs].X
        m = np.asarray(sub.mean(axis=1)).ravel()
        scores[ct] = pd.Series(m, index=adraw.obs_names).groupby(adata.obs["leiden"].values).mean()
    smat = pd.DataFrame(scores)
    cl2ct = smat.idxmax(axis=1)
    adata.obs["celltype"] = adata.obs["leiden"].map(cl2ct).astype(str)

    # 핵심: FN1/ALDH2 세포유형별 발현
    rows = []
    for g in ["FN1", "ALDH2"]:
        if g not in adraw.var_names:
            continue
        Xg = adraw[:, g].X
        expr = np.asarray(Xg.todense()).ravel() if hasattr(Xg, "todense") else np.asarray(Xg).ravel()
        s = pd.Series(expr, index=adraw.obs_names)
        for ct, idx in adata.obs.groupby("celltype").groups.items():
            v = s.loc[idx]
            rows.append({"gene": g, "celltype": ct, "mean_expr": float(v.mean()),
                         "pct_expressing": float((v > 0).mean() * 100)})
    tab = pd.DataFrame(rows).sort_values(["gene", "mean_expr"], ascending=[True, False])
    tab.to_csv(OUT / "FN1_ALDH2_by_celltype.csv", index=False)
    print(tab.to_string(index=False))
    for g in ["FN1", "ALDH2"]:
        top = tab[tab.gene == g].iloc[0]
        print(f">> {g} 최고발현: {top.celltype} (mean={top.mean_expr:.3f}, %={top.pct_expressing:.1f})")

    sc.pl.dotplot(adata, ["FN1", "ALDH2"], groupby="celltype", save="_FN1_ALDH2.pdf", show=False)


if __name__ == "__main__":
    main()
