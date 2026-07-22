# ============================================================================
# 01_scrna_process.R — scRNA 1/4: QC → 정규화 → PCA → Harmony → 클러스터링 → 세포주석
#   저자 "single-cell RNA-seq analysis 1.R" 재현 (DKD 파이프라인만, 잔여 코드 제외)
#   데이터: GSE209781 (신장, Control 3 + Late DKD 3 = 6샘플)
#   입력 : data/1-7. GSE209781_RAW/ (10x, 각 GSM tar.gz) 또는 병합 rda
#   출력 : output/Late_annotated.RDS (→ 02가 읽음), QC·UMAP·마커 그림
#   ※ Seurat 필요. 샌드박스 실행 불가 → 사용자 R에서 실행.
# ============================================================================
set.seed(123)
# 논문 버전: Seurat v5.3.0, harmony v1.2.3 (Cell Ranger 4 산출 10x)
suppressMessages({
  library(Seurat); library(harmony); library(dplyr)
  library(ggplot2); library(patchwork)
})

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
RAW  <- file.path(ROOT, "data/1-7. GSE209781_RAW")   # 각 GSM tar.gz → 압축해제 후 10x 폴더
MLDIR <- file.path(ROOT, "R-reproduce/06_scrna")
OUT  <- file.path(MLDIR, "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(MLDIR, "R/scrna_func.R"))   # annotate_by_markers, KIDNEY_MARKERS

## ── 0) 6샘플 읽어 병합 (10x) ────────────────────────────────────────────────
#   tar.gz 압축해제 시 각 폴더에 barcodes/features/matrix. 폴더명=GSM_라벨.
#   ※ 이미 병합된 late_merged_seurat.rda 있으면 아래 대신 readRDS 사용.
samples <- list(
  Control_NM01 = "NM01", Control_NM02 = "NM02", Control_NM03 = "NM03",
  Late_DKD01   = "DKD01", Late_DKD02   = "DKD02", Late_DKD03   = "DKD03")
# ★ 메모리 절약: 샘플 하나씩 읽어 즉시 QC로 축소 → 큰 원시행렬 바로 해제(rm+gc)
#   (raw 배럴 679만 → CreateSeuratObject+QC 후 수천 세포로 줄여 merge → RAM 폭발 방지)
MT_CUT <- 10   # 미토 컷 (논문 15% / 저자코드 10% — 10%라야 논문 세포수 18,816 재현)
objs <- list()
for (nm in names(samples)) {
  dir <- file.path(RAW, samples[[nm]])
  if (!dir.exists(dir)) { message("[!] 폴더 없음: ", dir, " (tar.gz 압축해제 필요)"); next }
  cnt <- Read10X(dir)
  o <- CreateSeuratObject(cnt, project = nm, min.cells = 3, min.features = 200)
  o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
  o <- subset(o, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < MT_CUT)
  objs[[nm]] <- o
  rm(cnt, o); gc()                                    # 큰 원시행렬 즉시 해제
  message("[읽기] ", nm, " 완료")
}
Late <- merge(objs[[1]], objs[-1], add.cell.ids = names(objs))
Late <- JoinLayers(Late)                              # Seurat v5 레이어 병합
rm(objs); gc()
Late$group   <- ifelse(grepl("^Control", Late$orig.ident), "Control", "Late_DKD")
Late$patient <- Late$orig.ident

## ── 1) QC 그림 (필터는 위 로딩 루프에서 이미 샘플별로 적용됨) ────────────────
#   ⚠ 미토 기준 논문 15% vs 저자코드 10% → 10%라야 논문 세포수(18,816) 재현 (MT_CUT는 위)
pdf(file.path(OUT, "01.featureViolin.pdf"), width = 10, height = 5)
print(VlnPlot(Late, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)); dev.off()
cat("[QC] 통과 세포:", ncol(Late), "개\n")

## ── 2) 정규화 · HVG · 스케일 · PCA ──────────────────────────────────────────
Late <- NormalizeData(Late, normalization.method = "LogNormalize", scale.factor = 10000)
Late <- FindVariableFeatures(Late, selection.method = "vst", nfeatures = 2000)
Late <- ScaleData(Late)
Late <- RunPCA(Late, npcs = 30, verbose = FALSE)

## ── 3) Harmony 배치보정 (저자: group.by="patient") ──────────────────────────
Late <- RunHarmony(Late, "patient", plot_convergence = TRUE)
pdf(file.path(OUT, "01.after_harmony.pdf")); print(DimPlot(Late, reduction = "harmony", group.by = "group")); dev.off()

## ── 4) 클러스터링 · UMAP/TSNE (저자: dims=1:30, resolution=0.2) ──────────────
Late <- FindNeighbors(Late, reduction = "harmony", dims = 1:30)
Late <- FindClusters(Late, resolution = 0.2)          # → 12 클러스터(0~11)
Late <- RunUMAP(Late, reduction = "harmony", dims = 1:30)
Late <- RunTSNE(Late, reduction = "harmony", dims = 1:30)
pdf(file.path(OUT, "01.UMAP_clusters.pdf")); print(DimPlot(Late, reduction = "umap", label = TRUE)); dev.off()

## ── 5) 세포주석 — 마커 자동배정 (클러스터 번호 무관, 생물학 기준) ────────────
#   ※ 저자는 수동 RenameIdents(번호 하드코딩)였으나, 우리 클러스터 번호가 저자와 달라
#     하드코딩은 라벨이 뒤바뀜(예: ALDH2 높은데 dPCT로 오분류). → 마커 z-score 자동배정.
markers <- FindAllMarkers(Late, min.pct = 0.25, logfc.threshold = 0.5, only.pos = TRUE)
write.csv(markers, file.path(OUT, "01.cluster_markers.csv"), row.names = FALSE)
Late <- annotate_by_markers(Late, KIDNEY_MARKERS)
# 근위세관(PCT)을 손상마커(VCAM1/HAVCR2/SPP1)로 PCT/dPCT 세분
inj <- colMeans(GetAssayData(Late, layer = "data")[intersect(c("VCAM1","HAVCR2","SPP1"), rownames(Late)), , drop = FALSE])
pctcell <- which(Late$celltype == "PCT")
if (length(pctcell)) {
  thr <- quantile(inj[pctcell], 0.7); ct <- as.character(Late$celltype)
  ct[pctcell[inj[pctcell] >= thr]] <- "dPCT"; Late$celltype <- factor(ct); Idents(Late) <- Late$celltype
}
cat("[주석] 세포타입:", paste(names(table(Late$celltype)), collapse = ", "), "\n")
pdf(file.path(OUT, "01.UMAP_celltype.pdf")); print(DimPlot(Late, reduction = "umap", label = TRUE)); dev.off()

## ── 6) 논문 Figure 7 패널 그림 ──────────────────────────────────────────────
# Fig 7A: t-SNE 세포 아틀라스
pdf(file.path(OUT, "01.Fig7A_tSNE_celltype.pdf"), width = 7, height = 6)
print(DimPlot(Late, reduction = "tsne", label = TRUE, repel = TRUE)); dev.off()

# Fig 7D: 12 세포타입 표준 마커 버블(DotPlot) — 논문 Fig 7D 마커 순서
fig7d_markers <- c("CD3E","TRAC","NKG7","ALDOB","LRP2","CUBN","PECAM1","EMCN",
                   "CD14","CD163","PDGFRB","ACTA2","RGS5","UMOD","SLC12A1","WNK1",
                   "PVALB","CD79B","MS4A1","TPSAB1","CPA3","S100A8","S100A9",
                   "CD38","JCHAIN","GZMB","PTPRC","ALDH2","FN1")
fig7d_markers <- fig7d_markers[fig7d_markers %in% rownames(Late)]
pdf(file.path(OUT, "01.Fig7D_marker_dotplot.pdf"), width = 12, height = 5)
print(DotPlot(Late, features = fig7d_markers,
              cols = c("#ffffbf", "#d73027"), dot.scale = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))); dev.off()

# FN1·ALDH2 세포별 발현 (Fig 7G 보조)
pdf(file.path(OUT, "01.ALDH2_FN1_violin.pdf"), width = 9, height = 4)
print(VlnPlot(Late, c("ALDH2", "FN1"), group.by = "celltype", pt.size = 0)); dev.off()

saveRDS(Late, file.path(OUT, "Late_annotated.RDS"))
cat("\n★ [01] 완료 → Late_annotated.RDS (12 세포타입 주석). 다음: 02_scrna_proportion.R\n")
