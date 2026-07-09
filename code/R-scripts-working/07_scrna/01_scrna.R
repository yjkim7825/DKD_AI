# 14_step7_scrna.R ----------------------------------------------------------
# STEP 7 : 단일세포 RNA-seq (원본 scRNA.Seurat0.7.R 로직, Seurat 표준)
#   데이터 = GSE209781 (10x, 6샘플: NM01-03 = Control, DKD01-03 = DKD) — 원본과 동일 데이터셋.
#   흐름  : Read10X → QC → LogNormalize → HVG → ScaleData → PCA → Harmony 통합
#           → FindClusters → UMAP → 마커기반 세포주석 → FN1/ALDH2 세포유형별 발현.
#   핵심  : FN1/ALDH2 가 어느 신장 세포에서 높은지 (DotPlot/VlnPlot + 표).
#   출력  : results/step7_scrna/
#   ※ tar 은 scratchpad 로만 풀고 원본 미변경.
# ---------------------------------------------------------------------------

suppressMessages({ library(Seurat); library(harmony); library(dplyr); library(ggplot2); library(patchwork) })
library(here)
source(here::here("config.R"))
OUT <- file.path(RES_DIR, "step7_scrna"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SCR <- file.path(SCRATCH_DIR, "sc")
dir.create(SCR, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

RAW <- file.path(DATA_ROOT, "1-7. GSE209781_RAW")
tars <- list.files(RAW, pattern = "\\.tar\\.gz$", full.names = TRUE)
samp_of <- function(f) sub("^GSM[0-9]+_", "", sub("\\.tar\\.gz$", "", basename(f)))

## ---- 1) tar 풀기 + Read10X + Seurat 객체 ----
objs <- list()
for (f in tars) {
  s <- samp_of(f)                    # NM01..DKD03
  d <- file.path(SCR, s)
  if (!dir.exists(file.path(d, s))) untar(f, exdir = d)
  tenx <- file.path(d, s)            # tar 내부 폴더명 = 샘플명
  mtx <- Read10X(data.dir = tenx)
  o <- CreateSeuratObject(counts = mtx, project = s, min.cells = 3, min.features = 200)
  o$orig.ident <- s
  o$group <- ifelse(grepl("^NM", s), "Control", "DKD")
  objs[[s]] <- o
  message("[read] ", s, ": ", ncol(o), " cells")
}

## ---- 2) merge + QC ----
sc <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
sc <- subset(sc, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)
message("[QC] 필터 후 세포: ", ncol(sc), " | 샘플별: ",
        paste(names(table(sc$orig.ident)), table(sc$orig.ident), sep="=", collapse=", "))

## ---- 3) 표준 전처리 + Harmony 통합 (체크포인트: 무거운 PCA/Harmony 1회만) ----
# 원본 scRNA.Seurat0.7.R 와 동일: LogNormalize(1e4) → vst 2000 → ScaleData → PCA30 → Harmony → dims 1:30.
# 원본은 Harmony 변수 "patient"(Control/DKD 2군), FindClusters resolution=0.2. 우리는:
#   · Harmony 변수 = orig.ident(6샘플 단위) — 표준 배치보정(원본 group단위보다 세분, DIFF 기록)
#   · resolution = 0.2 로 원본과 정렬.
RES  <- 0.2                                   # 원본과 동일
PCSEL <- 30
CKPT <- file.path(SCR, "GSE209781_harmony.rds")   # Harmony 까지만 캐시(클러스터링 전)
if (file.exists(CKPT)) {
  sc <- readRDS(CKPT); message("[ckpt] Harmony 전처리 결과 로드")
} else {
  sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 1e4)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
  sc <- ScaleData(sc, verbose = FALSE)
  sc <- RunPCA(sc, npcs = PCSEL, verbose = FALSE)
  sc <- RunHarmony(sc, group.by.vars = "orig.ident", verbose = FALSE)
  sc <- JoinLayers(sc)             # v5: 마커/발현 계산 위해 레이어 병합
  saveRDS(sc, CKPT)
}
# 클러스터링/UMAP (원본 resolution=0.2)
sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:PCSEL, verbose = FALSE)
sc <- FindClusters(sc, resolution = RES, verbose = FALSE)
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:PCSEL, verbose = FALSE)
message("[cluster] resolution=", RES, " 클러스터 수: ", length(levels(sc$seurat_clusters)))

## ---- 4) 마커기반 세포주석 (클러스터별 마커세트 평균발현 최대) ----
markers <- list(
  Podocyte    = c("NPHS1","NPHS2","PODXL","PTPRO"),
  PCT         = c("LRP2","CUBN","SLC34A1","SLC22A6","ALDOB","GPX3"),
  LOH         = c("UMOD","SLC12A1","CLDN16","KCNJ1"),
  DCT         = c("SLC12A3","WNK1","CALB1"),
  CD          = c("AQP2","AQP3","SCNN1G","SCNN1B"),
  Endothelial = c("PECAM1","FLT1","EMCN","VWF","PLVAP"),
  Mesangial   = c("PDGFRB","ACTA2","RGS5","NOTCH3","TAGLN"),
  Fibroblast  = c("COL1A1","COL3A1","DCN"),
  T_cell      = c("CD3E","TRAC","CD3D","IL7R"),
  NK          = c("GNLY","NKG7","GZMB","KLRD1"),
  Mono_Mac    = c("CD14","CD163","LYZ","C1QC","CSF1R"),
  B_cell      = c("CD79A","CD79B","MS4A1","BANK1"),
  Plasma      = c("JCHAIN","CD38","IGHG1","MZB1"),
  Mast        = c("TPSAB1","CPA3","KIT"),
  Neutrophil  = c("S100A8","S100A9","FCGR3B","CSF3R"))

avgByCluster <- AverageExpression(sc, assays = "RNA", slot = "data",
                                  group.by = "seurat_clusters")$RNA
colnames(avgByCluster) <- sub("^g", "", colnames(avgByCluster))   # v5 가 숫자 클러스터에 'g' 접두사 부여 → 제거
score <- sapply(markers, function(gs) {
  gs <- intersect(gs, rownames(avgByCluster))
  if (length(gs) == 0) return(rep(0, ncol(avgByCluster)))
  colMeans(avgByCluster[gs, , drop = FALSE])
})                                  # 행=클러스터, 열=세포타입
clAssign <- colnames(score)[max.col(score, ties.method = "first")]
names(clAssign) <- rownames(score)  # 클러스터ID -> 세포타입
sc$celltype <- unname(clAssign[as.character(sc$seurat_clusters)])  # unname: 위치기반 대입(바코드 매칭 회피)
message("[annot] 세포타입 분포:\n")
print(table(sc$celltype))
write.csv(data.frame(cluster = names(clAssign), celltype = clAssign,
                     n_cells = as.integer(table(sc$seurat_clusters)[names(clAssign)])),
          file.path(OUT, "cluster_annotation.csv"), row.names = FALSE)
# 마커 점수 행렬(클러스터 x 세포타입) 저장 — 주석 근거
write.csv(cbind(cluster = rownames(score), round(as.data.frame(score), 3)),
          file.path(OUT, "cluster_marker_scores.csv"), row.names = FALSE)

# FindAllMarkers (원본 scRNA.Seurat0.7.R 와 동일 인자) → 클러스터별 top 마커
allmk <- tryCatch(FindAllMarkers(sc, min.pct = 0.25, logfc.threshold = 0.5, only.pos = TRUE, verbose = FALSE),
                  error = function(e) { message("[FindAllMarkers 실패] ", conditionMessage(e)); NULL })
if (!is.null(allmk) && nrow(allmk)) {
  topmk <- allmk %>% dplyr::filter(p_val_adj < 0.05) %>% dplyr::group_by(cluster) %>%
    dplyr::slice_max(avg_log2FC, n = 20)
  write.csv(topmk, file.path(OUT, "cluster_markers_top20.csv"), row.names = FALSE)
  message("[markers] 클러스터별 top20 마커 저장")
}

## ---- 5) 핵심: FN1/ALDH2 세포유형별 발현 ----
genes <- c("FN1","ALDH2")
present <- genes[genes %in% rownames(sc)]
expr <- FetchData(sc, vars = c(present, "celltype"))
tab <- do.call(rbind, lapply(present, function(g) {
  agg <- tapply(expr[[g]], expr$celltype, function(v) c(mean = mean(v), pct = mean(v > 0)*100))
  data.frame(gene = g, celltype = names(agg),
             mean_expr = round(sapply(agg, `[`, "mean"), 3),
             pct_expressing = round(sapply(agg, `[`, "pct"), 1), row.names = NULL)
}))
tab <- tab[order(tab$gene, -tab$mean_expr), ]
write.csv(tab, file.path(OUT, "FN1_ALDH2_by_celltype.csv"), row.names = FALSE)
cat("\n== FN1/ALDH2 세포유형별 발현 (mean, %expressing) ==\n"); print(tab, row.names = FALSE)
for (g in present) {
  top <- tab[tab$gene == g, ][1, ]
  cat(sprintf(">> %s 최고발현 세포유형: %s (mean=%.3f, %%expr=%.1f)\n", g, top$celltype, top$mean_expr, top$pct_expressing))
}

## 논문 대조 CSV: 논문 scRNA — FN1 = 내피/系膜/손상PCT, ALDH2 = PCT(근위세뇨관)
paper_loc <- c(FN1 = "Endothelial/Mesangial/injured-PCT", ALDH2 = "PCT")
cmp7 <- do.call(rbind, lapply(present, function(g) {
  top <- tab[tab$gene == g, ][1, ]
  data.frame(gene = g, paper_top_celltype = paper_loc[[g]],
             ours_top_celltype = top$celltype,
             ours_mean = top$mean_expr, ours_pct = top$pct_expressing,
             consistent = grepl(top$celltype, paper_loc[[g]], ignore.case = TRUE) |
                          (g == "ALDH2" & top$celltype == "PCT") |
                          (g == "FN1" & top$celltype %in% c("Endothelial","Mesangial","PCT")),
             row.names = NULL)
}))
write.csv(cmp7, file.path(OUT, "compare_paper_vs_ours.scRNA.csv"), row.names = FALSE)
cat("\n== 논문 대조 (FN1/ALDH2 세포국소화) ==\n"); print(cmp7, row.names = FALSE)

## ---- 6) 그림 ----
Idents(sc) <- "celltype"
pdf(file.path(OUT, "UMAP_celltype.pdf"), width = 8, height = 6)
print(DimPlot(sc, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend())
dev.off()
pdf(file.path(OUT, "DotPlot_FN1_ALDH2.pdf"), width = 8, height = 5)
print(DotPlot(sc, features = present, cols = c("#ffffbf","#d73027")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()
pdf(file.path(OUT, "Violin_FN1_ALDH2.pdf"), width = 11, height = 5)
print(VlnPlot(sc, features = present, group.by = "celltype", pt.size = 0) & NoLegend() & labs(x = ""))
dev.off()
pdf(file.path(OUT, "FeaturePlot_FN1_ALDH2.pdf"), width = 11, height = 5)
print(FeaturePlot(sc, features = present, reduction = "umap"))
dev.off()

saveRDS(sc, file.path(SCR, "GSE209781_seurat.rds"))
message("\n[STEP7] 완료 -> ", OUT)
