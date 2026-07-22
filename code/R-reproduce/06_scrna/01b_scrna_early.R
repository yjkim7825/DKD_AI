# ============================================================================
# 01b_scrna_early.R — 초기 DKD (GSE131882) 를 말기와 같은 파이프라인으로 처리
#   저자 코드는 말기(GSE209781)만 → 초기는 우리가 동일 파라미터로 추가 처리.
#   목적: Fig 6H / 단계별(Control·Early·Late) 비교용 Early 객체 생성.
#   입력 : data/1-6. GSE131882_RAW/*.dgecounts.rds.gz (control s1~3 + diabetes s1~3)
#   출력 : output/Early_annotated.RDS
#   논문 버전: Seurat v5.3.0, harmony v1.2.3
# ============================================================================
set.seed(123)
suppressMessages({ library(Seurat); library(harmony); library(dplyr); library(ggplot2) })

ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MLDIR <- file.path(ROOT, "R-reproduce/06_scrna")
RAW   <- file.path(ROOT, "data/1-6. GSE131882_RAW")
OUT   <- file.path(MLDIR, "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(MLDIR, "R/scrna_func.R"))

## ── 6샘플 읽어 병합 (dropEst dgecounts) ─────────────────────────────────────
samples <- list(
  Control_s1 = "GSM3823939_control.s1",  Control_s2 = "GSM3823940_control.s2",
  Control_s3 = "GSM3823941_control.s3",  Early_s1   = "GSM3823942_diabetes.s1",
  Early_s2   = "GSM3823943_diabetes.s2", Early_s3   = "GSM3823944_diabetes.s3")
# GSE131882은 유전자 이름이 Ensembl ID → GSE209781 features.tsv로 심볼 변환
ens2sym <- load_ens2sym(file.path(ROOT, "data/1-7. GSE209781_RAW/NM01/features.tsv.gz"))
objs <- list()
for (nm in names(samples)) {
  cm <- read_dge(file.path(RAW, paste0(samples[[nm]], ".dgecounts.rds.gz")), ens2sym = ens2sym)
  objs[[nm]] <- CreateSeuratObject(cm, project = nm, min.cells = 3, min.features = 200)
  message("[읽기] ", nm, " 완료 (", nrow(objs[[nm]]), "유전자 × ", ncol(objs[[nm]]), "세포)")
}
Early <- merge(objs[[1]], objs[-1], add.cell.ids = names(objs)); Early <- JoinLayers(Early)
Early$group   <- ifelse(grepl("Control", Early$orig.ident), "Control", "Early_DKD")
Early$patient <- Early$orig.ident

## ── QC → 정규화 → PCA → Harmony → 클러스터 (말기와 동일 파라미터) ────────────
Early[["percent.mt"]] <- PercentageFeatureSet(Early, pattern = "^MT-")
Early <- subset(Early, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10)  # mt<10 (재현일치)
Early <- NormalizeData(Early, normalization.method = "LogNormalize", scale.factor = 10000)
Early <- FindVariableFeatures(Early, selection.method = "vst", nfeatures = 2000)
Early <- ScaleData(Early)
Early <- RunPCA(Early, npcs = 30, verbose = FALSE)
Early <- RunHarmony(Early, "patient")
Early <- FindNeighbors(Early, reduction = "harmony", dims = 1:30)
Early <- FindClusters(Early, resolution = 0.2)
Early <- RunTSNE(Early, reduction = "harmony", dims = 1:30)
Early <- RunUMAP(Early, reduction = "harmony", dims = 1:30)

## ── 세포주석 (마커 자동배정) + dPCT 세분 ────────────────────────────────────
Early <- annotate_by_markers(Early, KIDNEY_MARKERS)
# PCT를 손상마커로 PCT/dPCT 세분 (논문 dPCT 정의)
inj <- colMeans(GetAssayData(Early, layer = "data")[intersect(c("VCAM1","HAVCR2","SPP1"), rownames(Early)), , drop = FALSE])
pctcell <- which(Early$celltype == "PCT")
if (length(pctcell)) {
  thr <- quantile(inj[pctcell], 0.7)
  ct <- as.character(Early$celltype); ct[pctcell[inj[pctcell] >= thr]] <- "dPCT"
  Early$celltype <- factor(ct); Idents(Early) <- Early$celltype
}
pdf(file.path(OUT, "01b.Early_tSNE_celltype.pdf"), width = 7, height = 6)
print(DimPlot(Early, reduction = "tsne", label = TRUE, repel = TRUE)); dev.off()

saveRDS(Early, file.path(OUT, "Early_annotated.RDS"))
cat("\n★ [01b] 초기 DKD(GSE131882) 처리 완료 → Early_annotated.RDS\n")
