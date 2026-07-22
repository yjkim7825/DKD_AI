# ============================================================================
# 02_scrna_proportion.R — scRNA 2/4: 세포비율 + FN1/ALDH2 발현정량 + 유전자세트 스코어
#   저자 "single-cell RNA-seq analysis 2.R" 재현
#   입력 : output/Late_annotated.RDS (01 산출),  data/4-1 Hallmark gmt(EMT 등)
#   출력 : 세포비율 누적막대, FN1/ALDH2 바이올린·FeaturePlot(+wilcox), EMT 등 ModuleScore
#          → output/Late_with_geneset_scores.RDS (→ 03이 읽음)
# ============================================================================
set.seed(123)
suppressMessages({
  library(Seurat); library(ggplot2); library(ggpubr); library(dplyr)
})

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
OUT  <- file.path(ROOT, "R-reproduce/06_scrna/output")
Late <- readRDS(file.path(OUT, "Late_annotated.RDS"))
Idents(Late) <- Late$celltype

## ── 1) 세포비율 누적막대 (Fig 7B: celltype × group) ─────────────────────────
prop <- as.data.frame(prop.table(table(Late$celltype, Late$group), margin = 2))
colnames(prop) <- c("Celltype", "Group", "Proportion")
pdf(file.path(OUT, "02.Fig7B_cell_proportion.pdf"), width = 6, height = 5)
print(ggplot(prop, aes(Group, Proportion, fill = Celltype)) +
        geom_bar(stat = "identity", position = "fill") + theme_classic() +
        labs(title = "Cell type proportion (Control vs Late DKD)")); dev.off()

# Fig 7C: 근위세관 아형(PCT/dPCT)만 비율
pt <- subset(Late, idents = c("PCT", "dPCT"))
propc <- as.data.frame(prop.table(table(pt$celltype, pt$group), margin = 2))
colnames(propc) <- c("Celltype", "Group", "Proportion")
pdf(file.path(OUT, "02.Fig7C_PCT_subtype_proportion.pdf"), width = 4, height = 5)
print(ggplot(propc, aes(Group, Proportion, fill = Celltype)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = c(PCT = "#00A087", dPCT = "#E69F00")) +
        theme_classic() + labs(title = "Proximal tubule subtypes")); dev.off()

# Fig 7E: PCT 아형별 주요 기능유전자 (ALDOB·CUBN·SLC34A1) 바이올린
pdf(file.path(OUT, "02.Fig7E_PCT_function_genes.pdf"), width = 8, height = 3.5)
print(VlnPlot(pt, c("ALDOB", "CUBN", "SLC34A1"), group.by = "celltype", pt.size = 0, ncol = 3)); dev.off()

## ── 2) FN1·ALDH2 발현 시각화 + 통계검정 ─────────────────────────────────────
# 세포타입별 바이올린
pdf(file.path(OUT, "02.ALDH2_FN1_violin_celltype.pdf"), width = 9, height = 4)
print(VlnPlot(Late, c("ALDH2", "FN1"), group.by = "celltype", pt.size = 0)); dev.off()
# tSNE FeaturePlot (Control vs DKD 분할)
pdf(file.path(OUT, "02.FN1_featureplot.pdf"), width = 9, height = 4)
print(FeaturePlot(Late, "FN1", reduction = "tsne", split.by = "group", order = TRUE)); dev.off()
pdf(file.path(OUT, "02.ALDH2_featureplot.pdf"), width = 9, height = 4)
print(FeaturePlot(Late, "ALDH2", reduction = "tsne", split.by = "group", order = TRUE)); dev.off()

# PCT vs dPCT 한정 발현 비교 (저자: wilcox)
pdpct <- subset(Late, idents = c("PCT", "dPCT"))
for (g in c("FN1", "ALDH2")) {
  ex <- FetchData(pdpct, vars = c(g, "celltype"))
  p  <- wilcox.test(ex[[g]] ~ ex$celltype)$p.value
  cat(sprintf("[PCT vs dPCT] %s wilcox p = %.3g\n", g, p))
}
pdf(file.path(OUT, "02.Fig7G_ALDH2_FN1_violin_PCT_dPCT.pdf"), width = 7, height = 4)  # Fig 7G
print(VlnPlot(pdpct, c("ALDH2", "FN1"), group.by = "celltype", pt.size = 0)); dev.off()

## ── 3) 유전자세트 ModuleScore (EMT·염증·산화·자가포식·노화·세포사) ──────────
#   저자: AddModuleScore(ctrl=…, seed=123). gmt는 data/4-1 Hallmark 사용 가능.
gmt <- file.path(ROOT, "data/4-1. h.all.v2026.1.Hs.symbols.gmt")
read_gmt <- function(path, key) {
  ln <- readLines(path); hit <- grep(key, ln, value = TRUE)
  if (!length(hit)) return(NULL)
  strsplit(hit[1], "\t")[[1]][-c(1, 2)]
}
sets <- list(EMT = "EPITHELIAL_MESENCHYMAL_TRANSITION",
             Inflam = "INFLAMMATORY_RESPONSE",
             Apoptosis = "APOPTOSIS")
for (nm in names(sets)) {
  genes <- read_gmt(gmt, sets[[nm]])
  if (is.null(genes)) next
  Late <- AddModuleScore(Late, features = list(intersect(genes, rownames(Late))),
                         name = paste0(nm, "_"), ctrl = 50, seed = 123)
}
# EMT 스코어 세포타입별 박스플롯
if ("EMT_1" %in% colnames(Late@meta.data)) {
  df <- FetchData(Late, vars = c("EMT_1", "celltype"))
  pdf(file.path(OUT, "02.EMT_score_celltype.pdf"), width = 7, height = 4)
  print(ggplot(df, aes(celltype, EMT_1, fill = celltype)) + geom_boxplot(outlier.size = 0.3) +
          theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "EMT module score by cell type")); dev.off()
}

## ── 4) 기능 농축: GSVA ssGSEA (논문: 186 KEGG, GSVA v2.2.0) — PCT vs dPCT ────
#   논문 Methods: "186개 KEGG 유전자세트, GSVA로 ssGSEA". 단일세포 전량은 무거워
#   PCT·dPCT 유사벌크(pseudobulk)로 경량 실행.
if (requireNamespace("GSVA", quietly = TRUE)) {
  suppressMessages(library(GSVA))
  kegg <- file.path(ROOT, "data/4-2. c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt")
  gs <- lapply(strsplit(readLines(kegg), "\t"), function(x) x[-c(1,2)])
  names(gs) <- sapply(strsplit(readLines(kegg), "\t"), `[`, 1)
  # PCT·dPCT 세포평균 발현(pseudobulk)
  expr <- as.matrix(GetAssayData(pdpct, layer = "data"))
  pb <- sapply(c("PCT","dPCT"), function(ct) rowMeans(expr[, pdpct$celltype==ct, drop=FALSE]))
  es <- tryCatch(gsva(gsvaParam(pb, gs)), error = function(e) gsva(pb, gs, method="ssgsea"))
  write.csv(es, file.path(OUT, "02.ssGSEA_KEGG_PCT_dPCT.csv"))
  cat("[GSVA] ssGSEA KEGG 저장 (PCT vs dPCT)\n")
} else cat("[GSVA] 패키지 없음 — 건너뜀 (BiocManager::install('GSVA'))\n")

saveRDS(Late, file.path(OUT, "Late_with_geneset_scores.RDS"))
cat("\n★ [02] 완료 → Fig7 B/C/E/G + 모듈스코어 + GSVA. 다음: 03_scrna_pseudotime.R\n")
