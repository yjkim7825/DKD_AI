# ============================================================================
# 04_scrna_cellchat.R — scRNA 4/4: CellChat 세포간 통신 (Control vs Late_DKD 비교)
#   저자 "single-cell RNA-seq analysis 4.R" 재현
#   입력 : output/Late_annotated.RDS (또는 with_geneset_scores) — group으로 분할
#   출력 : 그룹별 통신 net CSV, 비교 그림(수·강도·rankNet·heatmap·bubble)
#   ※ CellChat 필요. 무거움(computeCommunProb).
# ============================================================================
set.seed(123)
# 논문 버전: CellChat v2.2.0 (Secreted Signaling, 리간드-수용체 통신)
suppressMessages({
  library(CellChat); library(Seurat); library(patchwork); library(dplyr)
})

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
OUT  <- file.path(ROOT, "R-reproduce/06_scrna/output")
Late <- readRDS(file.path(OUT, "Late_annotated.RDS"))
Idents(Late) <- Late$celltype

## ── 그룹별 CellChat 객체 만들기 ─────────────────────────────────────────────
make_cc <- function(obj) {
  data.input <- GetAssayData(obj, assay = "RNA", layer = "data")
  meta <- data.frame(celltype = obj$celltype, row.names = colnames(obj))
  cc <- createCellChat(object = as.matrix(data.input), meta = meta, group.by = "celltype")
  cc <- setIdent(cc, ident.use = "celltype")
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")   # 분비신호
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)   # presto 없이 표준 Wilcoxon
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, raw.use = TRUE, type = "truncatedMean", trim = 0.1)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

ctrl <- make_cc(subset(Late, subset = group == "Control"))
late <- make_cc(subset(Late, subset = group == "Late_DKD"))
write.csv(subsetCommunication(ctrl), file.path(OUT, "04.net_Control.csv"), row.names = FALSE)
write.csv(subsetCommunication(late), file.path(OUT, "04.net_LateDKD.csv"), row.names = FALSE)
saveRDS(ctrl, file.path(OUT, "cellchat_control.RDS"))
saveRDS(late, file.path(OUT, "cellchat_late.RDS"))

## ── 두 그룹 병합 후 비교 ────────────────────────────────────────────────────
object.list <- list(Control = ctrl, LateDKD = late)
cc <- mergeCellChat(object.list, add.names = names(object.list))

pdf(file.path(OUT, "04.compare_interactions.pdf"), width = 8, height = 4)
print(compareInteractions(cc, group = c(1, 2), measure = "count") +
      compareInteractions(cc, group = c(1, 2), measure = "weight")); dev.off()

pdf(file.path(OUT, "04.rankNet.pdf"), width = 7, height = 6)
print(rankNet(cc, mode = "comparison", stacked = TRUE, do.stat = TRUE)); dev.off()

# dPCT(source)가 Late_DKD에서 증강시키는 신호 (저자 강조)
pdf(file.path(OUT, "04.dPCT_signaling_LateDKD.pdf"), width = 8, height = 7)
print(netVisual_bubble(cc, sources.use = "dPCT", targets.use = levels(Late$celltype),
                       comparison = c(1, 2), max.dataset = 2, thresh = 0.001,
                       title.name = "dPCT-increased signaling in Late DKD")); dev.off()

cat("\n★ [04] 완료 → CellChat 비교 (Control vs Late DKD). dPCT 신호 증강 확인.\n")
