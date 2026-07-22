# ============================================================================
# run_gsea_annotation_compare.R — GSEA 두 방식 비교
#   A안(현재/저자 GSEA.R): 심볼 gmt 직접 (TERM2GENE=심볼, org.Hs.eg.db 안 씀)
#   B안(논문 Methods 명시): org.Hs.eg.db 로 심볼→ENTREZ 변환 후 GSEA
#   → 두 결과(NES·p.adjust)가 같은지 경로별로 대조
#   입력: DEG diff (Late_vs_Control) + Hallmark gmt
#   출력: 03_pathway/output/annotation_compare.txt
# ============================================================================
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
})

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
DEG  <- file.path(ROOT, "R-reproduce/02_deg/output/GSE142025_3group/diff_Late_vs_Control.txt")
GMT  <- file.path(ROOT, "data/4-1. h.all.v2026.1.Hs.symbols.gmt")
OUT  <- file.path(ROOT, "R-reproduce/03_pathway/output")

# DEG diff → logFC (심볼)
d <- read.table(DEG, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
gmt_sym <- read.gmt(GMT)   # term, gene(심볼)

## ── A안: 심볼 직접 (현재 방식) ──────────────────────────────────────────────
lfc_sym <- sort(setNames(d$logFC, rownames(d)), decreasing = TRUE)
set.seed(123)
gseaA <- as.data.frame(GSEA(lfc_sym, TERM2GENE = gmt_sym, pvalueCutoff = 1))
cat("A안(심볼): 유의경로", sum(gseaA$p.adjust < 0.05), "개\n")

## ── B안: org.Hs.eg.db 로 ENTREZ 변환 후 ─────────────────────────────────────
# 1) DEG 심볼 → ENTREZ
map <- bitr(rownames(d), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
d2 <- merge(data.frame(SYMBOL = rownames(d), logFC = d$logFC), map, by = "SYMBOL")
lfc_ent <- sort(setNames(d2$logFC, d2$ENTREZID), decreasing = TRUE)
cat("심볼→ENTREZ 매핑:", nrow(d), "→", nrow(d2), " (매핑 실패", nrow(d)-nrow(d2), "개 탈락)\n")
# 2) gmt 심볼 → ENTREZ (같은 org.Hs.eg.db)
gmap <- bitr(unique(gmt_sym$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gmt_ent <- merge(gmt_sym, gmap, by.x = "gene", by.y = "SYMBOL")[, c("term","ENTREZID")]
colnames(gmt_ent) <- c("term","gene")
# 3) GSEA (ENTREZ)
set.seed(123)
gseaB <- as.data.frame(GSEA(lfc_ent, TERM2GENE = gmt_ent, pvalueCutoff = 1))
cat("B안(ENTREZ): 유의경로", sum(gseaB$p.adjust < 0.05), "개\n")

## ── 두 방식 경로별 대조 ─────────────────────────────────────────────────────
cmp <- merge(
  data.frame(ID = gseaA$ID, NES_A = round(gseaA$NES,3), padj_A = signif(gseaA$p.adjust,3), size_A = gseaA$setSize),
  data.frame(ID = gseaB$ID, NES_B = round(gseaB$NES,3), padj_B = signif(gseaB$p.adjust,3), size_B = gseaB$setSize),
  by = "ID", all = TRUE)
cmp$dNES <- round(abs(cmp$NES_A - cmp$NES_B), 3)   # NES 차이 절대값
cmp <- cmp[order(-cmp$dNES), ]
write.table(cmp, file.path(OUT, "annotation_compare.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n===== 두 방식 대조 (NES 차이 큰 순 상위) =====\n")
print(head(cmp, 10), row.names = FALSE)
cat("\nNES 평균 차이:", round(mean(cmp$dNES, na.rm=TRUE), 4),
    " / 최대 차이:", round(max(cmp$dNES, na.rm=TRUE), 3), "\n")
cat("유의(p<0.05) 경로 수 — A:", sum(gseaA$p.adjust<0.05), " B:", sum(gseaB$p.adjust<0.05), "\n")
cat("\n★ 완료 — output/annotation_compare.txt\n")
