# ============================================================================
# run_deg_3group_plots.R — 정상·초기·말기 3그룹 통합 시각화
#   ① 히트맵: Control/Early/Late 세 그룹 컬러바로 한 번에
#   ② 볼케이노: 말기vs정상 DEG를
#        - 초기부터 이미 변함(초기vs정상에도 있음) = 주황
#        - 말기에만 새로 변함                      = 빨강/파랑
#   입력: 전처리 발현행렬 + 02_deg 의 diff 결과
#   출력: 02_deg/output/GSE142025_3group/  (heatmap_3group.pdf, vol_3group.pdf)
# ============================================================================
suppressMessages({ library(limma); library(pheatmap); library(ggplot2); library(dplyr) })

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
IN   <- file.path(ROOT, "R-reproduce/01_preprocessing/bulk/output/02_normalized/GSE142025.normalized.txt")
DEGDIR <- file.path(ROOT, "R-reproduce/02_deg/output/GSE142025_3group")            # DEG 결과 읽는 곳
OUT    <- file.path(ROOT, "R-reproduce/02_deg/output/GSE142025_3group_combined")   # 통합 그림 저장 곳
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
LOGFC <- 0.585; ADJP <- 0.05

# 발현행렬 + 그룹
rt  <- read.table(IN, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
mat <- as.matrix(rt)
grp <- sapply(strsplit(colnames(mat), "_"), function(x) tail(x,1))
# 그룹 순서: Control → Early → Late (열 정렬)
ord <- order(factor(grp, levels=c("Control","Early","Late")))
mat <- mat[, ord]; grp <- grp[ord]

# DEG 결과 읽기
diffL <- read.table(file.path(DEGDIR,"diff_Late_vs_Control.txt"),  header=T, sep="\t", row.names=1, check.names=F)
diffE <- read.table(file.path(DEGDIR,"diff_Early_vs_Control.txt"), header=T, sep="\t", row.names=1, check.names=F)
allL  <- read.table(file.path(DEGDIR,"all_Late_vs_Control.txt"),   header=T, sep="\t", row.names=1, check.names=F)

## ── ① 3그룹 히트맵 ──────────────────────────────────────────────────────────
# 말기vs정상 DEG 중 top50↑ + top50↓ (유의성 순)
sig <- diffL[order(diffL$adj.P.Val), ]
up  <- head(rownames(sig[sig$logFC>0,]), 50)
dn  <- head(rownames(sig[sig$logFC<0,]), 50)
genes <- c(up, dn); genes <- genes[genes %in% rownames(mat)]
hm <- mat[genes, , drop=FALSE]

ann <- data.frame(Group = factor(grp, levels=c("Control","Early","Late")))
rownames(ann) <- colnames(mat)
ann_colors <- list(Group = c(Control="#2ca02c", Early="#ff7f0e", Late="#d62728"))  # 초록/주황/빨강

pdf(file.path(OUT, "heatmap_3group.pdf"), width=11, height=8)
pheatmap(hm,
         annotation_col = ann, annotation_colors = ann_colors,
         color = colorRampPalette(c("blue2","white","red2"))(50),
         cluster_cols = FALSE, show_colnames = FALSE, scale = "row",
         gaps_col = cumsum(table(factor(grp, levels=c("Control","Early","Late")))),  # 그룹 사이 간격
         main = "GSE142025 top DEG — Control / Early / Late",
         fontsize = 9, fontsize_row = 5.5)
dev.off()
cat("① 히트맵 저장: heatmap_3group.pdf  (top", length(genes), "유전자)\n")

## ── ② 볼케이노: 3구간(정상↔초기↔말기)으로 3분류 색 구분 ───────────────────
#   N→E = 초기vs정상(diffE), E→L = 말기vs초기(diffLE)
#   ① 초기에만(N→E만)  ② 둘다(N→E∩E→L)  ③ 말기에만(E→L만)
diffLE <- read.table(file.path(DEGDIR,"diff_Late_vs_Early.txt"), header=T, sep="\t", row.names=1, check.names=F)
NE <- rownames(diffE); EL <- rownames(diffLE)
only_early <- setdiff(NE, EL); both <- intersect(NE, EL); only_late <- setdiff(EL, NE)

# x축 = 최종 변화(말기vs정상 logFC) 로 놓고 3분류 색칠
v <- allL; v$gene <- rownames(v); v$cat <- "Not (유의X)"
v$cat[v$gene %in% only_early] <- "① 초기에만 (N→E)"
v$cat[v$gene %in% both]       <- "② 둘다 (N→E & E→L)"
v$cat[v$gene %in% only_late]  <- "③ 말기에만 (E→L)"

cols <- c("Not (유의X)"="grey85",
          "① 초기에만 (N→E)"  ="#ff7f0e",   # 주황 = 초기에만
          "② 둘다 (N→E & E→L)"="#2ca02c",   # 초록 = 계속 변함
          "③ 말기에만 (E→L)"  ="#d62728")   # 빨강 = 말기에만

p <- ggplot(v, aes(logFC, -log10(adj.P.Val), color=cat)) +
  geom_point(size=0.9, alpha=0.75) +
  scale_color_manual(values=cols, name="변화 구간") +
  geom_vline(xintercept=c(-LOGFC,LOGFC), linetype="dashed", color="grey50") +
  geom_hline(yintercept=-log10(ADJP), linetype="dashed", color="grey50") +
  labs(title="변화 구간별 DEG (x=말기vs정상 logFC)",
       subtitle="초기에만(주황) / 둘다-계속(초록) / 말기에만(빨강)") +
  theme_bw() + theme(plot.title=element_text(size=12,hjust=.5,face="bold"))
pdf(file.path(OUT, "vol_3stage.pdf"), width=8, height=6); print(p); dev.off()

# 벤 다이어그램 (N→E vs E→L 겹침) — 3분류를 벤으로도
suppressMessages(if(!requireNamespace("VennDiagram",quietly=TRUE)) install.packages("VennDiagram"))
if (requireNamespace("VennDiagram", quietly=TRUE)) {
  library(VennDiagram)
  v2 <- venn.diagram(list("정상→초기 (N→E)"=NE, "초기→말기 (E→L)"=EL),
                     filename=NULL, fill=c("#ff7f0e","#d62728"), alpha=0.5,
                     cex=1.5, cat.cex=1.1, main="변화 구간 겹침")
  pdf(file.path(OUT,"venn_stage.pdf"), width=6, height=6); grid::grid.draw(v2); dev.off()
}

cat("② 3분류 볼케이노 저장: vol_3stage.pdf + venn_stage.pdf\n")
cat("   ① 초기에만:", length(only_early),
    " / ② 둘다(계속):", length(both),
    " / ③ 말기에만:", length(only_late), "\n")

## ── ③ 가장 많이 변한 유전자 Top 테이블 (3분류별) ──────────────────────────
# 각 유전자에 말기vs정상 logFC 붙여서 |logFC| 큰 순
topTab <- function(genes, label, useDiff, n=15) {
  gg <- intersect(genes, rownames(useDiff))
  if (length(gg)==0) return(NULL)
  t <- useDiff[gg, c("logFC","adj.P.Val")]
  t <- t[order(-abs(t$logFC)), ]; t <- head(t, n)
  data.frame(구간=label, gene=rownames(t), logFC=round(t$logFC,2),
             adjP=signif(t$adj.P.Val,2), 방향=ifelse(t$logFC>0,"↑","↓"))
}
top_all <- rbind(
  topTab(only_early, "①초기에만", diffE),    # 초기 기준 logFC
  topTab(both,       "②둘다",     diffL),    # 말기vs정상 기준
  topTab(only_late,  "③말기에만", diffLE)    # 말기vs초기 기준
)
write.table(top_all, file.path(OUT,"TOP_3stage.txt"), sep="\t", quote=F, row.names=F)
cat("③ 3분류 Top 테이블 저장: TOP_3stage.txt\n")
cat("\n=== 각 구간 최다 변화 유전자 ===\n")
for (lab in c("①초기에만","②둘다","③말기에만")) {
  s <- top_all[top_all$구간==lab, ][1, ]
  if (!is.na(s$gene)) cat(sprintf("  %s: %s (logFC %.2f %s)\n", lab, s$gene, s$logFC, s$방향))
}
cat("\n★ 3그룹 시각화 완료 — output/GSE142025_3group_combined/ 확인\n")
