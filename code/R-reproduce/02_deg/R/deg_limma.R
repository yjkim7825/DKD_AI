# ============================================================================
# deg_limma.R — DEG 분석 함수 (두 그룹 비교, 2·3그룹 공용)  [저자 limma 방식]
#   저자 differential expression analysis.R 로직 그대로:
#     design → lmFit → makeContrasts → eBayes → topTable
#   기준: |logFC| > 0.585 (=1.5배)  &  adj.P.Val < 0.05
#   ※ 함수 정의만 — 실행은 run_deg.R 에서
# ============================================================================
suppressMessages({ library(limma); library(dplyr); library(pheatmap); library(ggplot2) })

LOGFC <- 0.585      # 논문 기준 (1.5배)
ADJP  <- 0.05

# 발현행렬 읽기 (전처리 결과 *.normalized.txt) → matrix + 그룹벡터
read_expr <- function(path) {
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  mat <- as.matrix(rt)
  grp <- sapply(strsplit(colnames(mat), "_"), function(x) tail(x, 1))   # 열이름 뒤 라벨
  list(mat = mat, grp = grp)
}

# 두 그룹 DEG (caseGrp vs ctrlGrp) — caseGrp가 '비교 기준(양수=case에서 증가)'
#   예) deg_two(mat, grp, "Late", "Control")  → Late vs Control
deg_two <- function(mat, grp, caseGrp, ctrlGrp, outDir, tag) {
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  keep <- grp %in% c(caseGrp, ctrlGrp)
  d <- mat[, keep]; g <- grp[keep]
  d <- d[, order(g)]; g <- sort(g)                       # 그룹 정렬
  cat("\n== ", tag, " ==  ", paste(names(table(g)), table(g), sep="=", collapse=" "), "\n")

  # limma
  design <- model.matrix(~0 + factor(g)); colnames(design) <- levels(factor(g))
  fit <- lmFit(d, design)
  cont <- makeContrasts(contrasts = paste0(caseGrp, " - ", ctrlGrp), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  allDiff <- topTable(fit2, adjust = "fdr", number = 200000)

  # 유의 DEG
  sig <- allDiff[abs(allDiff$logFC) > LOGFC & allDiff$adj.P.Val < ADJP, ]
  up  <- rownames(sig[sig$logFC > 0, ]); dn <- rownames(sig[sig$logFC < 0, ])
  cat("   DEG:", nrow(sig), " (↑", length(up), " / ↓", length(dn), ")\n")

  # 저장: 전체 / 유의 / 상하향 목록
  write.table(rbind(id = colnames(allDiff), allDiff),
              file.path(outDir, paste0("all_", tag, ".txt")), sep="\t", quote=F, col.names=F)

  # ★ Top 테이블: 가장 많이 변한 유전자 (|logFC| 큰 순 Top20)
  if (nrow(sig) > 0) {
    topN <- sig[order(-abs(sig$logFC)), c("logFC","adj.P.Val")]
    topN <- head(topN, 20); topN$방향 <- ifelse(topN$logFC>0, "↑Up", "↓Down")
    topN <- cbind(gene = rownames(topN), round(topN[,1:2],3), 방향=topN$방향)
    write.table(topN, file.path(outDir, paste0("TOP20_", tag, ".txt")),
                sep="\t", quote=F, row.names=F)
    cat("   Top1 변화:", rownames(sig)[which.max(abs(sig$logFC))],
        "(logFC", round(sig$logFC[which.max(abs(sig$logFC))],2), ")\n")
  }
  if (nrow(sig) > 0) {
    write.table(rbind(id = colnames(sig), sig),
                file.path(outDir, paste0("diff_", tag, ".txt")), sep="\t", quote=F, col.names=F)
    writeLines(up, file.path(outDir, paste0("up_", tag, ".txt")))
    writeLines(dn, file.path(outDir, paste0("down_", tag, ".txt")))
  }

  # 그림 (히트맵 top50↑+50↓ / 볼케이노)
  if (nrow(sig) > 0) {
    gu <- head(rownames(sig[sig$logFC>0,]), 50); gd <- head(rownames(sig[sig$logFC<0,]), 50)
    hm <- d[c(gu, gd), , drop=FALSE]
    ann <- data.frame(Group = g); rownames(ann) <- colnames(d)
    pdf(file.path(outDir, paste0("heatmap_", tag, ".pdf")), width=10, height=7)
    pheatmap(hm, annotation_col=ann, color=colorRampPalette(c("blue2","white","red2"))(50),
             cluster_cols=FALSE, show_colnames=FALSE, scale="row",
             fontsize=8, fontsize_row=5.5)
    dev.off()
  }
  vd <- allDiff
  vd$Sig <- ifelse(vd$adj.P.Val < ADJP & abs(vd$logFC) > LOGFC,
                   ifelse(vd$logFC > 0, "Up", "Down"), "Not")
  p <- ggplot(vd, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = Sig)) +
    scale_color_manual(values = c("Down"="blue2","Not"="grey","Up"="red2")) +
    labs(title = tag) + theme(plot.title = element_text(size=14, hjust=.5, face="bold"))
  pdf(file.path(outDir, paste0("vol_", tag, ".pdf")), width=5.5, height=4.5); print(p); dev.off()

  data.frame(비교=tag, DEG=nrow(sig), Up=length(up), Down=length(dn))
}
