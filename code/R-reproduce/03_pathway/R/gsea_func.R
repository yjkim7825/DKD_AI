# ============================================================================
# gsea_func.R — GSEA 함수 (저자 GSEA.R 로직 그대로)
#   입력: DEG diff 결과(diff_data) + 경로 유전자셋(gmt: TERM2GENE)
#   과정: logFC로 유전자 순위 → GSEA → p.adjust<0.05 경로 필터
#   ※ 함수 정의만 — 실행은 run_gsea.R 에서
# ============================================================================
suppressMessages({ library(clusterProfiler); library(enrichplot); library(ggplot2) })

ADJP <- 0.05   # 경로 유의성 기준 (보정 p)

# diff_data = read.table(diff_*.txt, header=T, row.names=1)  → logFC 컬럼 + 유전자 행이름
# gmt       = read.gmt(...)  (term, gene 2열)
# tag       = 비교이름 (예: "Late_vs_Control")
# db        = "Hallmark" 또는 "KEGG" (출력 파일명 구분)
run_gsea <- function(diff_data, gmt, tag, db, outDir, min_terms = 5) {   # 논문: top 5 경로 시각화
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  figDir <- file.path(outDir, "figures")                  # 이미지 전용 폴더
  dir.create(figDir, showWarnings = FALSE, recursive = TRUE)

  # ── ① 유전자 순위 만들기 (logFC 내림차순) ──────────────────────────────
  logFC <- diff_data$logFC
  names(logFC) <- rownames(diff_data)
  logFC <- sort(logFC, decreasing = TRUE)   # 큰 양수(↑말기) → 큰 음수(↓말기)
  cat("\n== GSEA:", db, "-", tag, "==  입력 유전자", length(logFC), "개\n")

  # ── ② GSEA 실행 ───────────────────────────────────────────────────────
  set.seed(123)
  kk <- GSEA(logFC, TERM2GENE = gmt, pvalueCutoff = 1)   # 다 뽑고 뒤에서 필터
  kkTab <- as.data.frame(kk)

  # 전체 결과 저장
  write.table(kkTab, file = file.path(outDir, paste0("GSEA.", db, ".", tag, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # ── ③ 유의 경로 필터 (p.adjust < 0.05) ────────────────────────────────
  sig <- kkTab[kkTab$p.adjust < ADJP, ]
  up   <- sig[sig$NES > 0, ]   # 말기에 켜진(증가) 경로
  down <- sig[sig$NES < 0, ]   # 말기에 꺼진(감소) 경로
  cat("   유의 경로:", nrow(sig), " (↑켜짐", nrow(up), " / ↓꺼짐", nrow(down), ")\n")
  if (nrow(sig) > 0) {
    cat("   ↑ 상위:", paste(head(up$ID[order(-up$NES)], 3), collapse=", "), "\n")
    cat("   ↓ 하위:", paste(head(down$ID[order(down$NES)], 3), collapse=", "), "\n")
  }

  # ── ④ 시각화 (gseaplot2 — 상위 NES>0 / NES<0 각 min_terms개) ──────────
  plot_terms <- function(kk, terms, title, suffix) {
    if (length(terms) == 0) return(invisible())
    g <- gseaplot2(kk, geneSetID = terms, base_size = 10, title = title, color = "firebrick")
    ggsave(file.path(figDir, paste0("GSEA.", db, ".", tag, ".", suffix, ".pdf")),
           plot = g, width = 8, height = 6)
  }
  if (nrow(up)   > 0) plot_terms(kk, head(rownames(up[order(-up$NES),]),   min_terms), paste(tag,"↑켜진경로"), "Up")
  if (nrow(down) > 0) plot_terms(kk, head(rownames(down[order(down$NES),]), min_terms), paste(tag,"↓꺼진경로"), "Down")

  data.frame(비교 = tag, DB = db, 유의경로 = nrow(sig), 켜짐 = nrow(up), 꺼짐 = nrow(down))
}
