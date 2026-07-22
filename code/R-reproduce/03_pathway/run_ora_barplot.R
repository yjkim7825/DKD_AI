# ============================================================================
# run_ora_barplot.R — Fig4 C/D/E/G 재현 (ORA 막대그래프)
#   논문 C/D/E/G = ORA(과대표현분석): DEG 목록이 어떤 Hallmark 경로에 몰렸나
#     x축 = Count(경로에 걸린 DEG 수), 색 = p.adjust
#   ※ 저자 GSEA.R엔 이 코드 없음(공유 안 함) → 논문 그림 보고 재현
#   입력: 02_deg/output/.../diff_*.txt (유의 DEG 목록) + Hallmark gmt
#   출력: 03_pathway/output/figures/ORA.Hallmark.<비교>.pdf
# ============================================================================
suppressMessages({ library(clusterProfiler); library(ggplot2); library(enrichplot) })

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
PATH <- file.path(ROOT, "R-reproduce/03_pathway")
DEG  <- file.path(ROOT, "R-reproduce/02_deg/output")
GMT_H <- file.path(ROOT, "data/4-1. h.all.v2026.1.Hs.symbols.gmt")
OUTF  <- file.path(PATH, "output/figures")            # ★ 이미지 전용 폴더
dir.create(OUTF, showWarnings = FALSE, recursive = TRUE)

gmtH <- read.gmt(GMT_H)

# diff 목록 읽기 → 유의 DEG 유전자 벡터 (up+down 전부)
read_genes <- function(path) {
  d <- read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  rownames(d)
}

# ORA 실행 + 막대그래프 (Count, p.adjust 색) — 논문 C/D/E/G 스타일
ora_bar <- function(genes, tag, outDir, showN = 15) {
  ora <- enricher(genes, TERM2GENE = gmtH, pvalueCutoff = 1, qvalueCutoff = 1)
  tab <- as.data.frame(ora)
  sig <- tab[tab$p.adjust < 0.05, ]
  cat("\n== ORA:", tag, "==  유의 경로", nrow(sig), "개\n")
  if (nrow(sig) == 0) { cat("   (유의 경로 없음)\n"); return(NULL) }

  write.table(tab, file.path(outDir, paste0("ORA.Hallmark.", tag, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  d <- sig[order(sig$Count, decreasing = TRUE), ]
  d <- head(d, showN)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  p <- ggplot(d, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_col() +
    scale_fill_gradient(low = "#d73027", high = "#4575b4") +   # 논문: 빨강(유의)→파랑
    labs(title = tag, x = "Count", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(outDir, paste0("ORA.Hallmark.", tag, ".pdf")),
         plot = p, width = 8, height = 5)
  cat("   상위:", paste(head(d$Description, 3), collapse = ", "), "\n")
  data.frame(비교 = tag, 유의경로 = nrow(sig))
}

# ── 실행: 142025 3비교 (C/D/E) + GSE30529 검증 (G) ──────────────────────────
jobs <- list(
  Early_vs_Control = file.path(DEG, "GSE142025_3group/diff_Early_vs_Control.txt"),  # C
  Late_vs_Control  = file.path(DEG, "GSE142025_3group/diff_Late_vs_Control.txt"),   # D
  Late_vs_Early    = file.path(DEG, "GSE142025_3group/diff_Late_vs_Early.txt"),     # E
  GSE30529_DKD_vs_Control = file.path(DEG, "GSE30529_2group/diff_GSE30529_DKD_vs_Control.txt")  # G(검증)
)
res <- list()
for (tag in names(jobs)) res[[tag]] <- ora_bar(read_genes(jobs[[tag]]), tag, OUTF)

summ <- do.call(rbind, res); rownames(summ) <- NULL
cat("\n===== ORA 유의 경로 요약 =====\n"); print(summ)
write.table(summ, file.path(OUTF, "ORA_summary.txt"), sep = "\t", quote = F, row.names = F)
cat("\n★ ORA 막대그래프 완료 — output/figures/ 확인\n")
