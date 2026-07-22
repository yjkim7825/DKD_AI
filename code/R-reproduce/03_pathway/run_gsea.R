# ============================================================================
# run_gsea.R — GSEA 실행 (GSE142025 3비교 × Hallmark/KEGG)
#   입력: 02_deg/output/GSE142025_3group/diff_*.txt  (경로 참조, 복사 X)
#   유전자셋: data/4-1 Hallmark(50) , 4-2 KEGG_legacy(186)
#   출력: 03_pathway/output/
#   기준: p.adjust < 0.05 (경로 유의)
# ============================================================================
for (p in c("clusterProfiler","enrichplot","ggplot2"))
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
PATH <- file.path(ROOT, "R-reproduce/03_pathway")
DEG  <- file.path(ROOT, "R-reproduce/02_deg/output/GSE142025_3group")   # DEG diff 읽는 곳
GMT_H <- file.path(ROOT, "data/4-1. h.all.v2026.1.Hs.symbols.gmt")        # Hallmark
GMT_K <- file.path(ROOT, "data/4-2. c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt")  # KEGG
OUT   <- file.path(PATH, "output")
source(file.path(PATH, "R/gsea_func.R"))

# 경로 유전자셋 로드 (term, gene 2열)
gmtH <- read.gmt(GMT_H);  cat("Hallmark 경로:", length(unique(gmtH$term)), "\n")
gmtK <- read.gmt(GMT_K);  cat("KEGG 경로:",     length(unique(gmtK$term)), "\n")

# DEG diff 파일 읽기 (저자와 동일: header=T, row.names=1 → logFC 컬럼 + 유전자 행이름)
read_diff <- function(tag)
  read.table(file.path(DEG, paste0("diff_", tag, ".txt")),
             header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

comparisons <- c("Early_vs_Control", "Late_vs_Control", "Late_vs_Early")

summary_all <- list()
for (tag in comparisons) {
  d <- read_diff(tag)
  summary_all[[paste0(tag,"_H")]] <- run_gsea(d, gmtH, tag, "Hallmark", OUT)
  summary_all[[paste0(tag,"_K")]] <- run_gsea(d, gmtK, tag, "KEGG",     OUT)
}

## ── 검증셋 GSE30529 (DKD vs Control) — 논문: "validated in GSE30529" ─────────
DEG30529 <- file.path(ROOT, "R-reproduce/02_deg/output/GSE30529_2group")
d30 <- read.table(file.path(DEG30529, "diff_GSE30529_DKD_vs_Control.txt"),
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
summary_all[["GSE30529_H"]] <- run_gsea(d30, gmtH, "GSE30529_DKD_vs_Control", "Hallmark", OUT)
summary_all[["GSE30529_K"]] <- run_gsea(d30, gmtK, "GSE30529_DKD_vs_Control", "KEGG",     OUT)

# ── 요약 ────────────────────────────────────────────────────────────────────
res <- do.call(rbind, summary_all); rownames(res) <- NULL
cat("\n===== GSEA 유의 경로 요약 =====\n"); print(res)
write.table(res, file.path(OUT, "GSEA_summary.txt"), sep="\t", quote=F, row.names=F)
cat("\n★ GSEA 완료 — output/ 확인\n")
# ============================================================================
# 출력(비교×DB 마다):
#   GSEA.<DB>.<비교>.txt        전체 경로 결과(NES, p.adjust 등)
#   GSEA.<DB>.<비교>.Up.pdf     ↑켜진 상위경로 gseaplot2
#   GSEA.<DB>.<비교>.Down.pdf   ↓꺼진 상위경로
#   GSEA_summary.txt            비교×DB 유의경로 개수 요약
# ============================================================================
