# ============================================================================
# run_deg.R — DEG 실행 (GSE142025 3비교 + GSE30529 2비교 검증)
#   입력: 전처리 결과 경로 직접 참조 (복사 X)
#     01_preprocessing/bulk/output/02_normalized/*.normalized.txt
#   출력: 02_deg/output/GSE142025_3group/ , GSE30529_2group/
#   기준: |logFC|>0.585 & adj.P<0.05  (논문 동일)
# ============================================================================
if (!requireNamespace("pheatmap", quietly=TRUE)) install.packages("pheatmap")

ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
DEG  <- file.path(ROOT, "R-reproduce/02_deg")
IN   <- file.path(ROOT, "R-reproduce/01_preprocessing/bulk/output/02_normalized")  # 입력 경로 참조
source(file.path(DEG, "R/deg_limma.R"))

summary_all <- list()

## ── 1) GSE142025 (3그룹) — 단계별 확인 : 3비교 ──────────────────────────────
cat("\n########## GSE142025 (3그룹: Control/Early/Late) ##########\n")
e <- read_expr(file.path(IN, "GSE142025.normalized.txt"))
out142 <- file.path(DEG, "output/GSE142025_3group")
summary_all[["Late_vs_Early"]]   <- deg_two(e$mat, e$grp, "Late",  "Early",   out142, "Late_vs_Early")
summary_all[["Late_vs_Control"]] <- deg_two(e$mat, e$grp, "Late",  "Control", out142, "Late_vs_Control")
summary_all[["Early_vs_Control"]]<- deg_two(e$mat, e$grp, "Early", "Control", out142, "Early_vs_Control")

## ── 2) GSE30529 (2그룹) — 검증 : DKD vs Control ────────────────────────────
cat("\n########## GSE30529 (2그룹: DKD/Control) — 검증 ##########\n")
e2 <- read_expr(file.path(IN, "GSE30529.normalized.txt"))
out305 <- file.path(DEG, "output/GSE30529_2group")
summary_all[["GSE30529_DKD_vs_Control"]] <- deg_two(e2$mat, e2$grp, "DKD", "Control", out305, "GSE30529_DKD_vs_Control")

## ── 요약 + 논문 대조 ────────────────────────────────────────────────────────
res <- do.call(rbind, summary_all); rownames(res) <- NULL
paper <- c(Late_vs_Early=2833, Late_vs_Control=3525, Early_vs_Control=390)
res$논문 <- paper[res$비교]
cat("\n===== DEG 요약 (우리 vs 논문) =====\n")
print(res)
write.table(res, file.path(DEG, "output/DEG_summary.txt"), sep="\t", quote=F, row.names=F)
cat("\n★ DEG 완료 — output/ 하위 폴더 확인\n")
# ============================================================================
# 출력 파일(각 비교마다):
#   all_*.txt      전체 유전자 통계
#   diff_*.txt     유의 DEG (|logFC|>0.585 & adj.P<0.05)
#   up_*.txt/down_*.txt   상향/하향 유전자 목록
#   heatmap_*.pdf  히트맵 (top50↑+50↓)
#   vol_*.pdf      볼케이노
# ============================================================================
