# ============================================================================
# run_mr.R — MR 전체 한 번에 (저자처럼 1→2→3 순서로 실행)
#   개별로 돌리려면 01_mr_run.R → 02_mr_filter.R → 03_mr_forest.R 각각 source
# ============================================================================
MRDIR <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/R-reproduce/04_mr"
source(file.path(MRDIR, "01_mr_run.R"))     # 노출→결과→MR→CSV
source(file.path(MRDIR, "02_mr_filter.R"))  # IVW 필터
source(file.path(MRDIR, "03_mr_forest.R"))  # forest plot
cat("\n===== MR 1·2·3 전체 완료 =====\n")
