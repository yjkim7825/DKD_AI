# config.R ----------------------------------------------------------------
# 경로/파라미터 중앙 설정. 원본 스크립트의 setwd("G:\\187geneMR\\...") 하드코딩 대체.
# 로컬 환경에 맞게 이 파일만 수정하면 전체 R 파이프라인이 동작한다.
#
# 사용법(각 단계 스크립트 맨 위에서):
#   source("config.R")            # R 작업 디렉터리를 pipeline_port/R 로 두고 실행하거나
#   BASE_DIR <- "C:/.../pipeline_port"; source(file.path(BASE_DIR,"R","config.R"))

# ---- 기본 디렉터리 ----
# BASE_DIR 를 미리 지정하지 않았으면 현재 작업 디렉터리의 상위를 사용
if (!exists("BASE_DIR")) {
  BASE_DIR <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
}
DATA_DIR   <- file.path(BASE_DIR, "data")      # GEO 다운로드/중간 파일 위치
RESULT_DIR <- file.path(BASE_DIR, "results")   # 산출물
dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 그룹 라벨 규약 (열 이름 = '{샘플ID}_{그룹}') ----
CONTROL_LABEL <- "Control"
CASE_LABEL    <- "DKD"
EARLY_LABEL   <- "Early"
LATE_LABEL    <- "Late"

# ---- DEG 필터 (differential expression analysis.R 동일) ----
LOGFC_FILTER <- 0.585   # |log2FC| (=1.5배)
ADJP_FILTER  <- 0.05    # BH adj.P

# ---- 머신러닝 (machine learning modeling 1.R 동일) ----
RANDOM_SEED   <- 123
CV_FOLDS      <- 10
SVM_RFE_SIZES <- c(2, 3, 4, 5, 6, 7, 8)

message("[config] BASE_DIR = ", BASE_DIR)
