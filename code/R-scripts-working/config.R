# config.R — 파이프라인 중앙 설정 (경로는 여기 한 곳에서만 정의) ------------------
# 기존 로컬 설정 파일 계승 + here 기반 루트 자동탐지.
# 루트 자동탐지: here 패키지가 프로젝트 루트(.git / .here 앵커)를 찾음 → 스크립트가
#   하위 STEP 폴더에 있어도 항상 같은 루트를 참조(../data 깊이 문제 없음).
# ▶ 사용자는 보통 수정 불필요. 데이터 위치가 다르면 아래 DATA_ROOT 한 줄만 바꾸세요.
# --------------------------------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")
suppressMessages(library(here))

# ---- 루트 (자동탐지; 필요 시 여기만 수정) ----
CODE_ROOT <- here::here()                                              # = R-scripts-working
DATA_ROOT <- normalizePath(file.path(CODE_ROOT, "..", "data"), mustWork = FALSE)  # = code/data

# ---- 산출물/작업 폴더 (원본과 분리) ----
OUT_DIR     <- file.path(DATA_ROOT, "processed")   # 중간 매트릭스(전처리 결과)
RES_DIR     <- file.path(CODE_ROOT, "results")     # 분석 결과(DEG/ROC 등)
SCRATCH_DIR <- file.path(tempdir(), "pipeline_scratch")  # 임시 작업(비출력, 하드코딩 없음)
dir.create(OUT_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(SCRATCH_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 번호 붙은 원본 RAW 폴더 (읽기 전용) ----
DIR_GSE96804  <- file.path(DATA_ROOT, "1-1. GSE96804_RAW")
DIR_GSE104948 <- file.path(DATA_ROOT, "1-2. GSE104948_RAW")
DIR_GSE104954 <- file.path(DATA_ROOT, "1-3. GSE104954_RAW")
DIR_GSE142025 <- file.path(DATA_ROOT, "1-4. GSE142025_RAW")
DIR_GSE30529  <- file.path(DATA_ROOT, "1-5. GSE30529_RAW")

# ---- 그룹 라벨 규약 ----
GROUP_CONTROL <- "Control"
GROUP_EARLY   <- "Early"
GROUP_LATE    <- "Late"
GROUP_DKD     <- "DKD"      # early+late 통합 2군 비교 시

# ---- DEG 필터 (논문 동일) ----
LOGFC_FILTER <- 0.585
ADJP_FILTER  <- 0.05

message("[config] CODE_ROOT = ", CODE_ROOT)
message("[config] DATA_ROOT = ", DATA_ROOT)
message("[config] 산출물 -> ", OUT_DIR, " , ", RES_DIR)
