# ============================================================================
# run_preprocess.R  — 실행 코드 (여기서 데이터 지정 + prep1→2→3 순서대로)
#   흐름: 원본 → prep1(유전자행렬) → prep2(정규화·라벨) → prep3(검증셋 ComBat)
#   출력: output/01_genematrix , 02_normalized , 03_merged
#   ※ data/processed 는 안 건드림. 결과는 output/ 에만.
# ============================================================================

## ── 0) 패키지 확인/설치 ────────────────────────────────────────────────────
need <- c("affy","oligo","limma","GEOquery","sva","makecdfenv","R.utils",
          "hgu133a2.db","org.Hs.eg.db","hta20transcriptcluster.db","pd.hta.2.0")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(p, update = FALSE, ask = FALSE)
}

## ── 1) 경로 + 모듈 불러오기 ────────────────────────────────────────────────
ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
RUN  <- file.path(ROOT, "R-reproduce/01_preprocessing/bulk")
DATA <- file.path(ROOT, "data")
OUT1 <- file.path(RUN, "output/01_genematrix")
OUT2 <- file.path(RUN, "output/02_normalized")
OUT3 <- file.path(RUN, "output/03_merged")
for (d in c(OUT1, OUT2, OUT3)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

source(file.path(RUN, "R/prep1_make_matrix.R"))
source(file.path(RUN, "R/prep2_normalize_label.R"))
source(file.path(RUN, "R/prep3_merge_combat.R"))

save_txt <- function(mat, dir, name) {
  out <- cbind(geneNames = rownames(mat), as.data.frame(mat))
  write.table(out, file.path(dir, paste0(name, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("   저장:", name, " (", nrow(mat), "유전자 ×", ncol(mat), "샘플 )\n")
}

## ── 2) 데이터셋 정의 (★ 여기만 바꾸면 됨) ──────────────────────────────────
#   method: "affy"  = 구형 칩(U133A/A2) → ReadAffy + cdf + db
#           "oligo" = 신형 칩(HTA-2.0)  → read.celfiles + pd패키지 + db
#   method: "affy"/"oligo"/"brainarray"
arrays <- list(
  GSE30529  = list(dir = "1-5. GSE30529_RAW",  method = "affy",  cdf = "hgu133a2cdf", db = "hgu133a2.db"),
  GSE96804  = list(dir = "1-1. GSE96804_RAW",  method = "oligo", db = "hta20transcriptcluster.db"),
  GSE104948 = list(dir = "1-2. GSE104948_RAW", method = "brainarray",
                   chipmap = "GSE104948_chipmap.csv",
                   cdfmap = list("712"  = "GPL24120_HGU133A_Hs_ENTREZG.cdf.gz",
                                 "1164" = "GPL22945_HGU133Plus2_Hs_ENTREZG.cdf.gz")),
  GSE104954 = list(dir = "1-3. GSE104954_RAW", method = "brainarray",
                   chipmap = "GSE104954_chipmap.csv",
                   cdfmap = list("712"  = "GPL24120_HGU133A_Hs_ENTREZG.cdf.gz",
                                 "1164" = "GPL22945_HGU133Plus2_Hs_ENTREZG.cdf.gz"))
)
labeled <- list()   # 정규화·라벨 끝난 행렬 보관 (prep3에서 씀)

## ── 3) 마이크로어레이: prep1 → prep2 ──────────────────────────────────────
for (gse in names(arrays)) {
  a <- arrays[[gse]]; cat("\n===== ", gse, " (", a$method, ") =====\n")
  res <- tryCatch(
    if (a$method == "oligo")      prep1_array_oligo(file.path(DATA, a$dir), a$db)
    else if (a$method == "brainarray") prep1_array_brainarray(file.path(DATA, a$dir),
                                          lapply(a$cdfmap, function(x) file.path(DATA, a$dir, x)),
                                          file.path(RUN, a$chipmap))
    else                          prep1_array(file.path(DATA, a$dir), a$cdf, a$db),
    error = function(e) { message("  prep1 실패: ", conditionMessage(e)); NULL })
  if (is.null(res)) next
  # brainarray 는 list(mat, grp), 나머지는 matrix
  if (is.list(res)) { m <- res$mat; grp_pre <- res$grp } else { m <- res; grp_pre <- NULL }
  save_txt(m, OUT1, paste0(gse, ".genematrix"))
  # ★ 옵션1: 마이크로어레이는 RMA가 이미 log2+정규화 완료 → prep2(정규화) 스킵, 라벨만
  grp <- if (!is.null(grp_pre)) grp_pre else label_by_geo(gse, colnames(m), file.path(DATA,"processed"))
  m <- apply_labels(m, grp)
  save_txt(m, OUT2, paste0(gse, ".normalized"))   # (어레이는 사실상 RMA결과+라벨)
  labeled[[gse]] <- m
}

## ── 4) RNA-seq(GSE142025): prep1 → prep2 ──────────────────────────────────
cat("\n===== GSE142025 (RNA-seq) =====\n")
m <- tryCatch(prep1_rnaseq(file.path(DATA, "1-4. GSE142025_RAW")),
              error = function(e) { message("  prep1 실패: ", conditionMessage(e)); NULL })
if (!is.null(m)) {
  save_txt(m, OUT1, "GSE142025.genematrix")
  m <- prep2_normalize(m)                              # RNA-seq는 여기서 log2+정규화 실제로 수행
  m <- apply_labels(m, label_by_prefix(colnames(m)))   # A=Late/B=Early/N=Control
  save_txt(m, OUT2, "GSE142025.normalized")
  labeled[["GSE142025"]] <- m
}

## ── 5) prep3: 검증셋(104948 + 104954) 병합 + ComBat ───────────────────────
cat("\n===== ComBat: 검증셋 (GSE104948 + GSE104954) =====\n")
if (all(c("GSE104948","GSE104954") %in% names(labeled))) {
  merged <- prep3_merge_combat(list(labeled$GSE104948, labeled$GSE104954))
  save_txt(merged, OUT3, "valid_104948_104954.combat")
} else {
  cat("   (검증셋 둘 중 하나가 없어 ComBat 생략)\n")
}

cat("\n★ 전처리 완료 — output/ 하위 3폴더 확인\n")
# ============================================================================
# 폴더 구조:
#   R-reproduce/01_preprocessing/
#     R/  prep1_make_matrix.R / prep2_normalize_label.R / prep3_merge_combat.R  (기능)
#     run_preprocess.R                                                          (실행)
#     output/ 01_genematrix / 02_normalized / 03_merged                        (결과)
# 실행:  source("C:/.../R-reproduce/01_preprocessing/run_preprocess.R")
# ============================================================================
