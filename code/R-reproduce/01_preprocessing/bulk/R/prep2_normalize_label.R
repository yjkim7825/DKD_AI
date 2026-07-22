# ============================================================================
# prep2_normalize_label.R  — 전처리 2: log2 자동판정 + 정규화 + 그룹 라벨
#   · log2 필요 여부 자동 판정(값이 크면 raw로 보고 log2)
#   · normalizeBetweenArrays 로 샘플간 분포 통일
#   · 그룹 라벨: GEO 메타(어레이) 또는 열이름 접두사(RNA-seq)
#   ※ 함수 정의만 — 실행은 run_preprocess.R 에서
# ============================================================================
suppressMessages({ library(limma); library(GEOquery) })

# log2 + 정규화
prep2_normalize <- function(mat) {
  qx <- as.numeric(quantile(as.matrix(mat), c(0, .25, .5, .75, .99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)   # 아직 raw면 TRUE
  if (LogC) mat <- log2(mat + 1)
  normalizeBetweenArrays(as.matrix(mat))
}

# 그룹 라벨 ① series_matrix 메타데이터로 (오프라인, 인터넷 불필요)
#   로컬 series_matrix.txt.gz 의 !Sample_title/source/characteristics 에서 diabet→DKD
label_by_geo <- function(gse, gsm_ids, seriesDir) {
  grp <- setNames(rep("NA", length(gsm_ids)), gsm_ids)
  smf <- file.path(seriesDir, paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(smf)) { message("  [", gse, "] series_matrix 없음 → 라벨 NA"); return(grp) }
  con <- gzfile(smf, "rt"); lines <- readLines(con); close(con)
  # GSM 순서(표 헤더 ID_REF 줄)
  hdr <- grep('^"?ID_REF', lines, value = TRUE)[1]
  gsm_order <- gsub('"', '', strsplit(hdr, "\t")[[1]])[-1]
  # 그룹 판단용 텍스트 줄들 합치기 (title/source/characteristics)
  metaLines <- grep("!Sample_title|!Sample_source_name|!Sample_characteristics", lines, value = TRUE)
  disease <- rep("", length(gsm_order))
  for (ln in metaLines) {
    vals <- gsub('"', '', strsplit(ln, "\t")[[1]])[-1]
    if (length(vals) == length(gsm_order)) disease <- paste(disease, tolower(vals))
  }
  g <- ifelse(grepl("diabet|nephropathy|\\bdn\\b|dkd", disease), "DKD", "Control")
  names(g) <- gsm_order
  common <- intersect(gsm_ids, names(g)); grp[common] <- g[common]
  grp
}

# 그룹 라벨 ② 열이름 접두사로 (GSE142025용): A=Late, B=Early, N=Control
label_by_prefix <- function(cols) {
  sapply(cols, function(c) {
    p <- substr(sub("^[0-9]*", "", toupper(c)), 1, 1)
    if (p == "A") "Late" else if (p == "B") "Early" else if (p == "N") "Control" else "NA"
  })
}

# 라벨을 열이름에 붙이기: GSM123 → GSM123_DKD
apply_labels <- function(mat, grp) { colnames(mat) <- paste0(colnames(mat), "_", grp); mat }
