# 01_preprocess.R ---------------------------------------------------------
# 원본: data preprocessing 1/2/3.R 를 config 기반 함수로 정리.
#   - probe2gene()   : GEO series matrix + 플랫폼 주석 -> gene x sample
#   - normalize_ds() : log2 자동판정 + normalizeBetweenArrays + 3군 라벨링
#   - combat_merge() : 교집합 유전자 병합 + ComBat 배치보정
# 하드코딩 경로 제거. 아래 예시 블록을 참고해 자신의 데이터 경로로 호출.

suppressMessages({
  library(readr); library(dplyr); library(stringr)
  library(limma); library(sva)
})

# BASE_DIR 미리 지정 후 source. 예:
# BASE_DIR <- "C:/Users/.../pipeline_port"; source(file.path(BASE_DIR,"R","config.R"))
if (!exists("DATA_DIR")) source("config.R")

## 1) probe -> gene ---------------------------------------------------------
# series      : GSExxxxx_series_matrix.txt(.gz)
# platform    : GPLxxxx 주석 파일
# skip_series : series matrix 표 시작 전 헤더 라인 수(자동 탐지도 지원)
# symbol_sep / symbol_index : 'AAA // BBB' 심볼에서 사용할 토큰(R 원본은 2번째=index 2)
probe2gene <- function(series, platform, out,
                       id_col = "ID", symbol_col = "Gene Symbol",
                       symbol_sep = " // ", symbol_index = 1,
                       series_skip = NULL, platform_skip = NULL) {
  Sys.setenv("VROOM_CONNECTION_SIZE" = 262144)
  # --- series matrix 표 자동 탐지 ---
  if (is.null(series_skip)) {
    con <- if (grepl("\\.gz$", series)) gzfile(series) else file(series)
    lines <- readLines(con, n = 200); close(con)
    hit <- grep("!series_matrix_table_begin", lines)
    series_skip <- if (length(hit)) hit[1] else 0
  }
  expr <- read_tsv(series, skip = series_skip, show_col_types = FALSE)
  expr <- as.data.frame(expr)
  colnames(expr)[1] <- "ID_REF"
  expr <- expr[!grepl("^!", expr$ID_REF), ]
  rownames(expr) <- expr$ID_REF
  expr$ID_REF <- NULL

  # --- 플랫폼 주석에서 probe -> symbol ---
  if (is.null(platform_skip)) {
    con <- if (grepl("\\.gz$", platform)) gzfile(platform) else file(platform)
    plines <- readLines(con, n = 200); close(con)
    hit <- grep(paste0("(^|\\t)", id_col, "(\\t|$)"), plines)
    platform_skip <- if (length(hit)) hit[1] - 1 else 0
  }
  probe <- read_tsv(platform, skip = platform_skip, show_col_types = FALSE)
  ids <- probe %>% select(all_of(c(id_col, symbol_col)))
  colnames(ids) <- c("ID", "Symbol")
  ids$Symbol <- vapply(strsplit(ids$Symbol, symbol_sep, fixed = TRUE),
                       function(x) if (length(x) >= symbol_index) trimws(x[symbol_index]) else trimws(x[1]),
                       character(1))
  ids <- ids[!is.na(ids$Symbol) & ids$Symbol != "", ]
  ids$ID <- as.character(ids$ID)

  expr$ID <- rownames(expr)
  merged <- inner_join(ids, expr, by = "ID")
  mat <- avereps(merged[, !(colnames(merged) %in% c("ID", "Symbol"))],
                 ID = merged$Symbol)
  mat <- as.data.frame(mat)
  out_tab <- cbind(geneNames = rownames(mat), mat)
  write.table(out_tab, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[probe2gene] ", nrow(mat), " genes x ", ncol(mat), " samples -> ", out)
  invisible(mat)
}

## 2) 정규화 + 3군 라벨링 ---------------------------------------------------
# control/early/late : 각 군 샘플명(한 줄 1개) 텍스트 파일 경로(없으면 NULL)
normalize_ds <- function(matrix_file, control = NULL, early = NULL, late = NULL,
                         geoid = "dataset", out = NULL) {
  rt <- read.table(matrix_file, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(rt) <- rt[, 1]; exp <- as.matrix(rt[, -1])
  data <- matrix(as.numeric(exp), nrow = nrow(exp),
                 dimnames = list(rownames(exp), colnames(exp)))
  data <- avereps(data)

  # 자동 log2
  qx <- as.numeric(quantile(data, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
  if ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)) {
    data[data < 0] <- 0; data <- log2(data + 1)
  }
  data <- normalizeBetweenArrays(data)

  rd <- function(p) if (is.null(p)) character(0) else trimws(readLines(p))
  cols <- c(); grp <- c()
  for (pl in list(c(control, "Control"), c(early, "Early"), c(late, "Late"))) {
    if (length(pl) < 2 || is.null(pl[[1]])) next
    for (s in rd(pl[[1]])) if (s %in% colnames(data)) { cols <- c(cols, s); grp <- c(grp, pl[[2]]) }
  }
  sub <- data[, cols, drop = FALSE]
  colnames(sub) <- paste0(cols, "_", grp)
  outData <- rbind(geneNames = colnames(sub), sub)
  if (is.null(out)) out <- file.path(DATA_DIR, paste0(geoid, "_threeGroups.normalize.txt"))
  write.table(outData, file = out, sep = "\t", quote = FALSE, col.names = FALSE)
  message("[normalize] ", geoid, ": ", ncol(sub), " samples -> ", out,
          "  (", paste(names(table(grp)), table(grp), sep="=", collapse=", "), ")")
  invisible(sub)
}

## 3) 교집합 병합 + ComBat --------------------------------------------------
# labeled_files : normalize_ds 산출물('{샘플}_{그룹}' 라벨) 경로 벡터
combat_merge <- function(labeled_files, out = file.path(DATA_DIR, "merge.normalize.txt")) {
  read_one <- function(f) {
    rt <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    as.matrix(rt)
  }
  mats <- lapply(labeled_files, read_one)
  inter <- Reduce(intersect, lapply(mats, rownames))
  allTab <- c(); batch <- c()
  for (i in seq_along(mats)) {
    m <- mats[[i]][inter, , drop = FALSE]
    colnames(m) <- paste0("B", i, "_", colnames(m))
    allTab <- if (i == 1) m else cbind(allTab, m)
    batch <- c(batch, rep(i, ncol(m)))
  }
  corrected <- ComBat(allTab, batch, par.prior = TRUE)
  outTab <- rbind(geneNames = colnames(corrected), corrected)
  write.table(outTab, file = out, sep = "\t", quote = FALSE, col.names = FALSE)
  message("[combat] ", nrow(corrected), " genes x ", ncol(corrected),
          " samples (", length(mats), " batches) -> ", out)
  invisible(corrected)
}

# ---- 실행 예시 (주석 해제 후 자신의 경로로 수정) --------------------------
# probe2gene(file.path(DATA_DIR,"GSE37263_series_matrix.txt.gz"),
#            file.path(DATA_DIR,"GPL5175.txt"),
#            out = file.path(DATA_DIR,"GSE37263_geneMatrix.txt"),
#            symbol_index = 2)
# normalize_ds(file.path(DATA_DIR,"GSE142025_geneMatrix.txt"),
#              control=file.path(DATA_DIR,"s1.txt"),
#              early  =file.path(DATA_DIR,"s2.txt"),
#              late   =file.path(DATA_DIR,"s3.txt"),
#              geoid="GSE142025")
# combat_merge(c(file.path(DATA_DIR,"GSE96804_labeled.txt"),
#                file.path(DATA_DIR,"GSE104948_labeled.txt")))
