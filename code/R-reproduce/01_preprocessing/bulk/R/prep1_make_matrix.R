# ============================================================================
# prep1_make_matrix.R  — 전처리 1: 원본 → 유전자 발현행렬
#   · 마이크로어레이: CEL → RMA → db매핑 → avereps (probe→유전자)
#   · RNA-seq       : 샘플별 txt(Symbol+값) 병합 (probe변환 없음)
#   ※ 이 파일은 "함수 정의만" — 실제 실행은 run_preprocess.R 에서
# ============================================================================
suppressMessages({ library(affy); library(limma) })

# 마이크로어레이 1개 → 유전자행렬 (행=유전자, 열=GSM)
prep1_array <- function(celDir, cdfname, dbpkg) {
  suppressMessages(library(dbpkg, character.only = TRUE))
  cel <- list.files(celDir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(cel) > 0)
  expr <- exprs(affy::rma(ReadAffy(filenames = cel, cdfname = cdfname)))  # probe×샘플 (RMA=log2+정규화 포함)
  colnames(expr) <- sub("^(GSM[0-9]+).*", "\\1", basename(colnames(expr)))
  db  <- get(dbpkg)
  map <- AnnotationDbi::select(db, keys = rownames(expr), columns = "SYMBOL", keytype = "PROBEID")
  map <- map[!is.na(map$SYMBOL), ]; map <- map[!duplicated(map$PROBEID), ]  # probe당 심볼 1개
  sym <- map$SYMBOL[match(rownames(expr), map$PROBEID)]
  keep <- !is.na(sym)
  avereps(expr[keep, , drop = FALSE], ID = sym[keep])   # 같은 심볼 평균 → 유전자×샘플
}

# 마이크로어레이(신형 Gene/HTA 칩) → 유전자행렬 : oligo 방식
#   affy 로 못 읽는 최신 칩(HTA-2.0 등). pd패키지 + oligo::rma 사용.
prep1_array_oligo <- function(celDir, dbpkg) {
  suppressMessages({ library(oligo); library(dbpkg, character.only = TRUE) })
  cel <- list.files(celDir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
  stopifnot(length(cel) > 0)
  raw  <- oligo::read.celfiles(cel)                       # pd패키지 자동 인식
  expr <- exprs(oligo::rma(raw, target = "core"))         # transcript(core) 수준
  colnames(expr) <- sub("^(GSM[0-9]+).*", "\\1", basename(colnames(expr)))
  db  <- get(dbpkg)
  map <- AnnotationDbi::select(db, keys = rownames(expr), columns = "SYMBOL", keytype = "PROBEID")
  map <- map[!is.na(map$SYMBOL), ]; map <- map[!duplicated(map$PROBEID), ]
  sym <- map$SYMBOL[match(rownames(expr), map$PROBEID)]
  keep <- !is.na(sym)
  avereps(expr[keep, , drop = FALSE], ID = sym[keep])
}

# 마이크로어레이(Brainarray 커스텀 CDF, 칩 2종 혼합) → 유전자행렬  [저자 방식]
#   · 폴더에 U133A(712) + Plus2(1164) CEL이 섞여 있음 → 칩 크기별로 나눠 각각 RMA 후 병합
#   · 196개 중 DN(당뇨→DKD) + LD(공여자→Control)만 필터
#   · probeset(ENTREZ) → org.Hs.eg.db 로 심볼
#   cdfMap: list("712"=U133A_cdf.gz, "1164"=Plus2_cdf.gz)
#   반환: list(mat=유전자행렬, grp=그룹벡터)
prep1_array_brainarray <- function(celDir, cdfMap, chipmapCsv) {
  suppressMessages({ library(affy); library(makecdfenv)
                     library(org.Hs.eg.db); library(AnnotationDbi); library(R.utils) })
  # 1) Python이 만든 칩맵 CSV 읽기 (file, gsm, group, chip) — R이 바이너리 파싱 안 함
  cm <- read.csv(chipmapCsv, stringsAsFactors = FALSE, colClasses = "character")
  cel  <- file.path(celDir, cm$file)
  grp  <- setNames(cm$group, cm$gsm)
  dims <- setNames(cm$chip, cel)
  cat("   필터: DKD", sum(cm$group=="DKD"), "+ Control", sum(cm$group=="Control"),
      " | 칩:", paste(names(table(cm$chip)), table(cm$chip), sep="=", collapse=" "), "\n")
  stopifnot(all(file.exists(cel)))

  # 3) 크기별로 각각 RMA + ENTREZ→심볼
  process_group <- function(files, cdfGz, tag) {
    tmpcdf <- tempfile(fileext = ".cdf")
    R.utils::gunzip(cdfGz, destname = tmpcdf, overwrite = TRUE, remove = FALSE)
    envName <- paste0("brainCdf_", tag)
    assign(envName, make.cdf.env(basename(tmpcdf), cdf.path = dirname(tmpcdf)), envir = globalenv())
    raw <- ReadAffy(filenames = files); raw@cdfName <- envName
    expr <- exprs(affy::rma(raw))
    colnames(expr) <- sub("^(GSM[0-9]+).*", "\\1", basename(colnames(expr)))
    eg  <- sub("_at$", "", rownames(expr))
    sym <- mapIds(org.Hs.eg.db, keys = eg, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    keep <- !is.na(sym)
    avereps(expr[keep, , drop = FALSE], ID = sym[keep])
  }
  mats <- list()
  for (d in names(table(dims))) {
    if (is.null(cdfMap[[d]])) { message("   [경고] 크기 ", d, " 용 CDF 없음, 건너뜀"); next }
    cat("   → 크기", d, "처리 (", sum(dims==d), "개 CEL )\n")
    mats[[d]] <- process_group(cel[dims == d], cdfMap[[d]], d)
  }
  # 4) 칩별 행렬을 공통 유전자로 병합
  genes <- Reduce(intersect, lapply(mats, rownames))
  mat <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE]))
  list(mat = mat, grp = grp[colnames(mat)])
}

# RNA-seq → 유전자행렬 (샘플별 txt 병합)
prep1_rnaseq <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  stopifnot(length(files) > 0)
  read_one <- function(f) {
    s <- sub(".*_(.+)\\.txt\\.gz$", "\\1", basename(f))               # 파일명에서 샘플ID
    d <- read.delim(gzfile(f), header = TRUE, check.names = FALSE)
    setNames(data.frame(d[[1]], d[[2]]), c("Symbol", s))             # 1열=심볼, 2열=값
  }
  lst <- lapply(files, read_one)
  mat <- Reduce(function(a, b) merge(a, b, by = "Symbol"), lst)      # 심볼 기준 병합
  rownames(mat) <- mat$Symbol; mat$Symbol <- NULL
  as.matrix(mat)
}
