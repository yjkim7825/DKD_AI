# ============================================================================
# run_all_geo_preprocess.R
#   목적: GEO bulk 5개 데이터셋 전처리 → 유전자×샘플 발현행렬 + 그룹 라벨
#   방식:
#     · 마이크로어레이(96804·104948·104954·30529): CEL → RMA → db매핑 → avereps
#     · RNA-seq(142025): 샘플별 txt(Symbol+값) 병합 (probe변환 없음)
#   출력: code/R-preprocessing-run/output/GSE####.genematrix.txt  (열 = GSM_그룹)
#   ※ data/processed 는 건드리지 않음. 결과는 output/ 에만 저장.
# ============================================================================

## ── 0) 패키지 ──────────────────────────────────────────────────────────────
need <- c("affy","limma","GEOquery",
          "hgu133a2.db","hgu133a.db","hgu133plus2.db","hta20transcriptcluster.db")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}
suppressMessages({ library(affy); library(limma); library(GEOquery) })

## ── 1) 경로 ────────────────────────────────────────────────────────────────
DATA_ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/data"
OUT_DIR   <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/R-preprocessing-run/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ── 2) 공통 함수 ──────────────────────────────────────────────────────────

# (A) 그룹 라벨: GEOquery 메타에서 diabet/nephropathy → DKD, 아니면 Control
label_groups <- function(gse, gsm_ids) {
  grp <- setNames(rep(NA_character_, length(gsm_ids)), gsm_ids)
  meta <- tryCatch(
    Biobase::pData(getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = OUT_DIR)[[1]]),
    error = function(e) { message("  [", gse, "] getGEO 실패(오프라인?): ", conditionMessage(e)); NULL })
  if (is.null(meta)) return(grp)
  cols <- grep("characteristics|source_name|title", names(meta), ignore.case = TRUE, value = TRUE)
  txt  <- apply(meta[, cols, drop = FALSE], 1, function(r) paste(tolower(r), collapse = " | "))
  names(txt) <- rownames(meta)
  isDKD <- grepl("diabet|nephropathy|\\bdn\\b|dkd", txt)
  g <- ifelse(isDKD, "DKD", "Control"); names(g) <- names(txt)
  common <- intersect(gsm_ids, names(g)); grp[common] <- g[common]
  grp
}

# (B) 마이크로어레이 1개 처리: CEL → RMA → db매핑 → avereps → 라벨 → 저장
process_array <- function(gse, celDir, cdfname, dbpkg) {
  cat("\n=====", gse, "(마이크로어레이) =====\n")
  suppressMessages(library(dbpkg, character.only = TRUE))
  cel <- list.files(celDir, pattern="\\.CEL\\.gz$", full.names=TRUE, ignore.case=TRUE)
  cat("■ 입력 CEL:", length(cel), "\n"); if(length(cel)==0){cat("  (CEL 없음, 건너뜀)\n");return(invisible())}
  raw  <- ReadAffy(filenames = cel, cdfname = cdfname)
  expr <- exprs(affy::rma(raw))
  colnames(expr) <- sub("^(GSM[0-9]+).*", "\\1", basename(colnames(expr)))
  cat("■ RMA 후 probe×샘플:", nrow(expr), "x", ncol(expr), "\n")
  db  <- get(dbpkg)
  map <- AnnotationDbi::select(db, keys=rownames(expr), columns="SYMBOL", keytype="PROBEID")
  map <- map[!is.na(map$SYMBOL), ]; map <- map[!duplicated(map$PROBEID), ]
  sym <- map$SYMBOL[match(rownames(expr), map$PROBEID)]
  keep <- !is.na(sym); expr <- expr[keep,,drop=FALSE]; sym <- sym[keep]
  exprSym <- avereps(expr, ID = sym)
  cat("■ 유전자×샘플:", nrow(exprSym), "x", ncol(exprSym), "\n")
  grp <- label_groups(gse, colnames(exprSym))
  colnames(exprSym) <- paste0(colnames(exprSym), "_", ifelse(is.na(grp),"NA",grp))
  save_matrix(gse, exprSym)
}

# (C) RNA-seq 처리: 샘플별 txt(Symbol+값) 병합 (probe변환 없음)
process_rnaseq <- function(gse, dir) {
  cat("\n=====", gse, "(RNA-seq) =====\n")
  files <- list.files(dir, pattern="\\.txt\\.gz$", full.names=TRUE)
  cat("■ 입력 txt:", length(files), "\n"); if(length(files)==0){cat("  (없음)\n");return(invisible())}
  read_one <- function(f){ s <- sub(".*_(.+)\\.txt\\.gz$","\\1",basename(f))
    d <- read.delim(gzfile(f), header=TRUE, check.names=FALSE)
    setNames(data.frame(d[[1]], d[[2]]), c("Symbol", s)) }
  lst <- lapply(files, read_one)
  mat <- Reduce(function(a,b) merge(a,b,by="Symbol"), lst)
  rownames(mat) <- mat$Symbol; mat$Symbol <- NULL
  cat("■ 유전자×샘플:", nrow(mat), "x", ncol(mat), "\n")
  # 142025 열이름(A11A 등) 앞글자로 그룹 추정: A=Late, B=Early, N=Control (데이터 관례)
  g <- sapply(colnames(mat), function(c){
        p <- substr(sub("^[0-9]*","",toupper(c)),1,1)
        if(p=="A") "Late" else if(p=="B") "Early" else if(p=="N") "Control" else "NA")})
  colnames(mat) <- paste0(colnames(mat), "_", g)
  save_matrix(gse, as.matrix(mat))
}

save_matrix <- function(gse, m) {
  out <- cbind(geneNames = rownames(m), as.data.frame(m))
  f <- file.path(OUT_DIR, paste0(gse, ".genematrix.txt"))
  write.table(out, f, sep="\t", quote=FALSE, row.names=FALSE)
  cat("■ 저장:", f, "\n")
}

## ── 3) 5개 데이터셋 실행 ───────────────────────────────────────────────────
# 마이크로어레이 4개 (플랫폼별 cdf/db 지정)
process_array("GSE30529",  file.path(DATA_ROOT,"1-5. GSE30529_RAW"),  "hgu133a2cdf", "hgu133a2.db")
process_array("GSE96804",  file.path(DATA_ROOT,"1-1. GSE96804_RAW"),  "hta20cdf",    "hta20transcriptcluster.db")
process_array("GSE104948", file.path(DATA_ROOT,"1-2. GSE104948_RAW"), "hgu133plus2cdf","hgu133plus2.db")
process_array("GSE104954", file.path(DATA_ROOT,"1-3. GSE104954_RAW"), "hgu133plus2cdf","hgu133plus2.db")
# RNA-seq 1개
process_rnaseq("GSE142025", file.path(DATA_ROOT,"1-4. GSE142025_RAW"))

cat("\n★ 완료 — 결과는 output/ 폴더 확인\n")
# ============================================================================
# 주의:
#  · GSE96804(HTA-2.0)는 CDF가 'hta20cdf' 또는 pd패키지가 필요할 수 있음.
#    hta20cdf 설치 안 되면 96804만 별도 처리 필요(주석 참고).
#  · GSE104948/104954는 플랫폼이 U133A+Plus2 섞여있을 수 있음 → CEL이 어느 칩인지
#    확인 후 hgu133a.db 도 병행해야 할 수 있음(우선 Plus2로 시도).
# ============================================================================
