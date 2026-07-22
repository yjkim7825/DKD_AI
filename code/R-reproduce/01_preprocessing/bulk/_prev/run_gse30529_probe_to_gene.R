# ============================================================================
# run_gse30529_probe_to_gene.R
#   목적: GSE30529(마이크로어레이) CEL 원본 → RMA 정규화 → probe→유전자심볼
#         변환 → 심볼당 평균 → 발현행렬(유전자×샘플) 저장
#   방식: CEL → affy::rma → hgu133a2.db 매핑 → limma::avereps
#         (GPL 주석표 다운 불필요 — 패키지가 매핑표를 내장)
#   입력: code/data/1-5. GSE30529_RAW/*.CEL.gz
#   출력: code/R-preprocessing-run/output/GSE30529.genematrix.txt
# ============================================================================

## ── 0) 필요한 패키지 (없으면 설치) ─────────────────────────────────────────
need <- c("affy", "limma", "hgu133a2.db")   # RMA / avereps / probe매핑
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}
suppressMessages({ library(affy); library(limma); library(hgu133a2.db) })

## ── 1) 경로 설정 (★ 유정님 PC 경로) ────────────────────────────────────────
DATA_DIR <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/data/1-5. GSE30529_RAW"
OUT_DIR  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/R-preprocessing-run/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ── 2) CEL 파일 읽기 ───────────────────────────────────────────────────────
celFiles <- list.files(DATA_DIR, pattern = "\\.CEL\\.gz$",
                       full.names = TRUE, ignore.case = TRUE)
cat("■ [입력] CEL 파일 수:", length(celFiles), "\n")
stopifnot(length(celFiles) > 0)

raw <- ReadAffy(filenames = celFiles, cdfname = "hgu133a2cdf")

## ── 3) RMA 정규화 (probe 수준) ─────────────────────────────────────────────
eset <- affy::rma(raw)
expr <- exprs(eset)                       # 행 = probe, 열 = 파일명
cat("■ [RMA 후] probe×샘플:", nrow(expr), "x", ncol(expr), "\n")
cat("── RMA 결과 미리보기 (probe 5개 × 샘플 3개) ──\n")
print(round(expr[1:5, 1:3], 3))

# 열이름 정리: 파일명 → GSM ID
colnames(expr) <- sub("^(GSM[0-9]+).*", "\\1", basename(colnames(expr)))

## ── 4) probe → 유전자 심볼 (hgu133a2.db) ───────────────────────────────────
map <- AnnotationDbi::select(hgu133a2.db, keys = rownames(expr),
                             columns = "SYMBOL", keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL), ]
map <- map[!duplicated(map$PROBEID), ]     # probe당 심볼 1개
sym <- map$SYMBOL[match(rownames(expr), map$PROBEID)]
keep <- !is.na(sym)
expr2 <- expr[keep, , drop = FALSE]
sym   <- sym[keep]
cat("\n■ [매핑] 심볼 붙은 probe:", nrow(expr2), "/", nrow(expr), "\n")
cat("── 매핑 예시 (probe → 심볼) ──\n")
print(head(data.frame(probe = rownames(expr2), symbol = sym), 5))

## ── 5) 같은 심볼 여러 probe → 평균 (avereps) ───────────────────────────────
exprSym <- avereps(expr2, ID = sym)
cat("\n■ [평균 후] 유전자×샘플:", nrow(exprSym), "x", ncol(exprSym), "\n")
cat("── 최종 발현행렬 미리보기 (유전자 5개 × 샘플 3개) ──\n")
print(round(exprSym[1:5, 1:3], 3))

## ── 6) 저장 ────────────────────────────────────────────────────────────────
outMat <- cbind(geneNames = rownames(exprSym), as.data.frame(exprSym))
outFile <- file.path(OUT_DIR, "GSE30529.genematrix.txt")
write.table(outMat, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n■ [출력] 저장 완료:", outFile, "\n")
cat("   → 행 =", nrow(exprSym), "유전자 / 열 =", ncol(exprSym), "샘플\n")

# ============================================================================
# 흐름 요약:
#   CEL 원본 → (RMA) probe행렬 → (hgu133a2.db) probe→심볼 → (avereps) 유전자행렬
#   결과: output/GSE30529.genematrix.txt  (다음 단계: 정규화·그룹라벨·병합)
# ============================================================================
