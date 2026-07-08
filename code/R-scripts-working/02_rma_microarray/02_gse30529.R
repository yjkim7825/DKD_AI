# 04_GSE30529_to_matrix.R ---------------------------------------------------
# STEP 2 (CEL -> RMA) : GSE30529 (DEG 검증 bulk, 세뇨관/tubulointerstitium)
#   플랫폼: Affymetrix HG-U133A_2 (3'-IVT, GPL571)
#   방식  : affy::ReadAffy -> rma() -> probeset ID -> hgu133a2.db 로 심볼 매핑
#           -> limma::avereps 로 심볼당 평균
#   그룹  : GEOquery::getGEO 메타데이터로 DKD / control(living donor) 구분
#   출력  : ../data/processed/GSE30529.labeled.txt  (열이름 = {GSM}_{group})
# 원본 저장소에 CEL->RMA 단계 없음 — 마이크로어레이 관례(affy RMA)로 새로 구현.
# ---------------------------------------------------------------------------

suppressMessages({
  library(affy)
  library(limma)
  library(hgu133a2.db)
  library(GEOquery)
})

library(here)
source(here::here("config.R"))

celDir <- DIR_GSE30529
celFiles <- list.files(celDir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
message("[GSE30529] CEL 파일 수: ", length(celFiles))
stopifnot(length(celFiles) > 0)

# ---- 1) CEL 읽기 + RMA ----
raw <- ReadAffy(filenames = celFiles, cdfname = "hgu133a2cdf")
eset <- affy::rma(raw)
expr <- exprs(eset)                       # 행=probeset, 열=파일명
message("[GSE30529] RMA 후 차원: ", nrow(expr), " x ", ncol(expr))

# 열이름 -> GSM ID (예: GSM757014_KS1-HG_U133A_2-2765.CEL.gz -> GSM757014)
gsmFromFile <- function(x) sub("^(GSM[0-9]+).*", "\\1", basename(x))
colnames(expr) <- gsmFromFile(colnames(expr))

# ---- 2) probeset -> gene symbol (hgu133a2.db) ----
probe2sym <- AnnotationDbi::select(hgu133a2.db, keys = rownames(expr),
                                   columns = "SYMBOL", keytype = "PROBEID")
probe2sym <- probe2sym[!is.na(probe2sym$SYMBOL), ]
probe2sym <- probe2sym[!duplicated(probe2sym$PROBEID), ]   # probe당 심볼 1개
symVec <- probe2sym$SYMBOL[match(rownames(expr), probe2sym$PROBEID)]
keep <- !is.na(symVec)
expr2 <- expr[keep, , drop = FALSE]
symVec <- symVec[keep]
message("[GSE30529] 심볼 매핑된 프로브: ", nrow(expr2), " (전체 ", nrow(expr), ")")

# ---- 3) 심볼당 평균 ----
exprSym <- avereps(expr2, ID = symVec)
message("[GSE30529] 유니크 심볼 수: ", nrow(exprSym))

# ---- 4) 그룹 라벨 (GEOquery) ----
grp <- rep(NA_character_, ncol(exprSym)); names(grp) <- colnames(exprSym)
gseMeta <- tryCatch({
  pData(getGEO("GSE30529", GSEMatrix = TRUE, getGPL = FALSE, destdir = OUT_DIR)[[1]])
}, error = function(e) { message("[GSE30529] getGEO 실패: ", conditionMessage(e)); NULL })

if (!is.null(gseMeta)) {
  charCols <- grep("characteristics|source_name|title", names(gseMeta), ignore.case = TRUE, value = TRUE)
  txt <- apply(gseMeta[, charCols, drop = FALSE], 1, function(r) paste(tolower(r), collapse = " | "))
  names(txt) <- rownames(gseMeta)
  # DKD 양성 정의(diabet/dn), 나머지(control/living donor)는 Control
  isDKD <- grepl("diabet|nephropathy|\\bdn\\b|dkd", txt)
  metaGrp <- ifelse(isDKD, GROUP_DKD, GROUP_CONTROL)
  names(metaGrp) <- names(txt)     # grepl/ifelse 가 이름을 버리므로 재부여
  common <- intersect(colnames(exprSym), names(metaGrp))
  grp[common] <- metaGrp[common]
  message("[GSE30529] 메타 매칭 샘플: ", length(common), " / ", ncol(exprSym))
  # 참고: source_name 원값 분포 확인용
  message("[GSE30529] source_name 분포: ",
          paste(names(table(gseMeta[["source_name_ch1"]])), table(gseMeta[["source_name_ch1"]]),
                sep = "=", collapse = " | "))
}
grp[is.na(grp)] <- "NA"
message("[GSE30529] 그룹 분포: ", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "))

# ---- 5) 출력 ----
colnames(exprSym) <- paste0(colnames(exprSym), "_", grp)
outTab <- cbind(geneNames = rownames(exprSym), as.data.frame(exprSym, check.names = FALSE))
outFile <- file.path(OUT_DIR, "GSE30529.labeled.txt")
write.table(outTab, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
message("[GSE30529] 저장 완료 -> ", outFile)
