# 03_GSE96804_to_matrix.R ---------------------------------------------------
# STEP 2 (CEL -> RMA) : GSE96804 (훈련 bulk, 사구체)
#   플랫폼: Affymetrix HTA-2_0 (whole-transcript, GPL17586)
#   방식  : oligo::read.celfiles -> rma(target="core") -> transcript cluster ID
#           -> getNetAffx 로 gene symbol 주석 -> limma::avereps 로 심볼당 평균
#   그룹  : 파일명이 전부 'DN' 이라 GEOquery::getGEO 메타데이터로 대조군/DKD 구분
#   출력  : ../data/processed/GSE96804.labeled.txt  (열이름 = {GSM}_{group})
# 원본 저장소에는 CEL->RMA 단계가 없음(이미 만들어진 매트릭스부터 시작). 이 스크립트는
# 마이크로어레이 관례(oligo RMA)에 따라 새로 구현한 것 — 원본 1:1 대조 대상 아님.
# ---------------------------------------------------------------------------

suppressMessages({
  library(oligo)
  library(limma)
  library(GEOquery)
})

library(here)
source(here::here("config.R"))

celDir <- DIR_GSE96804
celFiles <- list.files(celDir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
message("[GSE96804] CEL 파일 수: ", length(celFiles))
stopifnot(length(celFiles) > 0)

# ---- 1) CEL 읽기 + RMA (core = gene/transcript cluster 레벨) ----
raw <- read.celfiles(celFiles)
eset <- oligo::rma(raw, target = "core")
expr <- exprs(eset)                       # 행=transcript cluster ID, 열=파일명
message("[GSE96804] RMA 후 차원: ", nrow(expr), " x ", ncol(expr))

# 열이름을 GSM ID 로 정리 (예: GSM2544275_DN_01.CEL.gz -> GSM2544275)
gsmFromFile <- function(x) sub("^(GSM[0-9]+).*", "\\1", basename(x))
colnames(expr) <- gsmFromFile(colnames(expr))

# ---- 2) transcript cluster -> gene symbol (netaffx 주석) ----
anno <- getNetAffx(eset, "transcript")    # pd.hta.2.0 내장 주석
fd <- pData(anno)
# geneassignment 형식: "NM_xxx // SYMBOL // desc // ..."  (첫 심볼만 사용)
extractSymbol <- function(ga) {
  ga <- as.character(ga)
  out <- rep(NA_character_, length(ga))
  ok <- !is.na(ga) & ga != "---" & ga != ""
  parts <- strsplit(ga[ok], "\\s*//\\s*")
  out[ok] <- vapply(parts, function(p) if (length(p) >= 2) p[2] else NA_character_, character(1))
  out
}
sym <- extractSymbol(fd$geneassignment)
names(sym) <- rownames(fd)

# expr 행이름(transcript cluster ID)에 맞춰 심볼 매핑
symVec <- sym[rownames(expr)]
keep <- !is.na(symVec) & symVec != ""
expr2 <- expr[keep, , drop = FALSE]
symVec <- symVec[keep]
message("[GSE96804] 심볼 매핑된 프로브: ", nrow(expr2), " (전체 ", nrow(expr), ")")

# ---- 3) 심볼당 평균 (avereps) ----
exprSym <- avereps(expr2, ID = symVec)
message("[GSE96804] 유니크 심볼 수: ", nrow(exprSym))

# ---- 4) 그룹 라벨 (GEOquery 메타데이터) ----
grp <- rep(NA_character_, ncol(exprSym))
names(grp) <- colnames(exprSym)
gseMeta <- tryCatch({
  gse <- getGEO("GSE96804", GSEMatrix = TRUE, getGPL = FALSE, destdir = OUT_DIR)
  pData(gse[[1]])
}, error = function(e) { message("[GSE96804] getGEO 실패: ", conditionMessage(e)); NULL })

if (!is.null(gseMeta)) {
  # 특성(characteristics) 텍스트를 모아 DKD/대조군 판정
  #   DKD  = source_name "diabetic nephropathy" (41)
  #   대조 = "unaffected portion of tumor nephrectomies" (20)
  # → 'diabet' 포함 여부로 양성 정의(대조군은 정상/nephrectomy 라 이름이 다양함)
  charCols <- grep("characteristics|source_name|title", names(gseMeta), ignore.case = TRUE, value = TRUE)
  txt <- apply(gseMeta[, charCols, drop = FALSE], 1, function(r) paste(tolower(r), collapse = " | "))
  names(txt) <- rownames(gseMeta)
  isDKD <- grepl("diabet", txt)
  metaGrp <- ifelse(isDKD, GROUP_DKD, GROUP_CONTROL)
  names(metaGrp) <- names(txt)      # grepl/ifelse 가 이름을 버리므로 재부여
  common <- intersect(colnames(exprSym), names(metaGrp))
  grp[common] <- metaGrp[common]
  message("[GSE96804] 메타 매칭 샘플: ", length(common), " / ", ncol(exprSym))
}
# 메타 실패/누락 시 표시
grp[is.na(grp)] <- "NA"
message("[GSE96804] 그룹 분포: ", paste(names(table(grp)), table(grp), sep="=", collapse=", "))

# ---- 5) 출력 (열이름 = GSM_group) ----
colnames(exprSym) <- paste0(colnames(exprSym), "_", grp)
outTab <- cbind(geneNames = rownames(exprSym), as.data.frame(exprSym, check.names = FALSE))
outFile <- file.path(OUT_DIR, "GSE96804.labeled.txt")
write.table(outTab, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
message("[GSE96804] 저장 완료 -> ", outFile)
