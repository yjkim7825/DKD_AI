# 05_GSE104948_104954_to_matrix.R -------------------------------------------
# STEP 2 (CEL -> RMA) : GSE104948(사구체) / GSE104954(세뇨관) 검증 bulk
#   플랫폼: 한 GSE 안에 HG-U133A + HG-U133_Plus_2 두 종 혼재
#   CDF   : 제공된 brainarray ENTREZG 커스텀 CDF (프로브 -> Entrez 유전자 직접 매핑)
#           HG-U133A       -> GPL24120_HGU133A_Hs_ENTREZG.cdf.gz
#           HG-U133_Plus_2 -> GPL22945_HGU133Plus2_Hs_ENTREZG.cdf.gz
#   방식  : 플랫폼별로 make.cdf.env + affy::rma -> Entrez(_at) 레벨 발현
#           -> 두 플랫폼을 공통 Entrez 교집합으로 cbind 병합
#           -> Entrez -> SYMBOL(org.Hs.eg.db) -> avereps
#   샘플  : 파일명 질환코드에서 DN=DKD, LD=Control 만 사용(그 외 질환 제외)
#           ※ TN(tumor nephrectomy 정상부)은 제외 — 논문 Supp Table 1 Control=21 = LD(21) 정확 일치
#             (2026-07-09 정렬: 이전엔 LD+TN=26 이었음. 근거는 DECISIONS STEP2 참조)
#   출력  : ../data/processed/GSE104948.labeled.txt , GSE104954.labeled.txt
# 원본 저장소에 CEL->RMA 단계 없음 — 마이크로어레이 관례로 새로 구현.
# 주의: CDF .gz 는 원본을 수정하지 않고 scratchpad 로 gunzip 복사 후 사용.
# ---------------------------------------------------------------------------

suppressMessages({
  library(affy)
  library(makecdfenv)
  library(limma)
  library(org.Hs.eg.db)
})

library(here)
source(here::here("config.R"))

SCRATCH <- file.path(SCRATCH_DIR, "cdf")
dir.create(SCRATCH, showWarnings = FALSE, recursive = TRUE)

# CDF 파일명 -> 플랫폼(CEL 헤더 cdfName) 매핑
CDF_FOR_CHIP <- c(
  "HG-U133A"       = "GPL24120_HGU133A_Hs_ENTREZG.cdf.gz",
  "HG-U133_Plus_2" = "GPL22945_HGU133Plus2_Hs_ENTREZG.cdf.gz"
)
ENVNAME_FOR_CHIP <- c("HG-U133A" = "brainarrayu133a", "HG-U133_Plus_2" = "brainarrayplus2")

# .gz 를 plain 으로 복사(원본 미변경) -- 이미 있으면 재사용
gunzip_copy <- function(gz, outDir) {
  out <- file.path(outDir, sub("\\.gz$", "", basename(gz)))
  if (file.exists(out) && file.info(out)$size > 0) return(out)
  ci <- gzfile(gz, "rb"); co <- file(out, "wb")
  on.exit({ close(ci); close(co) })
  repeat { b <- readBin(ci, "raw", 1e7); if (length(b) == 0) break; writeBin(b, co) }
  out
}

# 커스텀 CDF 환경을 (한 번만) 만들어 캐시
.cdfEnvCache <- new.env()
get_cdf_env <- function(chip, gseDir) {
  if (exists(chip, envir = .cdfEnvCache)) return(get(chip, envir = .cdfEnvCache))
  gz <- file.path(gseDir, CDF_FOR_CHIP[[chip]])
  message("  [CDF] ", chip, " 준비: ", basename(gz))
  plain <- gunzip_copy(gz, SCRATCH)
  env <- make.cdf.env(basename(plain), cdf.path = dirname(plain), verbose = FALSE)
  assign(chip, env, envir = .cdfEnvCache)
  env
}

# 파일명 -> 질환코드 (예: GSM2810645_H1-Glom-DN1 -> DN)
dcode_of <- function(file) {
  lab <- sub("^GSM[0-9]+_", "", sub("\\.CEL\\.gz$", "", basename(file)))
  sub("[0-9_].*$", "", sub("^[^-]+-[^-]+-", "", lab))
}
gsm_of <- function(file) sub("^(GSM[0-9]+).*", "\\1", basename(file))

process_gse <- function(gseDir, gseId, keepCodes = c("DN", "LD")) {   # TN 제외(논문 Control=21=LD)
  message("\n===== ", gseId, " =====")
  cel <- list.files(gseDir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
  dcode <- dcode_of(cel)
  sel <- dcode %in% keepCodes
  cel <- cel[sel]; dcode <- dcode[sel]
  message("[", gseId, "] DN/LD/TN CEL: ", length(cel),
          " (", paste(names(table(dcode)), table(dcode), sep = "=", collapse = ", "), ")")

  # 각 CEL 의 플랫폼(chip) 판정
  chip <- vapply(cel, function(f) affyio::read.celfile.header(f)$cdfName, character(1))

  # 플랫폼별 RMA
  exprList <- list()
  for (ch in unique(chip)) {
    idx <- chip == ch
    message("[", gseId, "] ", ch, " : ", sum(idx), " arrays RMA")
    env <- get_cdf_env(ch, gseDir)
    envName <- ENVNAME_FOR_CHIP[[ch]]
    raw <- ReadAffy(filenames = cel[idx])
    raw@cdfName <- envName
    assign(envName, env, envir = .GlobalEnv)   # getCdfInfo 가 .GlobalEnv 에서 찾음
    e <- affy::rma(raw)
    m <- exprs(e)
    colnames(m) <- gsm_of(colnames(m))
    exprList[[ch]] <- m
  }

  # 공통 Entrez(_at) 교집합으로 병합
  common <- Reduce(intersect, lapply(exprList, rownames))
  message("[", gseId, "] 공통 Entrez probeset: ", length(common),
          " (플랫폼별 ", paste(vapply(exprList, nrow, 0L), collapse = "/"), ")")
  merged <- do.call(cbind, lapply(exprList, function(m) m[common, , drop = FALSE]))

  # Entrez -> SYMBOL
  entrez <- sub("_at$", "", rownames(merged))
  sym <- mapIds(org.Hs.eg.db, keys = entrez, column = "SYMBOL",
                keytype = "ENTREZID", multiVals = "first")
  keep <- !is.na(sym) & sym != ""
  merged <- merged[keep, , drop = FALSE]; sym <- sym[keep]
  exprSym <- avereps(merged, ID = sym)
  message("[", gseId, "] 유니크 심볼: ", nrow(exprSym), " x ", ncol(exprSym), " 샘플")

  # 그룹: DN -> DKD, LD -> Control (TN 은 keepCodes 에서 이미 제외)
  gsmAll <- gsm_of(cel)
  codeByGsm <- setNames(dcode, gsmAll)
  grp <- ifelse(codeByGsm[colnames(exprSym)] == "DN", GROUP_DKD, GROUP_CONTROL)
  message("[", gseId, "] 그룹 분포: ",
          paste(names(table(grp)), table(grp), sep = "=", collapse = ", "))

  colnames(exprSym) <- paste0(colnames(exprSym), "_", grp)
  outTab <- cbind(geneNames = rownames(exprSym), as.data.frame(exprSym, check.names = FALSE))
  outFile <- file.path(OUT_DIR, paste0(gseId, ".labeled.txt"))
  write.table(outTab, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[", gseId, "] 저장 완료 -> ", outFile)
  invisible(exprSym)
}

process_gse(DIR_GSE104948, "GSE104948")
process_gse(DIR_GSE104954, "GSE104954")
message("\n[STEP2] GSE104948/104954 완료")
