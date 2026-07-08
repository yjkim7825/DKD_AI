# 06_merge_combat.R ---------------------------------------------------------
# STEP 3 : 데이터셋 병합 + ComBat 배치보정 (원본 `data preprocessing 3.R` 로직)
#   훈련(data.train.txt) = ComBat(GSE96804 + GSE104948)  ← 둘 다 사구체(glomeruli)
#   검증(data.test.txt)  = ComBat(GSE104954 + GSE30529)  ← 둘 다 세뇨관(tubule)
#   근거: 원본 ML 스크립트 작업경로가 'GSE96804_104948'(train), 입력이 data.train/test.txt.
#   방식: 공통 유전자 교집합 -> 데이터셋명 접두사로 cbind -> ComBat(batch=데이터셋).
#         열이름 = {GSE}_{GSM}_{Control|DKD} (뒤 _Control/_DKD 접미사로 STEP4가 라벨 추출).
#   출력: ../data/processed/  data.train.txt / data.test.txt (+ *.preNorm.txt 참고용)
# 원본과 차이: 원본은 폴더 내 모든 txt 자동 수집. 여기선 세트별 파일 목록을 명시(재현 명확성).
# ---------------------------------------------------------------------------

suppressMessages({
  library(limma)
  library(sva)
})

library(here)
source(here::here("config.R"))

# STEP 2 산출물 (processed/*.labeled.txt) 을 세트로 묶음
TRAIN_SETS <- c(GSE96804 = "GSE96804.labeled.txt",  GSE104948 = "GSE104948.labeled.txt")  # 사구체
TEST_SETS  <- c(GSE104954 = "GSE104954.labeled.txt", GSE30529  = "GSE30529.labeled.txt")   # 세뇨관

read_labeled <- function(path) {
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  m <- as.matrix(rt)
  storage.mode(m) <- "double"
  avereps(m)   # 유전자 중복 시 평균(이미 유니크지만 안전차원)
}

merge_combat <- function(sets, outTag) {
  message("\n===== ", outTag, " : ", paste(names(sets), collapse = " + "), " =====")
  mats <- lapply(names(sets), function(id) {
    m <- read_labeled(file.path(OUT_DIR, sets[[id]]))
    colnames(m) <- paste0(id, "_", colnames(m))   # 데이터셋 접두사 부여
    message("  [", id, "] ", nrow(m), " genes x ", ncol(m), " samples")
    m
  })
  names(mats) <- names(sets)

  # 공통 유전자 교집합
  common <- Reduce(intersect, lapply(mats, rownames))
  message("  공통 유전자(교집합): ", length(common))
  allTab <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
  batch  <- rep(seq_along(mats), vapply(mats, ncol, 0L))
  message("  병합 후: ", nrow(allTab), " genes x ", ncol(allTab), " samples, batch 수 ", length(unique(batch)))

  # ComBat 전(preNorm) 저장
  preOut <- rbind(geneNames = colnames(allTab), allTab)
  write.table(preOut, file = file.path(OUT_DIR, paste0(outTag, ".preNorm.txt")),
              sep = "\t", quote = FALSE, col.names = FALSE)

  # 배치 내 분산 0 유전자 제거(ComBat 안정화)
  keep <- apply(allTab, 1, function(v) all(tapply(v, batch, function(z) stats::sd(z) > 0)))
  if (any(!keep)) message("  배치 내 무분산 유전자 제거: ", sum(!keep))
  allTab <- allTab[keep, , drop = FALSE]

  # ComBat 배치보정
  cb <- ComBat(dat = allTab, batch = batch, par.prior = TRUE)
  cbOut <- rbind(geneNames = colnames(cb), cb)
  outFile <- file.path(OUT_DIR, paste0(outTag, ".txt"))
  write.table(cbOut, file = outFile, sep = "\t", quote = FALSE, col.names = FALSE)

  grp <- sub(".*_(Control|DKD)$", "\\1", colnames(cb))
  message("  그룹 분포: ", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "))
  message("  저장 -> ", outFile)
  invisible(cb)
}

merge_combat(TRAIN_SETS, "data.train")
merge_combat(TEST_SETS,  "data.test")
message("\n[STEP3] 병합 + ComBat 완료")
