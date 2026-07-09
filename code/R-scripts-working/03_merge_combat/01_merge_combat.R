# 01_merge_combat.R ---------------------------------------------------------
# STEP 3 : 데이터셋 병합 + ComBat 배치보정.
#   ▶ 병합 로직은 원본 `data preprocessing 3.R` 를 1:1 로 따름:
#       공통 유전자 교집합 → 데이터셋명 접두사로 cbind → batch=데이터셋 → ComBat(par.prior=TRUE)
#       열이름 = {GSE}_{GSM}_{Control|DKD}, 출력 = *.preNorm.txt(전) / *.txt(후).
#
#   ▶ 세트 구성은 **논문 Figure 4 설계로 통일**(step4_paper 가 정본):
#       · 훈련(paper)   = GSE96804 단독(사구체, 단일 데이터셋 → ComBat 불필요) = data.train.paper.txt
#       · 검증(paper)   = ComBat(GSE104948 + GSE104954)                       = data.valid.paper.txt  ← STEP4 사용
#     (원본 ML 스크립트는 data.train/test.txt 를 씀 → 원본 네이밍 호환으로
#       data.train=ComBat(96804+104948), data.test=ComBat(104954+30529) 도 참고 생성.)
#     논문(Figure4) vs 원본(ML 경로 GSE96804_104948) 이 서로 다른 부분 → 논문을 따름(DIFF 표에 명시).
#
#   중간산출물(ComBat 전후 요약·그룹수·대조 CSV) → results/step3_merge/.
#   출력 매트릭스 → ../data/processed/. 원본 데이터·스크립트 무수정.
# ---------------------------------------------------------------------------

suppressMessages({ library(limma); library(sva) })
library(here)
source(here::here("config.R"))

S3_DIR <- file.path(RES_DIR, "step3_merge"); dir.create(S3_DIR, showWarnings = FALSE, recursive = TRUE)

read_labeled <- function(path) {
  rt <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  m <- as.matrix(rt); storage.mode(m) <- "double"; avereps(m)
}
grp_of  <- function(cn) sub(".*_(Control|DKD)$", "\\1", cn)
sdist   <- function(m) apply(m, 2, median, na.rm = TRUE)   # 샘플별 중앙값(정규화 확인)

SUMM <- list(); SAMP <- list()

# 원본 data preprocessing 3.R 와 동일한 병합+ComBat. single=TRUE 면 단일 데이터셋(ComBat 생략).
merge_combat <- function(sets, outTag, single = FALSE) {
  message("\n===== ", outTag, " : ", paste(names(sets), collapse = " + "), " =====")
  mats <- lapply(names(sets), function(id) {
    m <- read_labeled(file.path(OUT_DIR, sets[[id]]))
    colnames(m) <- paste0(id, "_", colnames(m))
    message("  [", id, "] ", nrow(m), " x ", ncol(m)); m
  })
  names(mats) <- names(sets)
  common <- Reduce(intersect, lapply(mats, rownames))
  allTab <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
  batch  <- rep(seq_along(mats), vapply(mats, ncol, 0L))
  message("  병합: ", nrow(allTab), " x ", ncol(allTab), " (batch ", length(unique(batch)), ")")

  # ComBat 전(preNorm) 저장 + 요약
  write.table(rbind(geneNames = colnames(allTab), allTab),
              file.path(OUT_DIR, paste0(outTag, ".preNorm.txt")), sep = "\t", quote = FALSE, col.names = FALSE)
  pre_med <- sdist(allTab)

  if (single) {
    out <- allTab                       # 단일 데이터셋 → ComBat 생략(배치 없음)
  } else {
    # 배치 내 분산 0 유전자 제거(ComBat 안정화). batch 는 열(샘플) 기준이라 행 제거와 무관.
    keep <- apply(allTab, 1, function(v) all(tapply(v, batch, function(z) stats::sd(z) > 0)))
    if (any(!keep)) message("  배치 내 무분산 유전자 제거: ", sum(!keep))
    out <- ComBat(dat = allTab[keep, , drop = FALSE], batch = batch, par.prior = TRUE)
  }

  outFile <- file.path(OUT_DIR, paste0(outTag, ".txt"))
  write.table(rbind(geneNames = colnames(out), out), outFile, sep = "\t", quote = FALSE, col.names = FALSE)

  grp <- grp_of(colnames(out)); post_med <- sdist(out)
  message("  그룹: ", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), " -> ", basename(outFile))

  SUMM[[outTag]] <<- data.frame(
    set = outTag, datasets = paste(names(sets), collapse = "+"),
    genes = nrow(out), samples = ncol(out),
    control = sum(grp == "Control"), dkd = sum(grp == "DKD"), batches = length(sets),
    combat = !single,
    FN1_mean = if ("FN1" %in% rownames(out)) round(mean(out["FN1", ]), 3) else NA,
    ALDH2_mean = if ("ALDH2" %in% rownames(out)) round(mean(out["ALDH2", ]), 3) else NA,
    log2_min = round(min(out), 2), log2_max = round(max(out), 2), stringsAsFactors = FALSE)
  SAMP[[outTag]] <<- data.frame(set = outTag, sample = colnames(out), group = grp,
    median_preComBat = round(pre_med, 3), median_postComBat = round(post_med, 3), stringsAsFactors = FALSE)
  invisible(out)
}

## ---- 논문 Figure 4 설계 (정본) ----
# 훈련 = GSE96804 단독 (ComBat 생략)
merge_combat(c(GSE96804 = "GSE96804.labeled.txt"), "data.train.paper", single = TRUE)
# 검증 = ComBat(GSE104948 + GSE104954)
merge_combat(c(GSE104948 = "GSE104948.labeled.txt", GSE104954 = "GSE104954.labeled.txt"), "data.valid.paper")

## ---- 원본 ML 네이밍 호환 (참고: data.train/test.txt) ----
merge_combat(c(GSE96804 = "GSE96804.labeled.txt",  GSE104948 = "GSE104948.labeled.txt"), "data.train")  # 사구체
merge_combat(c(GSE104954 = "GSE104954.labeled.txt", GSE30529  = "GSE30529.labeled.txt"),  "data.test")   # 세뇨관

## ---- 중간산출물 + 대조 CSV ----
summTab <- do.call(rbind, SUMM)
write.csv(summTab, file.path(S3_DIR, "merge_combat_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, SAMP), file.path(S3_DIR, "sample_median_pre_post_ComBat.csv"), row.names = FALSE)

# 논문 대조: 검증셋(paper) 그룹수. 논문 = 104948(C21/D12)+104954(C21/D17) = Control 42 / DKD 29
vp <- summTab[summTab$set == "data.valid.paper", ]
cmp <- data.frame(
  metric = c("valid(paper) Control", "valid(paper) DKD", "valid(paper) total",
             "train(paper) samples"),
  paper  = c(42, 29, 71, 61),
  ours   = c(vp$control, vp$dkd, vp$samples, summTab$samples[summTab$set == "data.train.paper"]))
cmp$match <- cmp$paper == cmp$ours
write.csv(cmp, file.path(S3_DIR, "compare_paper_vs_ours.design.csv"), row.names = FALSE)
cat("\n== 논문 설계 대조 ==\n"); print(cmp, row.names = FALSE)
message("\n[STEP3] 완료. 중간산출물 -> ", S3_DIR)
