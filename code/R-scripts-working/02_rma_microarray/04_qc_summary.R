# 04_qc_summary.R ----------------------------------------------------------
# STEP2 중간산출물 저장: 4개 RMA 발현매트릭스(processed/)의 QC·요약을 results/step2_rma/ 로.
#   전체 매트릭스는 processed/ 에 있으므로(무거운 RMA 재실행 불필요) 여기서는 요약만 생성:
#   차원·log2 스케일·그룹수·FN1/ALDH2 존재/발현·샘플별 분포 + 논문 Supp1 대조 CSV.
# ※ 원본 불변, 결과는 results/ 로만.

library(here)
source(here::here("config.R"))
suppressMessages(library(limma))

S2_DIR <- file.path(RES_DIR, "step2_rma"); dir.create(S2_DIR, showWarnings = FALSE, recursive = TRUE)

mats <- list(
  GSE96804  = "GSE96804.labeled.txt",
  GSE30529  = "GSE30529.labeled.txt",
  GSE104948 = "GSE104948.labeled.txt",
  GSE104954 = "GSE104954.labeled.txt")

# 논문 Supp Table 1 (Control / DKD)
paperN <- list(GSE96804 = c(20, 41), GSE30529 = c(12, 10),
               GSE104948 = c(21, 12), GSE104954 = c(21, 17))

summ <- list(); sampleStats <- list()
for (nm in names(mats)) {
  f <- file.path(OUT_DIR, mats[[nm]])
  if (!file.exists(f)) { message("[skip] 없음: ", f); next }
  m <- as.matrix(read.table(f, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  grp <- sub(".*_", "", colnames(m))
  nCtrl <- sum(grp == GROUP_CONTROL); nDKD <- sum(grp == GROUP_DKD)
  foc <- intersect(c("FN1", "ALDH2"), rownames(m))
  fn1  <- if ("FN1"  %in% rownames(m)) round(mean(m["FN1", ]), 3)  else NA
  aldh <- if ("ALDH2" %in% rownames(m)) round(mean(m["ALDH2", ]), 3) else NA

  summ[[nm]] <- data.frame(
    dataset = nm, genes = nrow(m), samples = ncol(m),
    control = nCtrl, dkd = nDKD,
    paper_control = paperN[[nm]][1], paper_dkd = paperN[[nm]][2],
    ctrl_match = nCtrl == paperN[[nm]][1], dkd_match = nDKD == paperN[[nm]][2],
    log2_min = round(min(m, na.rm = TRUE), 2), log2_max = round(max(m, na.rm = TRUE), 2),
    FN1_present = "FN1" %in% rownames(m), ALDH2_present = "ALDH2" %in% rownames(m),
    FN1_mean = fn1, ALDH2_mean = aldh, stringsAsFactors = FALSE)

  # 샘플별 분포(정규화 확인용): median / Q1 / Q3
  q <- apply(m, 2, quantile, probs = c(.25, .5, .75), na.rm = TRUE)
  sampleStats[[nm]] <- data.frame(dataset = nm, sample = colnames(m), group = grp,
                                  Q1 = round(q[1, ], 3), median = round(q[2, ], 3),
                                  Q3 = round(q[3, ], 3), stringsAsFactors = FALSE)
  cat(sprintf("[%s] %d genes x %d (C%d/DKD%d)  log2 %.2f-%.2f  FN1=%s ALDH2=%s\n",
              nm, nrow(m), ncol(m), nCtrl, nDKD, min(m), max(m), fn1, aldh))
}

summTab <- do.call(rbind, summ)
write.csv(summTab, file.path(S2_DIR, "RMA_matrix_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, sampleStats), file.path(S2_DIR, "RMA_sample_distribution.csv"), row.names = FALSE)

# 논문 Supp1 대조 CSV (표본수)
cmp <- summTab[, c("dataset", "control", "paper_control", "dkd", "paper_dkd", "ctrl_match", "dkd_match")]
write.csv(cmp, file.path(S2_DIR, "compare_paper_vs_ours.samples.csv"), row.names = FALSE)

cat("\n[저장] ", S2_DIR, "/ (RMA_matrix_summary.csv, RMA_sample_distribution.csv, compare_paper_vs_ours.samples.csv)\n", sep = "")
