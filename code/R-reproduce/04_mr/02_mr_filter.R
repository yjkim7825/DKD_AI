# ============================================================================
# 02_mr_filter.R  (= 저자 Mendelian randomization 2.R)
#   01이 만든 CSV를 읽어 인과 유전자 필터:
#     IVW p<0.05  &  3방법(IVW·Egger·median) OR 방향 일치  &  다면발현 p>0.05
#   입력 : output/table.MRresult.<tag>.csv , output/table.pleiotropy.<tag>.csv
#   출력 : output/IVW.filter.<tag>.csv
# ============================================================================
ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
OUT   <- file.path(ROOT, "R-reproduce/04_mr/output")

filter_one <- function(tag) {
  mrFile  <- file.path(OUT, paste0("table.MRresult.",  tag, ".csv"))
  pleFile <- file.path(OUT, paste0("table.pleiotropy.", tag, ".csv"))
  if (!file.exists(mrFile)) { message("[", tag, "] MRresult 없음 — 01 먼저 실행"); return(invisible()) }
  rt <- read.csv(mrFile, header = TRUE, check.names = FALSE)

  # ① 유전자별: 3방법 다 있고 IVW p<0.05 & OR 방향 일치
  ivw <- data.frame()
  for (g in unique(rt$exposure)) {
    gd <- rt[rt$exposure == g, ]
    if (nrow(gd) == 3) {
      if (gd[gd$method == "Inverse variance weighted", "pval"] < 0.05) {
        if (sum(gd$or > 1) == nrow(gd) | sum(gd$or < 1) == nrow(gd)) ivw <- rbind(ivw, gd)
      }
    }
  }
  # ② 다면발현 p>0.05 인 유전자만 (있을 때)
  if (file.exists(pleFile)) {
    ple <- read.csv(pleFile, header = TRUE, check.names = FALSE)
    ple <- ple[ple$pval > 0.05, ]
    ivw <- ivw[ivw$exposure %in% ple$exposure, ]
  }
  if (nrow(ivw)) write.csv(ivw, file.path(OUT, paste0("IVW.filter.", tag, ".csv")), row.names = FALSE)
  cat("[IVW필터 ", tag, "] 통과 유전자: ",
      if (nrow(ivw)) paste(unique(ivw$exposure), collapse = ", ") else "(없음)", "\n", sep = "")
}

for (tag in c("FinnGen","GCST")) filter_one(tag)
cat("\n★ [02] 필터 완료 → IVW.filter.*.csv (다음: 03_mr_forest.R)\n")
