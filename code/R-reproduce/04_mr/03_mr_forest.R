# ============================================================================
# 03_mr_forest.R  (= 저자 Mendelian randomization 3.R)
#   MR 결과 CSV(IVW 행)를 읽어 forest plot
#   입력 : output/table.MRresult.<tag>.csv  (IVW 방법만 사용)
#   출력 : output/forest_MR_<tag>.png
# ============================================================================
ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
OUT  <- file.path(ROOT, "R-reproduce/04_mr/output")
for (p in c("forestploter","grid")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressMessages({ library(forestploter); library(grid) })

forest_one <- function(tag) {
  f <- file.path(OUT, paste0("table.MRresult.", tag, ".csv"))
  if (!file.exists(f)) { message("[", tag, "] MRresult 없음 — 01 먼저 실행"); return(invisible()) }
  rt <- read.csv(f, header = TRUE, check.names = FALSE)
  d <- rt[rt$method == "Inverse variance weighted", ]           # 저자: IVW만
  if (!nrow(d)) return(invisible())

  # 표 문자열 컬럼 (저자 MR3.R 형식)
  d$` ` <- paste(rep(" ", 12), collapse = " ")
  d$`OR (95% CI)` <- ifelse(is.na(d$or), "", sprintf("%.3f (%.3f to %.3f)", d$or, d$or_lci95, d$or_uci95))
  d$pval <- ifelse(d$pval < 0.001, "<0.001", sprintf("%.3f", d$pval))
  show <- d[, c("exposure","nsnp"," ","OR (95% CI)","pval")]

  tm <- forest_theme(base_size = 14, ci_pch = 16, ci_lwd = 1.5,
                     refline_lty = "dashed", refline_col = "grey20")
  p <- forest(show, est = d$or, lower = d$or_lci95, upper = d$or_uci95,
              ci_column = 3, ref_line = 1,
              xlim = c(0, max(3, ceiling(max(d$or_uci95, na.rm = TRUE)))), theme = tm)
  png(file.path(OUT, paste0("forest_MR_", tag, ".png")), width = 1050, height = 120 + 45*nrow(d), res = 120)
  print(p); dev.off()
  cat("[forest ", tag, "] 저장: forest_MR_", tag, ".png (", nrow(d), "유전자)\n", sep = "")
}

for (tag in c("GCST","FinnGen")) forest_one(tag)

## ── 논문 Figure 3D 재현: 필터 통과 유전자를 "유의 코호트" 결과로 합쳐 forest ──
#   FinnGen 통과 7개 + GCST 통과 3개 = 논문 D의 10개
forest_figureD <- function() {
  parts <- list()
  for (tag in c("FinnGen","GCST")) {
    f <- file.path(OUT, paste0("IVW.filter.", tag, ".csv"))
    if (!file.exists(f)) next
    d <- read.csv(f, header = TRUE, check.names = FALSE)
    d <- d[d$method == "Inverse variance weighted", ]   # 유의 통과 유전자의 IVW
    if (nrow(d)) { d$cohort <- tag; parts[[tag]] <- d }
  }
  if (!length(parts)) return(invisible())
  D <- do.call(rbind, parts)
  D <- D[order(D$exposure), ]                            # 알파벳순 (논문 D와 동일)

  D$` ` <- paste(rep(" ", 12), collapse = " ")
  D$`OR (95% CI)` <- sprintf("%.3f (%.3f to %.3f)", D$or, D$or_lci95, D$or_uci95)
  D$pval_s <- ifelse(D$pval < 0.001, "<0.001", sprintf("%.3f", D$pval))
  show <- D[, c("exposure","nsnp","cohort"," ","OR (95% CI)","pval_s")]
  colnames(show) <- c("Exposure","Nsnp","Cohort"," ","OR (95% CI)","Pval")

  tm <- forest_theme(base_size = 13, ci_pch = 16, ci_lwd = 1.5,
                     refline_lty = "dashed", refline_col = "grey20")
  p <- forest(show, est = D$or, lower = D$or_lci95, upper = D$or_uci95,
              ci_column = 4, ref_line = 1, xlim = c(0, 2), theme = tm)
  png(file.path(OUT, "forest_MR_FigureD.png"), width = 1150, height = 120 + 42*nrow(D), res = 120)
  print(p); dev.off()
  cat("[forest D] 저장: forest_MR_FigureD.png (", nrow(D), "유전자 = FinnGen+GCST 통과 조합)\n", sep = "")
}
forest_figureD()

cat("\n★ [03] forest 완료 → forest_MR_*.png (+ FigureD 조합)\n")
