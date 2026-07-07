# 04_roc.R ----------------------------------------------------------------
# 원본: machine learning modeling 2.R (pROC 유전자별 ROC/AUC + bootstrap CI)
suppressMessages(library(pROC))
if (!exists("DATA_DIR")) source("config.R")

# expr  : 발현행렬('{샘플}_{그룹}' 라벨)
# genes : 유전자 리스트 파일(LASSO.gene.txt)
# title : "Training" / "Validation"
run_roc <- function(expr, genes, title = "Training",
                    outdir = file.path(RESULT_DIR, "roc")) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  rt <- read.table(expr, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  y <- ifelse(grepl(paste0("_", CONTROL_LABEL, "$"), colnames(rt)), 0, 1)
  geneRT <- readLines(genes)

  summ <- data.frame()
  for (g in geneRT) {
    if (!(g %in% rownames(rt))) next
    roc1 <- roc(y, as.numeric(rt[g, ]), quiet = TRUE)
    ci1 <- as.numeric(ci.auc(roc1, method = "bootstrap"))
    pdf(file.path(outdir, paste0(title, "_ROC.", g, ".pdf")), width = 3.5, height = 3.5)
    plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE,
         main = paste0(g, " of ", title))
    text(0.34, 0.38, paste0("95% CI: ", sprintf("%.3f", ci1[1]), "-", sprintf("%.3f", ci1[3])),
         col = "red"); dev.off()
    summ <- rbind(summ, data.frame(gene = g, AUC = round(as.numeric(roc1$auc), 3),
                                   CI_low = round(ci1[1], 3), CI_high = round(ci1[3], 3)))
  }
  write.csv(summ, file.path(outdir, paste0("roc_auc_summary_", title, ".csv")), row.names = FALSE)
  print(summ)
  invisible(summ)
}

# ---- 실행 예시 ----
# run_roc(file.path(DATA_DIR,"data.train.txt"), file.path(RESULT_DIR,"ml","LASSO.gene.txt"), "Training")
# run_roc(file.path(DATA_DIR,"data.test.txt"),  file.path(RESULT_DIR,"ml","LASSO.gene.txt"), "Validation")
