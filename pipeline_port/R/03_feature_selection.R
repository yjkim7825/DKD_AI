# 03_feature_selection.R --------------------------------------------------
# 원본: machine learning modeling 1.R (LASSO + SVM-RFE + Venn 교집합)
suppressMessages({
  library(glmnet); library(e1071); library(kernlab); library(caret); library(VennDiagram)
})
if (!exists("DATA_DIR")) source("config.R")

# train : '{샘플}_{그룹}' 라벨 학습 발현행렬
# genes : 후보 유전자 리스트 파일(없으면 전체 유전자 사용)
run_feature_selection <- function(train, genes = NULL,
                                  outdir = file.path(RESULT_DIR, "ml")) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  set.seed(RANDOM_SEED)
  rt <- read.table(train, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  if (!is.null(genes)) {
    g <- scan(genes, what = "character", quiet = TRUE)
    rt <- rt[rownames(rt) %in% g, ]
  }
  x <- t(as.matrix(rt))
  y <- factor(gsub(paste0(".*_(", CONTROL_LABEL, "|", CASE_LABEL, ")$"), "\\1", rownames(x)),
              levels = c(CONTROL_LABEL, CASE_LABEL))

  ## LASSO (cv.glmnet, lambda.min) --------------------------------
  fit <- glmnet(x, y, family = "binomial", alpha = 1)
  cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1,
                     type.measure = "deviance", nfolds = CV_FOLDS)
  pdf(file.path(outdir, "cvfit.pdf"), width = 6, height = 5.5); plot(cvfit); dev.off()
  coef <- coef(fit, s = cvfit$lambda.min)
  lassoGene <- rownames(coef)[which(coef != 0)]
  lassoGene <- setdiff(lassoGene, "(Intercept)")
  writeLines(lassoGene, file.path(outdir, "LASSO.gene.txt"))
  message("[LASSO] ", length(lassoGene), " genes")

  ## SVM-RFE (caret::rfe) -----------------------------------------
  Profile <- rfe(x = x, y = as.numeric(y),
                 sizes = SVM_RFE_SIZES,
                 rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
                 methods = "svmRadial",
                 preProcess = c("center", "scale"))
  pdf(file.path(outdir, "SVM-RFE.pdf"), width = 6, height = 5.5)
  plot(Profile$results$Variables, Profile$results$RMSE, type = "o",
       xlab = "Variables", ylab = "RMSE (CV)", col = "darkgreen"); dev.off()
  svmGene <- Profile$optVariables
  writeLines(svmGene, file.path(outdir, "SVM-RFE.gene.txt"))
  message("[SVM-RFE] ", length(svmGene), " genes")

  ## Venn 교집합 --------------------------------------------------
  geneList <- list(LASSO = unique(lassoGene), `SVM-RFE` = unique(svmGene))
  vp <- venn.diagram(geneList, filename = NULL,
                     fill = c("cornflowerblue", "darkorchid1"), alpha = .5,
                     cat.col = c("cornflowerblue", "darkorchid1"))
  pdf(file.path(outdir, "venn.pdf"), width = 5, height = 5); grid::grid.draw(vp); dev.off()
  inter <- Reduce(intersect, geneList)
  writeLines(inter, file.path(outdir, "interGenes.txt"))
  message("[Venn] LASSO ∩ SVM-RFE = ", length(inter), " genes: ",
          paste(inter, collapse = ", "))
  invisible(inter)
}

# ---- 실행 예시 ----
# run_feature_selection(file.path(DATA_DIR, "data.train.txt"),
#                       genes = file.path(DATA_DIR, "interGenes.List.txt"))
