# 08_step4_ml.R -------------------------------------------------------------
# STEP 4 : LASSO + SVM-RFE 교집합 -> ROC (원본 machine learning modeling 1·2.R 로직)
#   두 갈래를 분리 실행/저장:
#     A(재현)  = interGenes 63 (Supp4 DEG ∩ MR-union)     -> results/step4_A_repro
#     B(탐색)  = 전체 DEG union(Supp3/4/5), MR필터 없음    -> results/step4_B_explore
#   각 갈래: LASSO(cv.glmnet, lambda.min) ∩ SVM-RFE(caret rfe, svmRadial) = interGenes.txt
#            그 후 train/test ROC-AUC 계산. FN1/ALDH2 는 교집합 여부와 무관하게 AUC 별도 산출.
# 핵심 재현 포인트: 교집합에 FN1·ALDH2 가 살아남는가.
# ---------------------------------------------------------------------------

suppressMessages({
  library(glmnet); library(e1071); library(kernlab); library(caret)
  library(pROC); library(VennDiagram); library(limma)
})
library(here)
source(here::here("config.R"))

TRAIN <- read.table(file.path(OUT_DIR, "data.train.txt"), header = TRUE, sep = "\t",
                    check.names = FALSE, row.names = 1); TRAIN <- as.matrix(TRAIN)
TEST  <- read.table(file.path(OUT_DIR, "data.test.txt"),  header = TRUE, sep = "\t",
                    check.names = FALSE, row.names = 1); TEST  <- as.matrix(TEST)
label_of <- function(cn) ifelse(grepl("_Control$", cn), "Control", "DKD")

auc_of <- function(mat, gene) {
  if (!(gene %in% rownames(mat))) return(NA_real_)
  y <- ifelse(grepl("_Control$", colnames(mat)), 0, 1)
  as.numeric(pROC::auc(pROC::roc(y, as.numeric(mat[gene, ]), quiet = TRUE)))
}

run_branch <- function(tag, listFile, outDir) {
  message("\n######## BRANCH ", tag, " ########")
  cand <- scan(listFile, what = "character", quiet = TRUE)
  feats <- intersect(cand, rownames(TRAIN))
  message("[", tag, "] 후보 ", length(cand), " | train 피처 ", length(feats))

  rt <- t(TRAIN[feats, , drop = FALSE])          # 행=샘플, 열=유전자
  grp <- label_of(rownames(rt))

  ## ---- LASSO ----
  set.seed(123)
  x <- as.matrix(rt); y <- factor(grp)
  fit <- glmnet(x, y, family = "binomial", alpha = 1)
  cvfit <- cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = "deviance", nfolds = 10)
  co <- coef(fit, s = cvfit$lambda.min)
  lassoGene <- rownames(co)[which(co != 0)]; lassoGene <- setdiff(lassoGene, "(Intercept)")
  writeLines(lassoGene, file.path(outDir, "LASSO.gene.txt"))
  message("[", tag, "] LASSO 선택 ", length(lassoGene), " genes")

  ## ---- SVM-RFE (caret rfe, svmRadial) ----
  set.seed(123)
  grpF <- factor(grp, levels = c("Control", "DKD"))
  ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10, verbose = FALSE)
  Profile <- tryCatch(
    rfe(x = rt, y = as.numeric(grpF),
        sizes = c(2, 3, 4, 5, 6, 7, 8),
        rfeControl = ctrl, method = "svmRadial",
        preProcess = c("center", "scale")),
    error = function(e) { message("[", tag, "] SVM-RFE 오류: ", conditionMessage(e)); NULL })
  svmGene <- if (is.null(Profile)) character(0) else Profile$optVariables
  writeLines(svmGene, file.path(outDir, "SVM-RFE.gene.txt"))
  message("[", tag, "] SVM-RFE 선택 ", length(svmGene), " genes")

  ## ---- 교집합 ----
  interGenes <- intersect(lassoGene, svmGene)
  writeLines(interGenes, file.path(outDir, "interGenes.txt"))
  fn1_in <- "FN1" %in% interGenes; aldh2_in <- "ALDH2" %in% interGenes
  message("[", tag, "] 교집합 ", length(interGenes), " genes | FN1:", fn1_in, " ALDH2:", aldh2_in)
  message("[", tag, "] 교집합 목록: ", paste(interGenes, collapse = ", "))

  ## Venn
  if (length(lassoGene) > 0 && length(svmGene) > 0) {
    vp <- venn.diagram(list(LASSO = lassoGene, `SVM-RFE` = svmGene), filename = NULL,
                       fill = c("cornflowerblue", "darkorchid1"), alpha = 0.5,
                       scaled = FALSE, cat.cex = 1.1, cex = 1.2, margin = 0.1)
    pdf(file.path(outDir, "venn.pdf"), width = 5, height = 5); grid::grid.draw(vp); dev.off()
  }

  ## ---- ROC-AUC (train/test) : 교집합 유전자 + FN1/ALDH2 항상 ----
  rocGenes <- unique(c(interGenes, "FN1", "ALDH2"))
  aucTab <- data.frame(
    gene = rocGenes,
    inIntersect = rocGenes %in% interGenes,
    AUC_train = vapply(rocGenes, function(g) auc_of(TRAIN, g), 0.0),
    AUC_test  = vapply(rocGenes, function(g) auc_of(TEST,  g), 0.0),
    row.names = NULL)
  write.csv(aucTab, file.path(outDir, "ROC_AUC.csv"), row.names = FALSE)

  # FN1/ALDH2 개별 ROC pdf (train/test)
  for (g in c("FN1", "ALDH2")) {
    for (nm in c("train", "test")) {
      mat <- if (nm == "train") TRAIN else TEST
      if (!(g %in% rownames(mat))) next
      yy <- ifelse(grepl("_Control$", colnames(mat)), 0, 1)
      r <- pROC::roc(yy, as.numeric(mat[g, ]), quiet = TRUE)
      pdf(file.path(outDir, paste0(nm, "_ROC.", g, ".pdf")), width = 3.5, height = 3.5)
      plot(r, print.auc = TRUE, col = "red", legacy.axes = TRUE,
           main = paste0(g, " (", tag, "-", nm, ")"))
      dev.off()
    }
  }
  message("[", tag, "] AUC(FN1) train=", round(auc_of(TRAIN,"FN1"),3), " test=", round(auc_of(TEST,"FN1"),3),
          " | AUC(ALDH2) train=", round(auc_of(TRAIN,"ALDH2"),3), " test=", round(auc_of(TEST,"ALDH2"),3))

  list(tag = tag, nCand = length(cand), nFeat = length(feats),
       nLasso = length(lassoGene), nSVM = length(svmGene),
       nInter = length(interGenes), interGenes = interGenes,
       FN1_in = fn1_in, ALDH2_in = aldh2_in,
       AUC = aucTab[aucTab$gene %in% c("FN1","ALDH2"), ])
}

A <- run_branch("A", file.path(RES_DIR, "step4_A_repro", "interGenes.List.txt"),
                file.path(RES_DIR, "step4_A_repro"))
B <- run_branch("B", file.path(RES_DIR, "step4_B_explore", "interGenes.List.txt"),
                file.path(RES_DIR, "step4_B_explore"))

## ---- A/B 비교표 ----
message("\n================ A/B 비교 ================")
cmp <- data.frame(
  branch = c("A(재현,63∩MR)", "B(탐색,DEG전체)"),
  cand = c(A$nCand, B$nCand), feat = c(A$nFeat, B$nFeat),
  LASSO = c(A$nLasso, B$nLasso), SVM = c(A$nSVM, B$nSVM),
  inter = c(A$nInter, B$nInter),
  FN1 = c(A$FN1_in, B$FN1_in), ALDH2 = c(A$ALDH2_in, B$ALDH2_in))
print(cmp, row.names = FALSE)
write.csv(cmp, file.path(RES_DIR, "step4_AB_compare.csv"), row.names = FALSE)
message("\nA 교집합: ", paste(A$interGenes, collapse=", "))
message("B 교집합: ", paste(B$interGenes, collapse=", "))
cat("\n== FN1/ALDH2 AUC (A) ==\n"); print(A$AUC, row.names = FALSE)
cat("\n== FN1/ALDH2 AUC (B) ==\n"); print(B$AUC, row.names = FALSE)
message("\n[STEP4] 완료")
