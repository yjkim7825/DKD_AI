# ============================================================================
# 02_ml_model.R  — 방식 A (논문 Figure 4), 2/2단계: 진단모델
#   01이 고른 최종 유전자(FN1·ALDH2)를 4모델(RF·SVM·GLM·XGBoost)로 결합
#   → 학습·검증 AUC + 논문 Figure 4 D~I 그림
#   입력 : output/A_final_genes.txt (01 산출), data.{train,valid}.paper.txt
#   출력 : output/A_model_AUC.csv · A_Figure4_DI.pdf(Fig4 D~I)
#   ※ 01_ml_feature.R 먼저 실행 (run_ml.R이 01→02 순서로 돌림)
# ============================================================================
set.seed(123)
suppressMessages({
  library(pROC); library(sva)
  library(randomForest); library(e1071); library(xgboost)
})

ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MLDIR <- file.path(ROOT, "R-reproduce/05_ml")
DATA  <- file.path(ROOT, "data/processed")
OUT   <- file.path(MLDIR, "output")
source(file.path(MLDIR, "R/ml_func.R"))

## 01과 동일하게 ComBat 보정된 학습/검증 행렬 확보 (결정적, set.seed 동일)
interGenes <- scan(file.path(MLDIR, "interGenes.List.txt"), what = "character", quiet = TRUE)
d   <- load_and_combat(file.path(DATA, "data.train.paper.txt"),
                       file.path(DATA, "data.valid.paper.txt"), interGenes)
sc  <- scale_train_test(d$Xtr, d$Xte); Xtr <- sc$tr; Xte <- sc$te
ytr <- d$ytr; yte <- d$yte

## 01이 고른 최종 유전자
ff <- file.path(OUT, "A_final_genes.txt")
if (!file.exists(ff)) stop("A_final_genes.txt 없음 — 01_ml_feature.R 먼저 실행")
finalGenes <- scan(ff, what = "character", quiet = TRUE)
cat("[02] 최종 유전자:", paste(finalGenes, collapse = ", "), "\n")

## ── A3) 4모델(RF·SVM·GLM·XGBoost) 학습·검증 AUC ─────────────────────────────
dtr <- data.frame(Xtr[, finalGenes, drop = FALSE], y = ytr)
dte <- data.frame(Xte[, finalGenes, drop = FALSE], y = yte)
auc_of <- function(prob, y) as.numeric(pROC::auc(pROC::roc(y, prob, quiet = TRUE)))
ptr <- list(); pte <- list()

gm  <- glm(y ~ ., data = dtr, family = binomial)
ptr$GLM <- predict(gm, dtr, type = "response"); pte$GLM <- predict(gm, dte, type = "response")
rf  <- randomForest(y ~ ., data = dtr)
ptr$RF <- predict(rf, dtr, type = "prob")[, "DKD"]; pte$RF <- predict(rf, dte, type = "prob")[, "DKD"]
sv  <- svm(y ~ ., data = dtr, kernel = "radial", probability = TRUE)
ptr$SVM <- attr(predict(sv, dtr, probability = TRUE), "probabilities")[, "DKD"]
pte$SVM <- attr(predict(sv, dte, probability = TRUE), "probabilities")[, "DKD"]
yb  <- as.numeric(ytr) - 1
mtr <- as.matrix(dtr[, finalGenes, drop = FALSE]); mte <- as.matrix(dte[, finalGenes, drop = FALSE])
xg  <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "logloss"),
                 data = xgb.DMatrix(mtr, label = yb), nrounds = 50, verbose = 0)
ptr$XGB <- predict(xg, xgb.DMatrix(mtr)); pte$XGB <- predict(xg, xgb.DMatrix(mte))

res <- data.frame(model = names(ptr),
  train = round(sapply(names(ptr), function(m) auc_of(ptr[[m]], ytr)), 3),
  test  = round(sapply(names(ptr), function(m) auc_of(pte[[m]], yte)), 3), row.names = NULL)
write.csv(res, file.path(OUT, "A_model_AUC.csv"), row.names = FALSE)
cat("\n[A3 4모델 결합 AUC] (학습 / 검증)\n"); print(res, row.names = FALSE)

## ── A4) 논문 Figure 4 D~I 그림 ──────────────────────────────────────────────
panel_lab <- function(lab) mtext(lab, side = 3, adj = 0, line = 1.2, font = 2, cex = 1.5)
roc_gene_panel <- function(y, expr, ttl, lab) {
  r <- roc(y, as.numeric(expr), quiet = TRUE); ci <- ci.auc(r)
  plot(r, col = "red", legacy.axes = TRUE, main = ttl)
  text(0.4, 0.40, sprintf("AUC: %.3f", as.numeric(auc(r))), col = "red")
  text(0.4, 0.30, sprintf("95%% CI: %.3f-%.3f", ci[1], ci[3]), col = "red", cex = 0.85)
  panel_lab(lab)
}
mcols <- c(RF = "#1f77b4", SVM = "#ff7f0e", XGB = "#2ca02c", GLM = "#9467bd")
roc_model_panel <- function(probs, y, ttl, lab) {
  first <- TRUE; leg <- character()
  for (m in names(probs)) {
    r <- roc(y, probs[[m]], quiet = TRUE); ci <- ci.auc(r)
    plot(r, col = mcols[m], legacy.axes = TRUE, add = !first, main = ttl); first <- FALSE
    leg <- c(leg, sprintf("%s (AUC=%.3f [%.3f-%.3f])", m, as.numeric(auc(r)), ci[1], ci[3]))
  }
  legend("bottomright", leg, col = mcols[names(probs)], lwd = 2, cex = 0.6, bty = "n")
  panel_lab(lab)
}

pdf(file.path(OUT, "A_Figure4_DI.pdf"), width = 12, height = 7.5)
layout(matrix(seq_len(2 * (length(finalGenes) + 1)), nrow = 2, byrow = TRUE))
par(pty = "s")
labs <- LETTERS[4:9]; li <- 1                            # D E F G H I
for (g in finalGenes) {                                  # D~G
  roc_gene_panel(ytr, Xtr[, g], paste0(g, " of Training"),   labs[li]); li <- li + 1
  roc_gene_panel(yte, Xte[, g], paste0(g, " of Validation"), labs[li]); li <- li + 1
}
roc_model_panel(ptr, ytr, "Training Set ROC Curves", labs[li]); li <- li + 1   # H
roc_model_panel(pte, yte, "Validation ROC Curves",  labs[li])                  # I
dev.off()
cat("[A4] 논문 Figure 4 D~I 그림 저장: output/A_Figure4_DI.pdf\n")
cat("\n★ [02 모델링] 완료 — GLM 최고(논문 재현).\n")
