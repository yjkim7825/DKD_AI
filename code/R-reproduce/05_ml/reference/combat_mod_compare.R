# ============================================================================
# combat_mod_compare.R (reference/)  — ComBat mod 옵션 비교
#   같은 파이프라인(LASSO→ROC>0.8→4모델)을 배치보정만 바꿔 2번 실행:
#     (a) mod = NULL        — 배치만 제거 (라벨 미사용, 무난·무누수)   [01의 기본]
#     (b) mod = ~그룹        — Control/DKD 생물학적 차이 명시 보존 (라벨 사용)
#   목적: 검증 AUC가 얼마나 달라지는지, mod 사용이 성능을 부풀리는지 확인
#   출력: output/combat_mod_compare.csv
# ============================================================================
set.seed(123)
suppressMessages({
  library(glmnet); library(pROC); library(sva)
  library(randomForest); library(e1071); library(xgboost)
})
ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MLDIR <- file.path(ROOT, "R-reproduce/05_ml")
DATA  <- file.path(ROOT, "data/processed")
OUT   <- file.path(MLDIR, "output")
source(file.path(MLDIR, "R/ml_func.R"))
interGenes <- scan(file.path(MLDIR, "interGenes.List.txt"), what = "character", quiet = TRUE)
auc_of <- function(p, y) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE)))

# 한 세팅(use_mod)에 대해 전체 파이프라인 돌려 결과 반환
run_setting <- function(use_mod) {
  d  <- load_and_combat(file.path(DATA, "data.train.paper.txt"),
                        file.path(DATA, "data.valid.paper.txt"), interGenes, use_mod = use_mod)
  sc <- scale_train_test(d$Xtr, d$Xte); Xtr <- sc$tr; Xte <- sc$te
  ytr <- d$ytr; yte <- d$yte; genes <- d$genes
  # LASSO
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = 10)
  co <- coef(cv, s = "lambda.min"); lassoGene <- setdiff(rownames(co)[which(as.numeric(co) != 0)], "(Intercept)")
  # ROC>0.8 필터
  final <- character()
  for (g in lassoGene) {
    at <- gene_auc(ytr, Xtr[, g]); av <- gene_auc(yte, Xte[, g])   # gene_auc: 방향 무관(max)
    if (at > 0.8 && av > 0.8) final <- c(final, g)
  }
  # 4모델 결합 검증 AUC
  dtr <- data.frame(Xtr[, final, drop = FALSE], y = ytr); dte <- data.frame(Xte[, final, drop = FALSE], y = yte)
  gm <- glm(y ~ ., data = dtr, family = binomial)
  rf <- randomForest(y ~ ., data = dtr)
  sv <- svm(y ~ ., data = dtr, kernel = "radial", probability = TRUE)
  yb <- as.numeric(ytr) - 1; mtr <- as.matrix(dtr[, final, drop = FALSE]); mte <- as.matrix(dte[, final, drop = FALSE])
  xg <- xgb.train(list(objective = "binary:logistic", eval_metric = "logloss"),
                  xgb.DMatrix(mtr, label = yb), nrounds = 50, verbose = 0)
  data.frame(
    mod = if (use_mod) "~그룹" else "NULL",
    final_genes = paste(final, collapse = "+"),
    GLM = round(auc_of(predict(gm, dte, type = "response"), yte), 3),
    RF  = round(auc_of(predict(rf, dte, type = "prob")[, "DKD"], yte), 3),
    SVM = round(auc_of(attr(predict(sv, dte, probability = TRUE), "probabilities")[, "DKD"], yte), 3),
    XGB = round(auc_of(predict(xg, xgb.DMatrix(mte)), yte), 3)
  )
}

cat("\n===== ComBat mod 비교 (검증셋 AUC) =====\n")
res <- rbind(run_setting(FALSE), run_setting(TRUE))
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "combat_mod_compare.csv"), row.names = FALSE)
cat("\n저장: output/combat_mod_compare.csv\n")
cat("해석: 두 세팅의 최종 유전자·AUC가 비슷하면 mod 영향 작음(=재현 안정).\n")
cat("      ~그룹 쪽이 크게 높으면 라벨 사용에 의한 성능 부풀림 의심 → 기본은 NULL 권장.\n")
