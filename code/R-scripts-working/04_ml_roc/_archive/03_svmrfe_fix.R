# 09_step4_svmrfe_fix.R -----------------------------------------------------
# STEP 4 보정: A(재현) 갈래의 SVM-RFE 를 실제 소수 피처로 축소되게 수정.
#   문제: 원본은 y 를 as.numeric() 로 바꿔 '회귀(RMSE)' RFE → 전체 45개가 최적으로 선택됨.
#   수정: (1) y 를 factor(Control/DKD) 로 두어 '분류' RFE 로 실행,
#         (2) 최적 크기를 tolerance 로 골라 성능 근접 시 '더 적은 피처' 선호(pickSizeTolerance).
#   LASSO 결과(12개)는 그대로 두고 교집합만 재계산 → FN1/ALDH2 생존·최종 유전자 수 보고.
#   산출: results/step4_A_repro/SVM-RFE.gene.v2.txt, interGenes.v2.txt (원본 v1 은 보존)
# ---------------------------------------------------------------------------

suppressMessages({ library(caret); library(kernlab); library(limma) })
library(here)
source(here::here("config.R"))

A_DIR <- file.path(RES_DIR, "step4_A_repro")
TRAIN <- as.matrix(read.table(file.path(OUT_DIR, "data.train.txt"), header = TRUE,
                              sep = "\t", check.names = FALSE, row.names = 1))

cand  <- scan(file.path(A_DIR, "interGenes.List.txt"), what = "character", quiet = TRUE)
feats <- intersect(cand, rownames(TRAIN))
rt <- t(TRAIN[feats, , drop = FALSE])
grp <- factor(ifelse(grepl("_Control$", rownames(rt)), "Control", "DKD"),
              levels = c("Control", "DKD"))

# 분류 RFE + 소수 피처 선호(tolerance 2%)
set.seed(123)
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 10, verbose = FALSE)
ctrl$functions$selectSize <- function(x, metric, maximize) {
  caret::pickSizeTolerance(x, metric = metric, tol = 2, maximize = maximize)
}
Profile <- rfe(x = rt, y = grp,
               sizes = c(2, 3, 4, 5, 6, 7, 8),
               metric = "Accuracy",
               rfeControl = ctrl,
               method = "svmRadial",
               preProcess = c("center", "scale"))

svmGene <- predictors(Profile)
cat("최적 크기(optsize):", Profile$optsize, " | 선택 피처:", length(svmGene), "\n")
cat("SVM-RFE(v2) 유전자:", paste(svmGene, collapse = ", "), "\n")
writeLines(svmGene, file.path(A_DIR, "SVM-RFE.gene.v2.txt"))

# 리샘플 성능 프로파일(크기별 Accuracy) 확인
print(Profile$results[, c("Variables", "Accuracy", "Kappa")])

# 교집합 재계산 (LASSO v1 그대로)
lasso <- readLines(file.path(A_DIR, "LASSO.gene.txt"))
inter2 <- intersect(lasso, svmGene)
writeLines(inter2, file.path(A_DIR, "interGenes.v2.txt"))

cat("\n== 교집합 재계산 (LASSO ∩ SVM-RFE.v2) ==\n")
cat("LASSO:", length(lasso), " | SVM-RFE.v2:", length(svmGene), " | 교집합:", length(inter2), "\n")
cat("교집합 목록:", paste(inter2, collapse = ", "), "\n")
cat("FN1 생존:", "FN1" %in% inter2, " | ALDH2 생존:", "ALDH2" %in% inter2, "\n")
