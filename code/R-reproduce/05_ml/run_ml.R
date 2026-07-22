# ============================================================================
# run_ml.R — 05_ml 채택 파이프라인 실행 (01 선별 → 02 모델링)
#   실행: setwd("C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/R-reproduce/05_ml")
#         source("run_ml.R")
# ============================================================================
MLDIR <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code/R-reproduce/05_ml"

cat("\n########## 01. 유전자 선별 (ComBat→LASSO→ROC>0.8) ##########\n")
source(file.path(MLDIR, "01_ml_feature.R"))

cat("\n########## 02. 진단모델 (4모델 + Figure 4 D~I) ##########\n")
source(file.path(MLDIR, "02_ml_model.R"))

cat("\n==================================================\n")
cat("★ 05_ml 완료 — 결과는 output/ 폴더\n")
cat("  01: A_PCA_batch.pdf, A_cvfit.pdf, A_lasso.gene.txt, A_roc_filter.csv, A_final_genes.txt\n")
cat("  02: A_model_AUC.csv, A_Figure4_DI.pdf\n")
cat("\n[참고] 저자 원본·대조 코드는 reference/ 폴더:\n")
cat("  - machine learning modeling 1.R / 2.R (저자 원본)\n")
cat("  - ml_author.R (저자 방식 재현) · combat_mod_compare.R (ComBat 옵션 비교)\n")
cat("  비교 문서 : ../docs/05_ML_두방식비교.md\n")
