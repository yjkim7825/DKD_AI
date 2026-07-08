# run_all.R — STEP 01~07 순차 실행 마스터 ---------------------------------------
# config.R 를 먼저 source 한 뒤 각 STEP 스크립트를 순서대로 실행한다.
# 경로는 here 기반이라 이 파일을 어디서 실행해도 동작(루트 = .git/.here 앵커).
#
# ⚠️ 주의: 일부 STEP 은 매우 무겁다(전체 재실행 시 시간·메모리 큼).
#   - STEP2 RMA(수 분), STEP6 MR(FinnGen 2.1GB 로드), STEP7 scRNA(수 분~십수 분).
#   필요한 STEP 만 골라 돌리려면 아래 STEPS 벡터에서 해당 줄만 남기세요.
#   중간 산출물(../data/processed, results/)이 이미 있으면 앞 STEP 을 건너뛰어도 됩니다.
# --------------------------------------------------------------------------------

library(here)
source(here::here("config.R"))

STEPS <- c(
  # STEP 1 — GSE142025 DEG
  "01_deg_gse142025/01_to_matrix.R",
  "01_deg_gse142025/02_deg.R",
  # STEP 2 — CEL → RMA 발현매트릭스 (무거움)
  "02_rma_microarray/01_gse96804.R",
  "02_rma_microarray/02_gse30529.R",
  "02_rma_microarray/03_gse104948_104954.R",
  # STEP 3 — 병합 + ComBat
  "03_merge_combat/01_merge_combat.R",
  # STEP 4 — LASSO+SVM-RFE → ROC (논문 본문 설계 확정본) + ROC 그림
  "04_ml_roc/01_paper_design_lasso_svm_roc.R",
  "04_ml_roc/02_roc_plots.R",
  # STEP 5 — Hallmark GSEA + KEGG ssGSEA
  "05_gsea/01_gsea.R",
  # STEP 6 — 2-표본 MR (FinnGen 2.1GB, 무거움)
  "06_mr/01_mr.R",
  # STEP 7 — scRNA Seurat (무거움)
  "07_scrna/01_scrna.R"
)

for (s in STEPS) {
  message("\n==================== RUN: ", s, " ====================")
  source(here::here(s))
}
message("\n[run_all] 전체 STEP 완료")
