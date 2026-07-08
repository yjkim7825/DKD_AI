# R-scripts-working — P25 (FN1·ALDH2 DKD multi-omics) 재현 파이프라인

논문 *FN1·ALDH2 DKD multi-omics* (Renal Failure 2025)의 R 파이프라인을 STEP 1~7로 구조화한 실행본.
전사체 → DEG → LASSO/SVM 바이오마커 → ROC → GSEA → MR → scRNA 를 순서대로 재현한다.

## 폴더 구조
```
R-scripts-working/
├─ config.R                # ★ 경로/파라미터 중앙 설정 (여기 한 곳만 수정). here 기반 루트 자동탐지.
├─ run_all.R               # 01~07 순차 실행 마스터
├─ .here                   # here 패키지 루트 앵커
├─ 01_deg_gse142025/       # STEP1  GSE142025 발현매트릭스 + limma DEG
│   ├─ 01_to_matrix.R
│   └─ 02_deg.R
├─ 02_rma_microarray/      # STEP2  CEL→RMA (GSE96804 / GSE30529 / GSE104948·104954)
│   ├─ 01_gse96804.R
│   ├─ 02_gse30529.R
│   └─ 03_gse104948_104954.R
├─ 03_merge_combat/        # STEP3  병합 + ComBat (train=사구체 / valid=세뇨관)
│   └─ 01_merge_combat.R
├─ 04_ml_roc/              # STEP4  LASSO∩SVM-RFE → ROC (논문 본문 설계 확정본)
│   ├─ 01_paper_design_lasso_svm_roc.R
│   ├─ 02_roc_plots.R
│   └─ _archive/           #   옛 탐색본(A/B 후보목록·SVM-RFE 보정) 보존
│       ├─ 01_candidate_lists.R
│       ├─ 02_ml_explore_AB.R
│       └─ 03_svmrfe_fix.R
├─ 05_gsea/                # STEP5  Hallmark GSEA + KEGG ssGSEA
│   └─ 01_gsea.R
├─ 06_mr/                  # STEP6  2-표본 MR (FinnGen → GWAS Catalog)
│   └─ 01_mr.R
├─ 07_scrna/               # STEP7  scRNA Seurat (GSE209781) → FN1/ALDH2 세포유형별 발현
│   └─ 01_scrna.R
├─ results/                # 분석 산출물 (config RES_DIR)
├─ _reference_original/    # 원본 저자 스크립트 복사본(G:\ 하드코딩) — 참고용, 파이프라인 제외
├─ CLAUDE.md  DECISIONS.md  RESTRUCTURE_REPORT.md
```
중간 매트릭스(전처리 결과)는 `../data/processed/` (config `OUT_DIR`) 에 쓴다. 원본 데이터 `../data/` 번호 폴더는 읽기 전용.

## 실행법
```r
# (1) 개별 STEP — 어느 위치에서 Rscript 를 돌려도 here 가 루트를 찾음
Rscript 01_deg_gse142025/01_to_matrix.R
Rscript 01_deg_gse142025/02_deg.R
...
# (2) 전체 순차 (무거운 STEP 포함 — 시간·메모리 주의)
Rscript run_all.R
```
- 모든 스크립트는 첫 부분에서 `source(here::here("config.R"))` 로 경로를 참조한다(하드코딩 없음).
- 경로가 다르면 **`config.R` 의 `DATA_ROOT` 한 줄만** 수정.

## 옛 파일 → 새 위치 매핑
| 옛 (평면 구조) | 새 (STEP 폴더) |
|---|---|
| 00_config_local.R | **config.R** (here 기반으로 재작성) |
| 01_GSE142025_to_matrix.R | 01_deg_gse142025/01_to_matrix.R |
| 02_GSE142025_DEG.R | 01_deg_gse142025/02_deg.R |
| 03_GSE96804_to_matrix.R | 02_rma_microarray/01_gse96804.R |
| 04_GSE30529_to_matrix.R | 02_rma_microarray/02_gse30529.R |
| 05_GSE104948_104954_to_matrix.R | 02_rma_microarray/03_gse104948_104954.R |
| 06_merge_combat.R | 03_merge_combat/01_merge_combat.R |
| 10_step34_paper_design.R | 04_ml_roc/01_paper_design_lasso_svm_roc.R |
| 12_step4_roc_plots.R | 04_ml_roc/02_roc_plots.R |
| 07_step4_candidate_lists.R | 04_ml_roc/_archive/01_candidate_lists.R |
| 08_step4_ml.R | 04_ml_roc/_archive/02_ml_explore_AB.R |
| 09_step4_svmrfe_fix.R | 04_ml_roc/_archive/03_svmrfe_fix.R |
| 11_step5_gsea.R | 05_gsea/01_gsea.R |
| 13_step6_mr.R | 06_mr/01_mr.R |
| 14_step7_scrna.R | 07_scrna/01_scrna.R |
| data preprocessing/ML/MR/single-cell *.R (원본 복사본) | _reference_original/ (파이프라인 제외) |

## 주요 재현 결과 (요약; 상세는 DECISIONS.md)
- STEP4: FN1 단일 ROC train 0.909 / valid 0.915 (논문 0.911/0.911). LASSO 6개 중 5/6 일치.
- STEP5: FN1 = EMT/ECM, ALDH2 = 산화적 인산화/대사 (유의).
- STEP6: MR IVW — FN1 위험(OR 2.78), ALDH2 보호(OR 0.67), GWAS Catalog. 저자 Supp9 와 정확 일치.
- STEP7: FN1 = 내피/系膜세포, ALDH2 = 근위세뇨관(PCT).

## 백업
정리 전 원본 스냅샷: `../R-scripts-working_backup_20260708/` (원상복구 보험).
