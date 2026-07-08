# RESTRUCTURE_REPORT.md — R 구조화 결과

일시: 2026-07-08 (자율 실행). 대상: `code/R-scripts-working/`.

## 1. 무엇을 했나
평면 구조(00~14 번호 스크립트 + 원본 복사본 혼재)를 **STEP 번호 정렬 폴더 구조**로 재편.
경로를 `config.R` 한 곳으로 집중하고 `here` 기반 루트 자동탐지로 하드코딩을 제거.

## 2. 안전 조치
- **백업**: `../R-scripts-working_backup_20260708/` 에 정리 전 스냅샷(106 파일, `.git`/`.claude` 제외) — robocopy.
- **원본 무수정**: `../R-scripts-for-pipeline-reproducibility/` 와 `../data/` 는 읽기만. 이 작업은 working 복제본에서만.
- working 안에 있던 **원본 저자 스크립트 복사본**(G:\ 하드코딩본)은 삭제하지 않고 `_reference_original/` 로 이동 보존.

## 3. 새 구조
```
config.R  run_all.R  .here
01_deg_gse142025/{01_to_matrix,02_deg}.R
02_rma_microarray/{01_gse96804,02_gse30529,03_gse104948_104954}.R
03_merge_combat/01_merge_combat.R
04_ml_roc/{01_paper_design_lasso_svm_roc,02_roc_plots}.R  + _archive/{01_candidate_lists,02_ml_explore_AB,03_svmrfe_fix}.R
05_gsea/01_gsea.R   06_mr/01_mr.R   07_scrna/01_scrna.R
results/   _reference_original/   README.md
```

## 4. 옛 → 새 매핑 (이동·리네임)
| 옛 | 새 |
|---|---|
| 00_config_local.R | config.R (here 기반 재작성, 옛 파일 삭제·백업보존) |
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
| data preprocessing 1~3 / differential expression analysis / GSEA / machine learning 1~2 / Mendelian 1~3 / single-cell 1~4 .R | _reference_original/ (14개, 파이프라인 제외) |
| Supplementary Table 3~7 *.csv, VennDiagram*.log | _reference_original/ (7개) |

## 5. 경로/하드코딩 처리
- 모든 스크립트의 `source("00_config_local.R")` → `library(here); source(here::here("config.R"))` 로 교체.
- 하드코딩 scratch 절대경로(05/13/14) → config 의 `SCRATCH_DIR`(= `tempdir()/pipeline_scratch`) 참조로 교체.
- `config.R`: `CODE_ROOT <- here::here()`, `DATA_ROOT <- normalizePath(file.path(CODE_ROOT,"..","data"))`.
  → 사용자는 데이터 위치가 다를 때 `DATA_ROOT` 한 줄만 수정.
- `.here` 앵커 추가(+ `.git` 존재) → 스크립트를 하위 폴더에서 실행해도 루트 정확 탐지.

## 6. 검증 (비파괴)
- **parse 문법검사: 0 에러 / 16 스크립트** (config, run_all, 전 STEP + _archive 포함).
- **하드코딩 경로 잔존: 0건** (grep `G:\`, `C:/Users`, `C:\Users`, `/scratchpad/`, `00_config_local` — _reference_original 제외).
- **config 경로 실존: 8/8 TRUE** (DATA_ROOT, OUT_DIR, RES_DIR, DIR_GSE96804/104948/104954/142025/30529).
- **하위폴더 스모크**: wd=`07_scrna/` 에서 `here::here()` → 루트 정확, `config.R`·타 STEP 스크립트 해석 성공.
- 무거운 전체 재실행은 하지 않음(규칙4). 문법·경로 정합만 확인.

## 7. 이슈/유의
- `_reference_original/` 의 원본 복사본은 **의도적으로 G:\ 하드코딩 유지**(저자 원본 대조용). 파이프라인·검증에서 제외.
- `SCRATCH_DIR` 를 `tempdir()` 기반으로 바꿔, STEP6/7 의 기존 캐시(bgzip/seurat rds)는 새 세션에서 재생성됨(무거운 재실행은 사용자 판단).
- 산출물 경로(`../data/processed`, `results/`)·임계값(LOGFC 0.585, ADJP 0.05)·그룹 라벨 규약은 그대로 계승.
- CLAUDE.md 는 옛 파일명을 언급(과거 기록) — 최신 구조는 본 리포트/README 기준.
