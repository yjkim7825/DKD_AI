# Python-scripts-working — P25 파이프라인 Python 이식본

R 파이프라인(`../R-scripts-working/`)을 STEP 1~7 그대로 Python 으로 이식.
임계값·대비방향·그룹 라벨 규약을 R 과 1:1 일치시킴(LOGFC 0.585, adjP 0.05, ComBat batch=데이터셋,
LASSO∩SVM-RFE, IVW 등). 경로는 `config.py` 한 곳(pathlib 자동탐지).

## 폴더 구조
```
Python-scripts-working/
├─ config.py               # ★ 경로/파라미터 (pathlib 루트 자동탐지). R config.R 와 1:1 값.
├─ run_all.py              # 01~07 순차 실행 마스터
├─ requirements.txt
├─ 01_deg_gse142025/{01_to_matrix,02_deg}.py
├─ 02_rma_microarray/{_rma_note,01_gse96804,02_gse30529,03_gse104948_104954}.py
├─ 03_merge_combat/01_merge_combat.py
├─ 04_ml_roc/{01_paper_design_lasso_svm_roc,02_roc_plots}.py
├─ 05_gsea/01_gsea.py   06_mr/01_mr.py   07_scrna/01_scrna.py
└─ results/               # Python 분석 결과 (R results 와 분리)
```
- 입력 매트릭스(RMA 산출물)는 `../data/processed/` (R 파이프라인과 공유) 에서 읽는다.
- 원본 데이터 `../data/` 번호 폴더는 읽기 전용.

## 실행법
```bash
py -m pip install -r requirements.txt
py 01_deg_gse142025/01_to_matrix.py     # 개별 STEP (config 자동탐지)
py run_all.py                            # 전체 순차 (무거운 STEP 포함)
```
- 모든 스크립트는 `sys.path` 에 루트를 넣고 `import config` 로 경로 참조(하드코딩 없음).
- 데이터 위치가 다르면 `config.py` 의 `DATA_ROOT` 한 줄만 수정.

## R ↔ Python 라이브러리 매핑
| STEP | R | Python |
|---|---|---|
| 1 DEG | limma (lmFit/eBayes) | scipy Welch t + BH (eBayes 근사) |
| 2 RMA | affy/oligo | **없음** — R 산출 매트릭스 재사용 (RMA 는 R 전용) |
| 3 ComBat | sva::ComBat | inmoose.pycombat_norm |
| 4 LASSO/SVM-RFE/ROC | glmnet / caret rfe / pROC | sklearn LogisticRegression(L1) / RFECV+SVC / roc_auc_score |
| 4 결합모델 | RF/e1071/xgboost | sklearn RF·SVC / xgboost |
| 5 GSEA/ssGSEA | clusterProfiler / GSVA | gseapy.prerank / gseapy.ssgsea |
| 6 MR IVW | TwoSampleMR | numpy 직접 구현 |
| 7 scRNA | Seurat + harmony | scanpy + harmonypy + leidenalg |

## 주의 (PORT_REPORT.md 상세)
- STEP2 RMA 는 Python 미대응 → R 선행 필요(입력 매트릭스 재사용).
- limma eBayes, MR LD clumping 등 R 특화 부분은 근사/저자 instrument 로 대체(주석 명시).
- 미설치 라이브러리(pycombat/gseapy/scanpy 등)는 함수+안내로 스텁 처리('TODO: 설치 후 실행').
