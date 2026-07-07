# PORT_REPORT.md — 숙제3 (Python 이식) 결과

일시: 2026-07-08 (자율 실행). 대상: `code/Python-scripts-working/` (신규).
기준: `../R-scripts-working/` (숙제2 구조화본)을 STEP 1~7 그대로 이식.

## 1. 구조 (R 과 대칭)
```
config.py  run_all.py  requirements.txt  README.md  PORT_REPORT.md
01_deg_gse142025/{01_to_matrix,02_deg}.py
02_rma_microarray/{_rma_note,01_gse96804,02_gse30529,03_gse104948_104954}.py
03_merge_combat/01_merge_combat.py
04_ml_roc/{01_paper_design_lasso_svm_roc,02_roc_plots}.py
05_gsea/01_gsea.py   06_mr/01_mr.py   07_scrna/01_scrna.py
results/
```
- `config.py`: `CODE_ROOT = Path(__file__).parent`, `DATA_ROOT = CODE_ROOT.parent/"data"` (pathlib 자동탐지).
  R config.R 와 **1:1 동일 값**(LOGFC 0.585, ADJP 0.05, 그룹 라벨, DIR_GSE*). 하위 STEP 에서 `import config`.
- 입력 매트릭스는 `../data/processed/` (R 산출물) 공유. Python 결과는 `results/` (R 과 분리, 비파괴).

## 2. R ↔ Python STEP·라이브러리 매핑
| STEP | R | Python | 이식 방식 |
|---|---|---|---|
| 1 DEG | limma lmFit/eBayes | scipy Welch t + BH | eBayes(moderated) 근사 — 개수 소폭 차이 가능 |
| 2 RMA | affy/oligo/makecdfenv | **없음** | RMA 는 R 전용 → R 산출 매트릭스 재사용(로더+TODO) |
| 3 ComBat | sva::ComBat | inmoose.pycombat_norm | batch=데이터셋, 교집합 병합 동일 |
| 4 LASSO | glmnet cv.glmnet(L1) | sklearn LogisticRegression(L1) | λ 선택 방식 차이(경계 유전자 다를 수 있음) |
| 4 SVM-RFE | caret rfe + svmRadial | sklearn RFECV + 선형 SVC | 커널/선택기준 차이 |
| 4 ROC | pROC(auto-direction) | sklearn roc_auc_score + `max(a,1-a)` | 자동 방향 보정 반영 |
| 4 결합모델 | RF/e1071/xgboost | sklearn RF·SVC / xgboost | 동일 4모델(XGBoost 없으면 GBM 대체) |
| 5 GSEA/ssGSEA | clusterProfiler / GSVA | gseapy.prerank / gseapy.ssgsea | logFC 랭킹·gmt 동일 |
| 6 MR IVW | TwoSampleMR | numpy 직접 IVW | 저자 instrument(Supp8)·harmonise·IVW 공식 동일 |
| 7 scRNA | Seurat + harmony | scanpy + harmonypy + leidenalg | QC/정규화/통합/클러스터/주석 동일 절차 |

## 3. 검증 (비파괴)
- **py_compile: 0 에러 / 14 파일** (config, run_all, 전 STEP).
- **하드코딩 경로 잔존: 0건** (grep `G:\`, `C:/Users`, `C:\Users`, `AppData`).
- **config 스모크**: 루트·하위폴더(`06_mr/`) 모두 `import config` 성공, DATA_ROOT/OUT_DIR/공유 매트릭스/Supp 폴더 실존 True.
- **읽기전용 파리티 스모크** (GSE96804, pandas+numpy 만; 매트릭스 비덮어씀):

| 지표 | Python | R | 판정 |
|---|---|---|---|
| FN1 train AUC | **0.909** | 0.909 | **정확 일치** |
| ALDH2 train AUC | 0.940 (=1−0.060, 자동방향) | 0.940 | 일치(방향보정 후) |

  → 동일 매트릭스·동일 판별력 확인. (ALDH2 는 보호유전자라 pROC 자동방향 관례 필요 — STEP4 에 `max(a,1-a)` 반영.)
- 무거운/대용량 실행은 하지 않음(규칙4): STEP2 RMA·STEP6 FinnGen 2.1GB·STEP7 scRNA 미실행.

## 4. 설치·실행 상태 (현 환경)
- 설치됨: **pandas, numpy** → STEP1(수정 없이)·파리티 스모크 실행 가능.
- 미설치(스텁+TODO 처리): **scipy, scikit-learn, matplotlib, xgboost, inmoose(pycombat), gseapy, scanpy, harmonypy, leidenalg**.
  → 각 STEP 은 미설치 시 `TODO: 설치 후 실행` 안내 후 스킵하도록 방어 코딩. `requirements.txt` 로 일괄 설치.

## 5. R↔Python 대조 (실행 가능 시)
| STEP | 대조 지표 | 상태 |
|---|---|---|
| 1 | DEG 개수(3대비) | 코드 완성. scipy 설치 후 실행 → R DEG 수와 대조 (eBayes 근사 주의) |
| 4 | FN1/ALDH2 단일 ROC-AUC | **파리티 확인**: FN1 0.909 일치 (train). 전체는 sklearn 설치 후. |
| 4 | LASSO∩SVM-RFE 교집합 | sklearn 설치 후 (R: FN1/ALDH2 포함 6개) |
| 6 | MR IVW OR/p | 코드 완성. FinnGen 로드는 무거움 → 실행 시 R(FN1 OR2.78/ALDH2 OR0.67)과 대조 |
| 5,7 | GSEA/scRNA | gseapy/scanpy 설치 후 (무거움) |

## 6. 미완/한계 (TODO)
- **STEP2 RMA**: Python 미대응(R 선행 필수). 로더만 제공.
- **limma eBayes**: Welch t 근사 → DEG 개수 R 과 완전 동일 보장 못함.
- **MR LD clumping**: 오프라인 불가 → 저자 도구변수(Supp8) 사용(R 과 동일 방침).
- **미설치 라이브러리**: 실제 수치 재현은 `requirements.txt` 설치 후 각 STEP 실행 필요.
- 기존 `../../pipeline_port/` 는 참고만 하고 재사용/수정하지 않음(신규 작성).
