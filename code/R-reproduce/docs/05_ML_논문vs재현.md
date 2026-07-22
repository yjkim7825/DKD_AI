# 05. ML — 논문 vs 우리 재현

## 📂 파일 경로
- **코드**: `R-reproduce/05_ml/` — 채택: `01_ml_feature.R`(선별) + `02_ml_model.R`(모델링) · `run_ml.R`(한번에) · `reference/`(저자 원본·대조)
- **결과**: `05_ml/output/` (A_lasso.gene.txt, A_roc_filter.csv, A_model_AUC.csv, A_PCA_batch.pdf, A_cvfit.pdf, A_Figure4_DI.pdf)
- **개념도·비교**: `05_ML_개념도.svg`, `05_ML_두방식비교.md`

---

## 0. 한눈에 — 논문 Figure 4 재현 결과

| 항목 | 논문 | 우리 재현 (`01_ml_feature.R`) | 판정 |
| --- | --- | --- | --- |
| 후보 유전자 | 10개 (DEG∩MR) | 동일 10개 | ✅ |
| 배치보정 | sva (ComBat) | **sva ComBat** (batch=학습/검증, mod=NULL) | ✅ 동일 도구 |
| LASSO 선택 | 6개 (Fig 4B) | 6개 (5개 일치) | ✅ 근접 |
| ROC>0.8 최종 | **FN1 · ALDH2** | **FN1 · ALDH2** | ✅ 완전일치 |
| 결합 최고모델 | GLM (검증 0.942) | **GLM (검증 0.972)** | ✅ 결론·수치 재현 |
| 방향 | FN1 위험↑ / ALDH2 보호↓ | 동일 | ✅ |

→ **최종 진단 바이오마커 FN1·ALDH2, "GLM 최고" 결론까지 정확히 재현**. 검증 AUC는 논문(0.942)을 상회.

---

## 1. 사용 데이터 (논문 코호트 그대로)
- **학습(training)**: `GSE96804` — 61명 (Control 20 / DKD 41)
- **검증(validation)**: `GSE104948 + GSE104954` — 71명 (Control 42 / DKD 29)
- 파일: `data/processed/data.train.paper.txt`, `data.valid.paper.txt`
- ※ 전처리(prep3)는 검증 내부(104948↔104954)만 ComBat 병합 → **학습↔검증은 미보정**이라 05에서 배치보정(아래 2단계)

---

## 2. 단계별 재현 결과

### ① 배치보정 (논문 Methods: "데이터셋 간 배치효과 sva로 보정")
- `sva::ComBat`, batch = 학습 vs 검증, mod = NULL (라벨 미사용, 무누수)
- Fig 4A(배치보정 후 PCA)도 재현 → `A_PCA_batch.pdf`

### ② LASSO (10겹 CV, λmin) — Fig 4B
| | 선택 유전자 |
| --- | --- |
| **논문 Fig4B (6개)** | CDKN1B, ALDH2, FN1, XAF1, TSPYL5, VNN2 |
| **우리 (6개)** | ALDH2, CDKN1B, FN1, IFI44L, VNN2, XAF1 |
- **5개 일치** — 우리는 TSPYL5 대신 IFI44L (정규화·probe 매핑 경계 차이, 이전 재현차이와 동일 원인)

### ③ ROC-AUC > 0.8 필터 (학습·검증 둘 다) — Fig 4C~G
| 유전자 | 학습 AUC | 검증 AUC | 통과 |
| --- | --- | --- | --- |
| **ALDH2** | 0.940 | 0.807 | ✅ |
| **FN1** | 0.909 | 0.871 | ✅ |
| CDKN1B | 0.967 | 0.718 | ✗ |
| VNN2 | 0.974 | 0.649 | ✗ |
| IFI44L | 0.627 | 0.640 | ✗ |
| XAF1 | 0.778 | 0.630 | ✗ |

→ **최종 FN1·ALDH2** — 논문과 정확 일치. 논문도 "검증에서 나머지 4개는 AUC<0.7, FN1·ALDH2만 >0.8"이라 보고(Fig 4C–G).

### ④ 4모델 결합 (FN1+ALDH2) — Fig 4H·I : 논문 vs 재현 (학습 / 검증 AUC)
| 모델 | 논문 (학습/검증) | 재현 (학습/검증) |
| --- | --- | --- |
| **GLM** | 0.978 / **0.942** | 0.980 / **0.972** |
| RF | 1.000 / 0.815 | 1.000 / 0.845 |
| SVM(RBF) | 0.977 / 0.807 | 0.973 / 0.875 |
| XGBoost | 0.999 / 0.844 | 0.997 / 0.825 |

→ **학습 AUC 4모델 다 거의 동일**, **GLM이 검증 최고**(논문 0.942 / 재현 0.972) — "GLM 우수" 결론 재현. RF·SVM·XGB 검증값은 ±0.03~0.07 차이(배치보정·시드·RF 랜덤성). 네 모델 모두 검증 AUC>0.8.

---

## 3. 차이 요약

| 항목 | 논문 | 우리 | 원인·판정 |
| --- | --- | --- | --- |
| LASSO 6개 중 1개 | TSPYL5 | IFI44L | probe 매핑·정규화 경계 (결론 무관) |
| 배치보정 배치정의 | 데이터셋별 | 학습/검증 2배치 | 검증은 이미 내부 병합됨 → 재분할 안 함 |
| GLM 검증 AUC | 0.942 | 0.972 | 재현이 논문 상회 (방법 동일) |
| 최종 유전자 | FN1·ALDH2 | **동일** | ✅ |

- 어느 차이도 **결론(FN1 위험·ALDH2 보호, GLM 최고)** 을 바꾸지 않음.
- **저자 코드 방식(LASSO∩SVM-RFE)** 은 별도 실행 시 VNN2 1개로 붕괴 → 논문 방식(ROC 필터)이 옳음을 실증 (상세 `05_ML_두방식비교.md`).

### ※ 저자 업로드 코드(ML1/2.R) 원본 대조 확인
업로드된 `machine learning modeling 1.R`·`2.R`을 우리 재현본 **`reference/ml_author.R`**과 직접 대조 — 로직 완전 일치:

| 원본 단계 | ml_author.R | 판정 |
| --- | --- | --- |
| LASSO (glmnet, λmin, 계수≠0) | B1 | ✅ 동일 |
| SVM-RFE (caret rfe·svmRadial·sizes 2~8) | B2 | ✅ 동일 |
| Venn 교집합 (Reduce intersect) | B3 | ✅ 동일 |
| 유전자별 ROC 학습+검증 (ci.auc bootstrap) | B4 | ✅ 동일 (검증 ROC도 생성하도록 보완) |

- ※ `01_ml_feature.R`은 ML1의 재현이 아니라 **논문 본문 방식**(의도�