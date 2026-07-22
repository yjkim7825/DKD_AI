# 05_ml — 머신러닝 진단모델 (FN1·ALDH2)

DEG∩MR로 좁힌 후보 유전자 10개에서 **진단 바이오마커**를 선별·검증한다.

## 두 방식을 다 구현해보고 논문 방식을 채택했다

논문 본문/그림과 저자가 올린 코드가 **선별 방식이 서로 달라서**, 둘 다 구현해 같은 입력(논문 코호트)으로 돌려 비교했다.

| | 방식 A — 논문 Figure 4 | 방식 B — 저자 업로드 코드 |
| --- | --- | --- |
| 스크립트 | `01_ml_feature.R` + `02_ml_model.R` | `reference/ml_author.R` |
| 선별 | LASSO → **ROC-AUC>0.8 필터** → 4모델 결합 | LASSO **∩ SVM-RFE** 교집합 |
| 최종 | **FN1 · ALDH2** | VNN2 1개로 붕괴 |
| 결과 | GLM 검증 AUC 0.972 (논문 0.942 재현·상회) | 덜 좁혀짐(SVM-RFE 구현 민감) |

**결론: 논문 방식(A)이 검증셋에서 무너지는 유전자를 ROC 이중필터로 정확히 걸러내 더 깔끔하고 재현성도 높았다(실행 결과 B는 VNN2 1개로 붕괴, A는 FN1·ALDH2 복원). 그래서 최종적으로 논문 방식(`01_ml_feature.R`+`02_ml_model.R`)을 채택했다.** 방식 B와 저자 원본은 `reference/`로 옮겨 대조용으로만 보관한다.

자세한 수치 비교: `../docs/05_ML_두방식비교.md` · 개념도: `../docs/05_ML_개념도.svg`

## 실행
```r
# 학습 GSE96804 / 검증 GSE104948+104954 (data/processed/*.paper.txt)
source("run_ml.R")                 # 방식 A 실행 (+ Figure 4 D~I 그림)
# source("reference/ml_author.R")  # 방식 B(참고) — caret·kernlab·VennDiagram 필요
```
필요 패키지: glmnet, pROC, sva, randomForest, e1071, xgboost (방식 B 추가: kernlab, caret, VennDiagram)

## 파일
| 경로 | 내용 |
| --- | --- |
| `01_ml_feature.R` | **채택 1/2** ComBat→LASSO→ROC>0.8 (유전자 선별) → `A_final_genes.txt` |
| `02_ml_model.R` | **채택 2/2** 4모델 결합 + Figure 4 D~I 그림 |
| `run_ml.R` | 01→02 한 번에 실행 |
| `R/ml_func.R` | 공용 함수(로딩·ComBat·표준화·AUC) |
| `interGenes.List.txt` | 후보 유전자 10개 (DEG∩MR) |
| `ml_compare.py` | 두 방식 numpy 비교(검증용) |
| `reference/machine learning modeling 1·2.R` | 저자 업로드 **원본** 코드 |
| `reference/ml_author.R` | 저자 방식 재현 (LASSO∩SVM-RFE + 유전자별 ROC) |
| `reference/combat_mod_compare.R` | ComBat mod=NULL vs ~그룹 비교 |
| `output/A_PCA_batch.pdf` · `A_Figure4_DI.pdf` | Fig 4A PCA · Fig 4 D~I |
| `output/` | 선택 유전자·ROC·모델 AUC·그림 |
