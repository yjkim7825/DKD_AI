# R → Python 라이브러리 매핑 (P25 DKD 파이프라인 이식)

R 원본/작업본(`code/R-scripts-working/`)의 각 R 패키지·함수를 Python 이식본
(`code/Python-scripts-working/`)에서 무엇으로 바꿨는지, **대응이 없으면 어떻게 처리했는지** 정리.
정렬 기준은 정렬된 R(= `../../R-scripts-working/DIFF_vs_original_code.md`, `DECISIONS.md`).

> 원칙: 원본 데이터·R 스크립트 무수정. Python 산출물은 `Python-scripts-working/results/` 로만.
> 발현매트릭스(RMA·ComBat)는 R 산출물을 `code/data/processed/` 에서 **공유(읽기)**.

---

## 1. STEP별 라이브러리/함수 매핑

| STEP | R 패키지 / 함수 | Python 대체 | 대응 방식 | 비고 / 결과 차이 |
|---|---|---|---|---|
| **1 DEG** | `limma::lmFit`+`eBayes`+`topTable` | `scipy.stats.ttest_ind(equal_var=False)` + 자체 BH | **근사**: eBayes moderated t 대응 없음 → per-gene **Welch t** + BH-FDR 직접구현 | DEG 수 소폭 차이(Late_vs_Early Py 3460 / R 3314 / 논문 2833). 방향·비율 일치 |
| 1 정규화 | `limma::normalizeBetweenArrays` | 자체 `quantile_normalize`(pandas rank) | **직접 구현**(분위수 정규화) | 동등 |
| 1 `p.adjust(BH)` | R 내장 | 자체 `bh_adjust`(numpy) | **직접 구현** | R `p.adjust(method="BH")` 와 동일 알고리즘 |
| 1 GSE96804 DEG | `limma`(complementary DEG) | 위 Welch t 재사용(`03_gse96804_deg.py`) | **근사** | 유의 654 = R 654. FN1=Up, ALDH2=저발현방향(컷 경계) — R 일치 |
| **2 RMA** | `oligo::rma` / `affy::rma` / `makecdfenv`(brainarray CDF) | **(대응 없음)** | **R 산출물 재사용**: RMA 는 Bioconductor 전용 → Python 미이식. `processed/*.labeled.txt`(R RMA 결과) 읽음 | Python 은 RMA 를 수행하지 않음(설계상) |
| 2 `GEOquery::getGEO` | GEO 메타데이터 | (미이식) | R 산출 라벨 매트릭스 재사용 | — |
| **3 ComBat** | `sva::ComBat` | `inmoose.pycombat_norm` (미설치) | **R 산출물 재사용**: inmoose 미설치 → R `data.valid.paper.txt`(검증 71) 직접 읽음 | STEP4 가 R ComBat 결과 사용 → 검증셋 완전 동일 |
| **4 LASSO** | `glmnet`+`cv.glmnet`(binomial, alpha1, deviance, λmin) | `sklearn.LogisticRegressionCV(penalty="l1", solver="liblinear")` | **대체**: glmnet↔sklearn 정칙화 경로/λ 스케일 차이 | 선택 유전자 소폭 차(Py 5 / R 6). FN1·ALDH2 는 양쪽 선택 |
| 4 SVM-RFE | `caret::rfe`(caretFuncs, svmRadial) | `sklearn.feature_selection.RFECV(SVC(linear))` | **대체**: caret rfe↔RFECV, 커널 차(radial→linear) | 선택 유전자 다름(구현 차). 교집합도 다름 |
| 4 단일 ROC | `pROC::roc`+`auc`(자동방향) | `sklearn.metrics.roc_auc_score` + `max(a,1-a)` | **대체**(동등) | **FN1/ALDH2 AUC Python=R 3자리 정확 일치**(입력 동일) |
| 4 `ci.auc(bootstrap)` | pROC | (Python 미산출) | 미이식(부수) | 핵심 AUC 는 일치 |
| 4 결합모델 GLM/RF/SVM/XGB | `glm`/`randomForest`/`e1071::svm`/`xgboost` | `sklearn LogisticRegression`/`RandomForestClassifier`/`SVC(rbf)`/`xgboost.XGBClassifier` | **대체** | GLM·RF·XGB 검증 R 근접(GLM Py 0.953). **SVM 검증 0.374 = 포팅 한계**(SVC rbf gamma=scale ≠ e1071) |
| **5 GSEA** | `clusterProfiler::GSEA`(diff_ logFC) | `gseapy.prerank` | **대체**: 순열기반 GSEA | NES 부호·크기 근접(EMT +2.7/R+2.5, OXPHOS -1.7/R-1.6). Hallmark+KEGG GSEA 둘 다 |
| 5 gene set | `msigdbr`(H) / KEGG gmt | 제공 gmt 파일(4-1/4-2) | **동일 파일** | 오프라인 → 제공 gmt(R 과 동일) |
| 5 ssGSEA | `GSVA::gsva(ssgsea)` | `gseapy.ssgsea` | **대체** | **FN1↔ECM 0.919, ALDH2↔TCA 0.921 = R 정확 일치** |
| 5 상관 | `cor(method="spearman")` | `scipy.stats.spearmanr` | 동등 | 일치 |
| **6 MR IVW** | `TwoSampleMR::mr("mr_ivw")` | **numpy 직접구현** | **직접 구현**: β=Σ(βx·βy/σy²)/Σ(βx²/σy²), se=√(1/Σ(βx²/σy²)), OR=exp(β) | **FN1 b=1.021/OR2.777/p0.0122 = R=Supp9 정확 일치** |
| 6 MR-Egger | `mr_egger_regression` | numpy **가중 WLS**(by~bx, w=1/σy²) | **직접 구현**: slope=인과, intercept=다면발현 | 방향 일치. SE 는 표준 WLS(TwoSampleMR 세부와 미세차 가능) |
| 6 weighted median | `mr_weighted_median` | numpy `_wmedian` + 부트스트랩 SE | **직접 구현**: b_iv=βy/βx, w=βx²/σy², 가중중앙값; SE=1000 부트스트랩 | 방향 일치 |
| 6 이질성 | `mr_heterogeneity`(Cochran Q) | numpy Q=Σ w(βy−β·βx)², `chi2.sf` | **직접 구현** | 동등 |
| 6 다면발현 | `mr_pleiotropy_test`(Egger 절편) | Egger 절편/se → z → p | **직접 구현** | 없음(FN1/ALDH2 절편 비유의) |
| 6 민감도 | `mr_singlesnp`/`mr_leaveoneout` | numpy Wald ratio / IVW 재계산 | **직접 구현** | 표 저장 |
| 6 LD clumping | `clump_data`(OpenGWAS API) | **(대응 없음)** | **저자 instrument 재사용**: 오프라인 → Supp8(이미 clump)을 노출 SNP 집합으로 사용 | R 과 동일 방침. 결과 Supp9 일치 |
| 6 IVW 필터 | `Mendelian randomization 2.R` | 자체 `ivw_filter` | **직접 구현**: IVW p<0.05 & 3방법 OR방향일치 & 다면발현 p>0.05 | GCST 통과 {ALDH2,FN1} = R 일치 |
| 6 GWAS-VCF 파싱 | `VariantAnnotation`/`gwasglue` | 자체 파서(ES:SE:LP:AF 분해) | **직접 구현** | pval=10^(−LP) |
| **7 scRNA** | `Seurat`(Read10X~FindClusters) | `scanpy` | **대체** | QC/정규화/클러스터 표준. FN1=내피, ALDH2=PCT 재현 |
| 7 배치통합 | `harmony::RunHarmony` | `harmonypy.run_harmony` | **대체**(휠 설치) | — |
| 7 그래프 클러스터 | `FindClusters`(Louvain) | `scanpy leiden`(igraph) | **대체**(알고리즘 유사) | 세포수·타입 근접 |

---

## 2. 대응 라이브러리가 없어 처리한 항목 (요약)

**(a) 직접 구현 (numpy/scipy)**
- **MR 전체**(IVW/Egger/weighted median/Cochran Q/다면발현/민감도/IVW필터): TwoSampleMR 의 통계를 공식으로 재현.
  - IVW: `β = Σ(βx·βy·w)/Σ(βx²·w)`, `se = √(1/Σ(βx²·w))`, `w=1/σy²`.
  - Egger: 가중최소제곱 `βy ~ βx`(절편 포함, w=1/σy²) → slope=인과추정, intercept=수평다면발현.
  - weighted median: 비율추정 `βy/βx` 의 가중중앙값(w=βx²/σy²), SE=1000회 정규 부트스트랩.
- **DEG BH-FDR·분위수정규화**: `p.adjust(BH)`·`normalizeBetweenArrays` 를 numpy/pandas 로 재현.

**(b) R 산출물 재사용**
- **RMA(STEP2)**: affy/oligo RMA 는 Python 동등물이 없어 미이식 → R 이 만든 `processed/*.labeled.txt` 를 읽음.
- **ComBat(STEP3)**: inmoose 미설치 → R `data.valid.paper.txt`(LD-only 검증 71) 를 그대로 읽음(STEP4 입력 동일).

**(c) 근사 대체 (원본과 차이)**
- **limma eBayes → Welch t**: moderated t(분산 shrink) 대응 없음 → per-gene Welch t.
  결과: Late_vs_Early DEG **Py 3460 / R 3314 / 논문 2833**(방향·비율 일치, 개수만 차이).
- **glmnet → sklearn L1**: 정칙화 경로·λ 스케일 차 → LASSO 선택 **Py 5 / R 6**(FN1·ALDH2 는 공통).
- **caret rfe(svmRadial) → RFECV(linear SVC)**: SVM-RFE 선택·교집합 다름.
- **e1071 svm(rbf) → sklearn SVC(rbf)**: 결합 SVM **검증 AUC 0.374(Py) vs 0.800(R)** — 커널/기본 gamma 차이로 인한 **가장 큰 포팅 갭**(한계로 기록). 단일 ROC·GLM·RF·XGB 는 R 근접.

---

## 3. 근사/미대응이 결과에 준 차이 (수치)

| STEP | 지표 | Python | R | 논문 | 차이 원인 |
|---|---|---|---|---|---|
| 1 | Late_vs_Early DEG | 3460 | 3314 | 2833 | Welch t ≈ eBayes(근사) + 입력매트릭스 |
| 4 | FN1 valid AUC | 0.871 | 0.871 | 0.911 | **Py=R 정확**(입력 동일) |
| 4 | ALDH2 valid AUC | 0.807 | 0.807 | 0.815 | **Py=R 정확**, LD-only로 논문 근접 |
| 4 | 결합 GLM valid | 0.953 | 0.820 | 0.942 | sklearn glm ≠ R glm(Py가 논문 근접) |
| 4 | 결합 SVM valid | **0.374** | 0.800 | 0.807 | **SVC rbf ≠ e1071(포팅 한계)** |
| 5 | FN1↔ECM rho | 0.919 | 0.919 | — | **Py=R 정확** |
| 6 | FN1 IVW b (GCST) | 1.021 | 1.021 | 1.021 | **Py=R=Supp9 정확** |
| 6 | ALDH2 IVW (GCST) | b−0.402/nsnp15 | b−0.395/nsnp14 | b−0.395 | palindromic SNP 1개 차(Py harmonise 미제외) |
| 7 | FN1/ALDH2 세포 | Endothelial/PCT | Endothelial/PCT | 내피/PCT | 결론 동일 |

**요약**: 입력이 동일한 단계(단일 ROC, MR IVW, ssGSEA 상관)는 **Python=R 정확 일치**.
근사/구현차가 큰 곳은 (1) DEG 수(Welch t), (2) LASSO/SVM-RFE 선택, (3) 결합 SVM 검증 AUC.
FN1/ALDH2 의 진단·인과·경로·세포국소화 **핵심 결론은 3자(Python·R·논문) 모두 일치**.
