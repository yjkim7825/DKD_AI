# DKD Multi-omics 재현 (P25: FN1 & ALDH2)

당뇨병성 신장질환(DKD)의 다중오믹스 바이오마커 논문 **P25**
*(Lin et al., "Multi-omics and machine learning identify FN1 and ALDH2 as diagnostic biomarkers
and therapeutic targets in early and late diabetic kidney disease", Renal Failure 2025)* 의
분석 파이프라인을 **R로 재현**하고 **Python으로 이식**한 재현연구 저장소입니다.
발현 → DEG → 인과(MR) → ML 바이오마커 → 경로(GSEA) → 단일세포로 이어지는 서사를 STEP 1~7로 구현하고,
**논문·R·Python 세 값을 나란히 대조**합니다.

## 재현 하이라이트

핵심 바이오마커 **FN1**(위험·섬유화 축)과 **ALDH2**(보호·산화/대사 축)의 진단 성능이
논문·R·Python에서 일치합니다. Python은 R과 동일한 입력 매트릭스를 공유하여 **소수점까지 같은 값**을 냅니다.

| 지표 | 논문 | R | Python |
|---|---|---|---|
| FN1 단일 ROC-AUC (train / valid) | 0.911 / 0.911 | 0.909 / 0.915 | 0.909 / 0.915 |
| ALDH2 단일 ROC-AUC (train / valid) | 0.912 / 0.815 | 0.940 / 0.784 | 0.940 / 0.784 |
| MR 인과 (GCST90435706) | FN1 위험 / ALDH2 보호 | FN1 OR 2.78 / ALDH2 OR 0.67 | 2.78 / 0.67 |
| GSEA 축 | FN1=EMT/ECM, ALDH2=대사 | rho 0.919 / 0.921 | 0.919 / 0.921 |
| scRNA 국소화 | 내피·메산지움 / 근위세뇨관 | Endothelial / PCT | Endothelial / PCT |

- **FN1은 논문과도 거의 정확 재현**(소수점 둘째 자리), **MR은 저자 Supp9와 완전 일치**.
- 전체 대조: [`pipeline_port/docs/3자_대조표_논문_R_Python.md`](pipeline_port/docs/3자_대조표_논문_R_Python.md)
  · `code/Python-scripts-working/results/python_vs_R_vs_paper.csv`

## 분석 파이프라인 (STEP 1~7)

```
GSE142025 ─▶ [1] DEG (limma) ─────────────────────────────────┐
CEL(96804 등) ─▶ [2] RMA 발현매트릭스 ─▶ [3] 병합 + ComBat ─▶ [4] LASSO + SVM-RFE ─▶ ROC
                                                              │        (FN1 · ALDH2)
eQTL + GWAS ─▶ [6] Mendelian Randomization (인과) ────────────┤
                                                DEG 랭킹 ─▶ [5] GSEA / ssGSEA (경로 축)
scRNA(209781) ─▶ [7] Seurat/scanpy (세포유형 국소화) ─────────┘
```

| STEP | 단계 | 핵심 산출 |
|---|---|---|
| 1 | GSE142025 DEG | Late/Early/Control 대비 차등발현 |
| 2 | CEL → RMA 발현매트릭스 (**R 전용**) | GSE96804 등 정규화 매트릭스 |
| 3 | 병합 + ComBat 배치보정 | train / valid 세트 |
| 4 | LASSO + SVM-RFE → ROC (**핵심**) | FN1·ALDH2 단일·결합 진단모델 |
| 5 | Hallmark GSEA + KEGG ssGSEA | FN1=EMT/ECM, ALDH2=OXPHOS/대사 |
| 6 | 2-표본 MR (IVW/Egger/WM) | FN1 위험·ALDH2 보호 인과 |
| 7 | scRNA (Seurat/scanpy) | 세포유형별 FN1·ALDH2 발현 |

## 폴더 구조

```
DKD_AI/
├── README.md
├── code/
│   ├── R-scripts-working/         # R 재현 파이프라인 (정본)
│   ├── Python-scripts-working/    # R→Python 이식본 (STEP 1:1 대칭)
│   ├── R-scripts-for-pipeline-reproducibility/   # 저자 원본 코드(참조)
│   └── data/
│       ├── processed/             # 축약 매트릭스(train/test/valid) — 공유 입력
│       └── 4-1..gmt / 4-2..gmt    # Hallmark / KEGG 유전자세트
└── pipeline_port/docs/            # 근거표 · 대조표 · 코드흐름 분석
```

### `code/R-scripts-working/` — R 재현 파이프라인
- **`config.R`** — 모든 경로·파라미터 중앙화(`here` 기반 루트 자동탐지). 전 스크립트가 이걸 `source`.
- **`run_all.R`** — STEP 1~7 순차 실행 엔트리포인트.
- **`install_packages.R`** — CRAN + Bioconductor + TwoSampleMR 일괄 설치.
- STEP 폴더:
  - `01_deg_gse142025/` — GSE142025 DEG(limma)
  - `02_rma_microarray/` — CEL → RMA 발현매트릭스
  - `03_merge_combat/` — 데이터셋 병합 + ComBat 배치보정
  - `04_ml_roc/` — LASSO + SVM-RFE → 단일·결합 ROC
  - `05_gsea/` — Hallmark GSEA + KEGG ssGSEA
  - `06_mr/` — Mendelian Randomization
  - `07_scrna/` — Seurat 단일세포
- `results/` — DEG·ROC·MR·GSEA·scRNA 산출물. `_reference_original/` — 저자 원본 코드 참조 보관.

### `code/Python-scripts-working/` — R→Python 이식본 (R과 대칭 구조)
- **`config.py`** — R `config.R`와 **1:1 동일 값**(pathlib 루트 자동탐지).
- **`run_all.py`** — STEP 순차 실행. **`requirements.txt`** — 의존성 일괄 설치.
- STEP 폴더 구성(`01_deg_gse142025` ~ `07_scrna`)이 R과 동일. `results/`에 Python 산출물 + 3자 대조표 CSV.
- ※ **STEP2(RMA)는 R 전용**(`affy/oligo`) — Python은 이식하지 않고 **R 산출 매트릭스를 재사용**합니다.
- venv 폴더 `DKD_AI/`는 git에서 제외됩니다.

### R ↔ Python 대응
두 파이프라인은 **STEP 번호로 1:1 대응**하며, 같은 `code/data/processed/` 매트릭스를 **공유 입력**으로 씁니다.
따라서 동일 입력·동일 판별식에서 결과가 일치합니다(단일 ROC·MR·GSEA·scRNA에서 소수점까지 동일).

| STEP | R (`R-scripts-working/`) | Python (`Python-scripts-working/`) |
|---|---|---|
| 1 DEG | limma lmFit/eBayes | scipy Welch t + BH (근사) |
| 2 RMA | affy/oligo RMA | **없음** → R 산출 매트릭스 재사용 |
| 3 ComBat | sva::ComBat | inmoose(미설치 시 R valid 재사용) |
| 4 ML/ROC | glmnet · caret · pROC | scikit-learn |
| 5 GSEA | clusterProfiler · GSVA | gseapy |
| 6 MR | TwoSampleMR | numpy 직접 IVW |
| 7 scRNA | Seurat + harmony | scanpy + harmonypy |

## 데이터셋

| 데이터셋 | 역할 | 출처 |
|---|---|---|
| GSE96804 | **훈련셋** (61: DKD 41 / Control 20) | GEO |
| GSE104948 / GSE104954 | **독립 검증셋** (사구체 / 세뇨관) | GEO |
| GSE142025 | early/late DKD **DEG 원천** (36) | GEO |
| GSE30529 | DEG 검증 코호트 (22) | GEO |
| GSE131882 / GSE209781 / GSE266146 | 단일세포(scRNA) 3종 | GEO |
| FN1 · ALDH2 cis-eQTL | MR **노출** | OpenGWAS (eqtl-a) |
| FinnGen R12 / GCST90435706 | MR **결과**(DKD outcome) | FinnGen / GWAS Catalog |
| Hallmark / KEGG gmt | GSEA 유전자세트 | MSigDB |

> **원본 RAW/GWAS(약 12GB)는 저장소에 미포함.** 논문 위치·다운로드 경로·샘플수 근거는
> [`pipeline_port/docs/데이터_근거표.md`](pipeline_port/docs/데이터_근거표.md) 참조.
> 저장소에는 **축약 데이터**(STEP4 ML 입력 매트릭스 + Hallmark·KEGG gmt)만 포함되며,
> 이 세트만으로 STEP4(ROC)·5(GSEA)·6(MR)의 핵심을 재현할 수 있습니다.

## 환경 · 설치

패키지 설치는 **R = `install_packages.R`, Python = `requirements.txt`** 두 가지로 일괄 처리합니다.

### R (버전 4.6.1)
```powershell
cd code\R-scripts-working
Rscript install_packages.R     # CRAN + Bioconductor + TwoSampleMR(GitHub) 한 번에 설치
```
> TwoSampleMR(STEP6)는 GitHub 설치라 rate limit(403)에 걸리면 약 1시간 뒤 그 부분만 다시 실행하세요.

### Python (3.13, venv 이름: DKD_AI)
```powershell
cd code\Python-scripts-working
python -m venv DKD_AI
.\DKD_AI\Scripts\Activate.ps1        # 실행정책 막히면: Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
python -m pip install --upgrade pip
pip install -r requirements.txt      # pandas numpy scipy scikit-learn matplotlib xgboost gseapy scanpy harmonypy leidenalg igraph anndata
```

## 실행 방법

### R
```r
setwd("code/R-scripts-working")
source("run_all.R")        # config.R 자동 source, STEP1~7
```

### Python
```powershell
cd code\Python-scripts-working
python run_all.py          # config.py 자동 import, STEP1~7 (STEP2는 R 매트릭스 재사용)
```

## 결과 상세

- **3자 대조표(논문 vs R vs Python)**: [`pipeline_port/docs/3자_대조표_논문_R_Python.md`](pipeline_port/docs/3자_대조표_논문_R_Python.md)
- **재현 정확도(우리 재현 vs 논문)**: [`pipeline_port/docs/재현정확도_대조표.md`](pipeline_port/docs/재현정확도_대조표.md)
- 원자료 CSV: `code/Python-scripts-working/results/python_vs_R_vs_paper.csv`,
  각 STEP 산출물은 `results/step4_paper · step5_gsea · step6_mr · step7_scrna/`.

## 한계 및 근사 구현 사항

- **STEP2 RMA (R 전용)**: `affy/oligo` RMA는 Python 미이식 → R 산출 발현매트릭스를 공유 입력으로 사용.
- **limma → Welch t 근사**: Python STEP1은 eBayes(moderated t)를 Welch t + BH로 근사 → DEG 개수 소폭 차이(방향·비율 일치).
- **ComBat 재사용**: `inmoose` 미설치 시 STEP3 검증셋은 R의 `data.valid.paper.txt`를 직접 사용.
- **결합모델 SVM 포팅 갭**: Python `SVC(rbf, gamma='scale')` ≠ R `e1071` → SVM 검증 AUC 0.405로 이탈(단일 ROC·타 모델은 일치).
- **결합 GLM 조직 이질성 갭**: 훈련(사구체) → 검증(사구체+세뇨관) 전이로 GLM 검증이 논문 0.942 대비 낮음(R 0.826 / Python 0.926).
- **MR LD clumping**: 오프라인 제약으로 독립 clumping 대신 저자 도구변수(Supp Table 8)를 사용(결과가 Supp9와 일치하여 타당).
- **FinnGen(2.1GB) outcome**: 무거워 스킵, GCST90435706으로 핵심 재현.

## 문서 인덱스 (`pipeline_port/docs/`)

| 문서 | 내용 |
|---|---|
| [데이터_근거표.md](pipeline_port/docs/데이터_근거표.md) | 각 데이터의 논문 위치·출처·다운로드·샘플수 근거 |
| [3자_대조표_논문_R_Python.md](pipeline_port/docs/3자_대조표_논문_R_Python.md) | 논문 vs R vs Python 전체 대조 |
| [재현정확도_대조표.md](pipeline_port/docs/재현정확도_대조표.md) | 우리 재현 vs 논문 정확도·한계 |
| [03_R코드_흐름분석.md](pipeline_port/docs/03_R코드_흐름분석.md) | 원본 R 코드 단계별 흐름·이유 |
| [01_필요데이터_목록.md](pipeline_port/docs/01_필요데이터_목록.md) · [02_Python_재구성_가이드.md](pipeline_port/docs/02_Python_재구성_가이드.md) | 데이터 목록 · Python 이식 가이드 |

## 인용 · 라이선스

- **논문**: Lin et al., *Multi-omics and machine learning identify FN1 and ALDH2 as diagnostic biomarkers
  and therapeutic targets in early and late diabetic kidney disease*, **Renal Failure**, 2025 (P25).
- **저자 원본 코드**: GitHub `ljw71865/R-scripts-for-pipeline-reproducibility`.
- **저자 결과·보충자료**: figshare **DOI 10.6084/m9.figshare.30190087** (CC BY 4.0).
- **유전자세트**: MSigDB (Hallmark / KEGG) — Broad Institute, MSigDB 라이선스 조건 준수.
- **공개 데이터**: GEO, OpenGWAS, FinnGen, GWAS Catalog 각 데이터의 이용 약관을 따릅니다.

본 저장소의 재현/이식 **코드**는 연구·교육 목적의 재현연구 산출물이며, 저자 원본 자료의 재사용은
위 CC BY 4.0 및 각 데이터 출처의 라이선스를 따릅니다. 원본 RAW 데이터와 개인/기관 문서는 포함되지 않습니다.
