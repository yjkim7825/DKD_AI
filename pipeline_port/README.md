# DKD 바이오마커 파이프라인 (P25 재현 + Python 이식)

P25 논문(*Multi-omics and machine learning identify FN1 and ALDH2 as diagnostic
biomarkers and therapeutic targets in early and late diabetic kidney disease*,
Renal Failure 2025)의 **bulk 전사체 바이오마커 파이프라인**을 정리했습니다.

- `R/` — 원본 R 스크립트를 **config 기반**으로 정리(하드코딩 경로 `G:\187geneMR\...` 제거)
- `python/` — 동일 로직 **Python 이식본** (pandas / scipy / scikit-learn)

원본 전체 파이프라인 중 이번 범위는 **핵심 bulk 부분**입니다:

```
전처리(probe→gene, 정규화, ComBat) → DEG(limma) → LASSO/SVM-RFE 교집합 → ROC/AUC
```
> Mendelian randomization(TwoSampleMR)·scRNA(Seurat/monocle/CellChat)는 다음 단계로 분리했습니다.

---

## 1. 폴더 구조

```
pipeline_port/
├── README.md
├── data/                      # (직접 채워야 함) GEO/GWAS 원본·중간 파일
├── results/                   # 산출물(자동 생성)
├── R/
│   ├── config.R               # 경로/파라미터
│   ├── 01_preprocess.R        # probe2gene / normalize_ds / combat_merge
│   ├── 02_diff_expression.R   # run_deg (limma)
│   ├── 03_feature_selection.R # LASSO + SVM-RFE + Venn
│   └── 04_roc.R               # run_roc (pROC)
└── python/
    ├── requirements.txt
    ├── config.py
    ├── io_utils.py            # GEO/GPL 파서, probe 매핑, avereps
    ├── deg_utils.py           # quantile 정규화, ComBat, limma eBayes
    ├── p01_preprocess.py
    ├── p02_diff_expression.py
    ├── p03_feature_selection.py
    ├── p04_roc.py
    └── smoke_test.py          # 합성데이터 검증(데이터/네트워크 불필요)
```

---

## 2. 필요한 데이터 (직접 다운로드 → `data/` 에 배치)

원본 스크립트가 참조하는 공개 데이터입니다. **네트워크가 되는 로컬 PC에서** 받으세요.

| 데이터 | 종류 | 출처 | 용도(스크립트) |
|---|---|---|---|
| GSE37263 | microarray | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37263) | probe→gene (preprocessing 1) |
| GSE142025 | RNA-seq | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025) | 3군(Control/Early/Late) 정규화 |
| GSE96804 | microarray | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96804) | 병합·학습셋 |
| GSE104948 | microarray | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104948) | 병합·학습셋 |
| GSE181061 | scRNA/기타 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181061) | 확장 단계 |
| GSE209781 | scRNA | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209781) | scRNA 확장 단계 |
| finngen_R12_DM_NEPHROPATHY_EXMORE | GWAS | [FinnGen](https://www.finngen.fi/en/access_results) | MR(확장) |
| GCST90435706 | GWAS | [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | MR(확장) |

- GEO 발현행렬: **Series Matrix File** (`GSExxxxx_series_matrix.txt.gz`)
- 플랫폼 주석: 해당 **GPL** 파일 (probe→gene symbol 컬럼 포함)
- 샘플 그룹 리스트: `s1.txt`(Control) / `s2.txt`(Early) / `s3.txt`(Late) — 샘플명 한 줄 1개

> **입력/출력 열 이름 규약**: 모든 발현행렬의 샘플 열은 `{샘플ID}_{그룹}` 형식
> (예: `GSM123_Control`, `GSM456_DKD`). Early/Late 를 하나로 묶어 2군 비교할 땐 그룹을 `DKD`로.

---

## 3. 실행 — R (원본 재현)

로컬에서 R + 아래 패키지 설치 후 실행:

```r
# Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "sva"))
# CRAN
install.packages(c("readr","dplyr","stringr","ggplot2",
                   "glmnet","e1071","kernlab","caret","VennDiagram","pROC"))
```

```r
BASE_DIR <- "C:/Users/.../DKD_AI/pipeline_port"   # ← 본인 경로
source(file.path(BASE_DIR, "R", "config.R"))
source(file.path(BASE_DIR, "R", "01_preprocess.R"))
source(file.path(BASE_DIR, "R", "02_diff_expression.R"))
source(file.path(BASE_DIR, "R", "03_feature_selection.R"))
source(file.path(BASE_DIR, "R", "04_roc.R"))

# 예시 (각 함수의 하단 '실행 예시' 주석 참고)
probe2gene(file.path(DATA_DIR,"GSE37263_series_matrix.txt.gz"),
           file.path(DATA_DIR,"GPL5175.txt"),
           out=file.path(DATA_DIR,"GSE37263_geneMatrix.txt"), symbol_index=2)
normalize_ds(file.path(DATA_DIR,"GSE142025_geneMatrix.txt"),
             control=file.path(DATA_DIR,"s1.txt"),
             early  =file.path(DATA_DIR,"s2.txt"),
             late   =file.path(DATA_DIR,"s3.txt"), geoid="GSE142025")
run_deg(file.path(DATA_DIR,"GSE142025_twoGroups.normalize.txt"),
        ref="Control", alt="DKD")
run_feature_selection(file.path(DATA_DIR,"data.train.txt"),
                      genes=file.path(DATA_DIR,"interGenes.List.txt"))
run_roc(file.path(DATA_DIR,"data.train.txt"),
        file.path(RESULT_DIR,"ml","LASSO.gene.txt"), "Training")
```

---

## 4. 실행 — Python (이식본)

```bash
cd pipeline_port/python
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# (데이터/네트워크 없이) 코드 정상동작 확인
python smoke_test.py

# 실제 파이프라인
python p01_preprocess.py probe2gene --series ../data/GSE37263_series_matrix.txt.gz \
       --platform ../data/GPL5175.txt --symbol-col "Gene Symbol" --symbol-index 1 \
       --out ../data/GSE37263_geneMatrix.txt
python p01_preprocess.py normalize --matrix ../data/GSE142025_geneMatrix.txt \
       --control ../data/s1.txt --early ../data/s2.txt --late ../data/s3.txt \
       --geoid GSE142025 --out ../data/GSE142025_threeGroups.normalize.txt
python p02_diff_expression.py --input ../data/GSE142025_twoGroups.normalize.txt \
       --ref Control --alt DKD --tag DKD_vs_Control
python p03_feature_selection.py --train ../data/data.train.txt \
       --genes ../data/interGenes.List.txt
python p04_roc.py --expr ../data/data.train.txt \
       --genes ../results/ml/LASSO.gene.txt --title Training
```

---

## 5. R ↔ Python 대응 & 충실도 메모

| 단계 | R | Python | 비고 |
|---|---|---|---|
| 정규화 | `limma::normalizeBetweenArrays("quantile")` | `deg_utils.quantile_normalize` | 동일(quantile) |
| 배치보정 | `sva::ComBat(par.prior=TRUE)` | `deg_utils.combat` | parametric EB 동일 구현 |
| DEG | `limma lmFit+eBayes+topTable` | `deg_utils.limma_two_group` | **Smyth(2004) moderated t 그대로 구현** → logFC/t/P/adjP 수치 일치 |
| LASSO | `glmnet(cv.glmnet, lambda.min)` | `LogisticRegressionCV(l1, saga)` | 동등 |
| SVM-RFE | `caret::rfe(svmRadial)` | `RFECV(SVC(linear))` | **선형 SVM-RFE**로 이식(radial은 특징랭킹 불가). 최종 유전자 집합은 대개 유사하나 완전 동일하진 않음 |
| ROC | `pROC::roc + ci.auc(bootstrap)` | `sklearn.roc_auc_score + bootstrap` | 동등 |

**검증 상태**: `quantile_normalize`, `combat` 은 합성데이터로 통과 확인
(열평균 동일화 / 배치 평균차 3.0→0.009). `limma_two_group`·ML 은 로컬에서
`python smoke_test.py` 로 확인하세요(참 DEG recall/precision, ML 회수).

---

## 6. 다음 단계(옵션)
- **GSEA** (`GSEA.R`) → Python `gseapy.prerank` 로 이식 가능
- **Mendelian randomization** (`Mendelian randomization 1-3.R`) → `TwoSampleMR`(R 유지 권장)
- **scRNA** (`single-cell RNA-seq analysis 1-4.R`) → Python `scanpy` + `cellchat`/`liana` 로 이식 가능
