# DKD Multi-omics 재현 (P25: FN1 & ALDH2)

## 1. 프로젝트 소개

당뇨병성 신장질환(DKD, Diabetic Kidney Disease) 다중오믹스 바이오마커 논문
*Lin et al., "Multi-omics and machine learning identify FN1 and ALDH2 as diagnostic biomarkers
and therapeutic targets in early and late diabetic kidney disease", Renal Failure 2025* (이하 **P25**)
의 분석 파이프라인을 **로컬에서 재현(R)** 하고 **Python으로 이식**한 재현연구 저장소입니다.
bulk 전사체 → DEG → Mendelian Randomization(인과) → LASSO/SVM-RFE 바이오마커 → ROC → GSEA →
단일세포(scRNA)로 이어지는 서사를 STEP 1~7로 구현하여, 논문·R·Python 세 값을 나란히 대조합니다.
핵심 결과(FN1·ALDH2 진단 성능·인과·경로·세포 국소화)는 R과 Python 양쪽에서 재현됩니다.

## 2. 폴더 구조

```
DKD_AI/
├── README.md
├── code/
│   ├── R-scripts-working/            # 재현 R 파이프라인 (STEP1~7, config 중앙화)
│   │   ├── config.R  run_all.R  .here
│   │   ├── 01_deg_gse142025/ … 07_scrna/
│   │   └── results/                  # R 산출물
│   ├── Python-scripts-working/       # Python 이식본 (R과 대칭 구조)
│   │   ├── config.py  run_all.py  requirements.txt  PORT_REPORT.md
│   │   ├── 01_deg_gse142025/ … 07_scrna/
│   │   └── results/                  # Python 산출물 (+ 3자 대조표 CSV)
│   ├── R-scripts-for-pipeline-reproducibility/   # 저자 원본 코드(참조용)
│   └── data/
│       ├── processed/                # 축약 매트릭스(재현용, train/test/valid)
│       └── 4-1..gmt / 4-2..gmt       # Hallmark / KEGG 유전자세트
└── pipeline_port/docs/               # 데이터 근거표 · 재현정확도/3자 대조표
```

원본 RAW 데이터(약 12GB)와 개인/기관 문서는 저장소에 포함되지 않습니다(§3, §7 참조).

## 3. 사용 데이터셋

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

> **원본 12GB는 저장소에 미포함.** 각 데이터의 논문 위치·다운로드 경로·샘플수 근거는
> [`pipeline_port/docs/데이터_근거표.md`](pipeline_port/docs/데이터_근거표.md) 참조.
> 저장소에는 **축약 데이터**(STEP4 ML 입력 매트릭스 train/test/valid + Hallmark·KEGG gmt)만 포함되며,
> 이 세트만으로 STEP4(ROC)·5(GSEA)·6(MR)의 핵심을 재현할 수 있습니다.

## 4. 설치

패키지 설치는 **R = `install_packages.R`, Python = `requirements.txt`** 두 가지로 일괄 처리합니다.

### R (버전 4.6.1) — `Rscript install_packages.R` 로 일괄 설치
```powershell
cd code\R-scripts-working
Rscript install_packages.R     # CRAN + Bioconductor + TwoSampleMR(GitHub) 한 번에 설치
```
> `install_packages.R`가 `BiocManager`/`remotes`부터 준비해 STEP1~7에 필요한 패키지를 모두 설치합니다.
> TwoSampleMR(STEP6)는 GitHub 설치라 rate limit(403)에 걸리면 약 1시간 뒤 그 부분만 다시 실행하세요.

### Python (venv 이름: DKD_AI) — `requirements.txt` 로 일괄 설치
```powershell
cd code\Python-scripts-working
python -m venv DKD_AI
.\DKD_AI\Scripts\Activate.ps1        # 실행정책 막히면: Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
python -m pip install --upgrade pip
pip install -r requirements.txt      # pandas numpy scipy scikit-learn matplotlib xgboost gseapy scanpy harmonypy leidenalg igraph anndata
```
> 검증 환경: Python 3.13. `inmoose`(ComBat)는 미설치 시 STEP3에서 R 산출 매트릭스를 재사용합니다(§7).

## 5. 실행 순서

### R
```r
setwd("code/R-scripts-working")
source("run_all.R")        # config.R 자동 source, STEP1~7
```

### Python
```powershell
cd code\Python-scripts-working
python run_all.py          # config.py 자동 import, STEP1~7
```

| STEP | 내용 |
|---|---|
| 1 | GSE142025 DEG (limma / Welch t) |
| 2 | CEL → RMA 발현매트릭스 (**R 전용**) |
| 3 | 병합 + ComBat 배치보정 |
| 4 | LASSO + SVM-RFE → 단일유전자·결합모델 ROC (핵심) |
| 5 | Hallmark GSEA + KEGG ssGSEA |
| 6 | 2-표본 MR (IVW / MR-Egger / weighted median) |
| 7 | scRNA (Seurat / scanpy) 세포유형별 발현 |

## 6. 핵심 결과 요약

전체 대조는 [`pipeline_port/docs/3자_대조표_논문_R_Python.md`](pipeline_port/docs/3자_대조표_논문_R_Python.md)
및 `code/Python-scripts-working/results/python_vs_R_vs_paper.csv` 참조.

**단일유전자 ROC-AUC (STEP4)**

| 유전자 | 세트 | 논문 | R | Python |
|---|---|---|---|---|
| FN1 | train / valid | 0.911 / 0.911 | 0.909 / 0.915 | 0.909 / 0.915 |
| ALDH2 | train / valid | 0.912 / 0.815 | 0.940 / 0.784 | 0.940 / 0.784 |

**인과·경로·세포 (STEP6·5·7)**

| 항목 | 결과 (R = Python, 논문 대조) |
|---|---|
| MR (GCST90435706) | FN1 OR **2.78**(위험), ALDH2 OR **0.67**(보호) — 저자 Supp9 일치 |
| GSEA 축 | FN1 = EMT/ECM(rho 0.919), ALDH2 = OXPHOS/대사(TCA 0.921) |
| scRNA 국소화 | FN1 = 내피/메산지움, ALDH2 = 근위세뇨관(PCT) |

→ **Python은 동일 입력에서 R과 소수점까지 일치**(단일 ROC·MR·GSEA·scRNA). FN1은 논문과도 거의 정확 재현.

## 7. 한계 및 근사 구현 사항

- **STEP2 RMA (R 전용)**: `affy/oligo` RMA는 Python 미이식 → R 산출 발현매트릭스를 공유 입력으로 사용.
- **limma → Welch t 근사**: Python STEP1은 eBayes(moderated t)를 Welch t + BH로 근사 → DEG 개수 소폭 차이(방향·비율 일치).
- **ComBat 재사용**: `inmoose` 미설치 시 STEP3 검증셋은 R의 `data.valid.paper.txt`를 직접 사용.
- **결합모델 SVM 포팅 갭**: Python `SVC(rbf, gamma='scale')` ≠ R `e1071` → SVM 검증 AUC 0.405로 이탈(단일 ROC·타 모델은 일치).
- **결합 GLM 조직 이질성 갭**: 훈련(사구체) → 검증(사구체+세뇨관) 전이로 GLM 검증이 논문 0.942 대비 낮음(R 0.826 / Python 0.926).
- **MR LD clumping**: 오프라인 제약으로 독립 clumping 대신 저자 도구변수(Supp Table 8)를 사용(결과가 Supp9와 일치하여 타당).
- **FinnGen(2.1GB) outcome**: 무거워 스킵, GCST90435706으로 핵심 재현.

## 8. 인용 및 라이선스

- **논문**: Lin et al., *Multi-omics and machine learning identify FN1 and ALDH2 as diagnostic biomarkers
  and therapeutic targets in early and late diabetic kidney disease*, **Renal Failure**, 2025 (P25).
- **저자 원본 코드**: GitHub `ljw71865/R-scripts-for-pipeline-reproducibility`.
- **저자 결과·보충자료**: figshare **DOI 10.6084/m9.figshare.30190087** (CC BY 4.0).
- **유전자세트**: MSigDB (Hallmark / KEGG) — Broad Institute, MSigDB 라이선스 조건 준수.
- **공개 데이터**: GEO, OpenGWAS, FinnGen, GWAS Catalog 각 데이터의 이용 약관을 따릅니다.

### 라이선스
본 저장소의 재현/이식 **코드**는 연구·교육 목적의 재현연구 산출물입니다. 저자 원본 자료·보충자료의
재사용은 위 CC BY 4.0 및 각 데이터 출처의 라이선스를 따릅니다. 원본 RAW 데이터와 개인/기관 문서는
포함되지 않습니다.
