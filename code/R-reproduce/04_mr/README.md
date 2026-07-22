# 04_mr — 멘델리안 무작위화 (MR)

DEG 후보 유전자(FN1·ALDH2)가 **진짜 DKD를 일으키는지(인과)** SNP으로 판별. 저자 `Mendelian randomization 1/2/3.R` 재현.

## 실행 (저자처럼 1→2→3, CSV로 데이터 전달)
```r
source("01_mr_run.R")     # MR1: 노출→결과→mr→CSV 저장
source("02_mr_filter.R")  # MR2: CSV 읽어 IVW 필터 → IVW.filter.csv
source("03_mr_forest.R")  # MR3: CSV 읽어 forest plot
# 또는 한 번에:
source("run_mr.R")        # 위 3개 순서대로
```

## 전 단계 연결
- **입력 후보**: `../02_deg/output/GSE142025_3group/diff_Late_vs_Control.txt` — FN1·ALDH2가 DEG인지 확인 후 MR
- 그 두 유전자의 eQTL(노출) + DKD GWAS(결과)로 인과 검증

## 입력 (`../../data/`)
| 종류 | 파일 |
| --- | --- |
| 노출(eQTL) | `2-1. eqtl-a...115414.vcf`(FN1), `2-2...111275.vcf.gz`(ALDH2) |
| 도구변수 | `6.../Supplementary Table 8.csv` (저자 clump된 IV) |
| 결과(GWAS) | `3-1. finngen_R12...gz`, `3-2. GCST90435706/...tsv.gz` |
| 논문 정답 | `6.../Supplementary Table 9.csv` |

## 과정 (저자처럼 파일 3개, CSV로 연결)
1. **`01_mr_run.R`** (MR1): `parse_vcf`(노출)→`harmonise`→`mr`(IVW·Egger·median)→민감도 → **CSV 저장**
2. **`02_mr_filter.R`** (MR2): 01의 `table.MRresult`·`table.pleiotropy` **읽어** IVW p<0.05 & 방향일치 & 다면발현 p>0.05 → `IVW.filter.csv`
3. **`03_mr_forest.R`** (MR3): `table.MRresult` **읽어** IVW forest plot
- `R/mr_func.R` = 01이 쓰는 공용 함수(parse_vcf·read_outcome·run_mr)

## 출력 (`output/`)
- `table.MRresult.{FinnGen,GCST}.csv` — MR 결과(OR)
- `IVW.filter.*.csv` — 인과 통과 유전자
- `table.{heterogeneity,pleiotropy,singleSNP,leaveoneout}.*.csv` — 민감도
- `forest_MR_GCST.png` — forest plot
- `compare_paper_vs_ours.MR_IVW.csv` — 논문 Supp9 대조

## 결과 (GCST outcome, 논문 완전 일치)
- **FN1**: OR 2.78, p=0.012 (위험↑) / **ALDH2**: OR 0.67, p=0.032 (보호↓)

## 구조
```
04_mr/
├── 01_mr_run.R      # MR1: 노출→결과→mr→CSV
├── 02_mr_filter.R   # MR2: IVW 필터
├── 03_mr_forest.R   # MR3: forest plot
├── R/mr_func.R      # 01 공용 함수
├── run_mr.R         # 1·2·3 한 번에
└── output/          # 결과 (CSV로 단계 간 전달)
```
※ eQTLGen 원본에서 IV를 새로 만드는 버전(OR 4.53)은 `04_MR_결과요약.md` 참고 (노출 소스 차이).
