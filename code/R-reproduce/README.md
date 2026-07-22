# R-reproduce — 논문 재현 (단계별 번호 폴더)

파이프라인 순서대로 번호를 붙였어요. 각 폴더는 그 단계만 담당해요.

| 폴더 | 단계 | 하는 일 | 입력 → 출력 |
| --- | --- | --- | --- |
| `01_preprocessing` | 전처리 | probe변환·정규화·라벨·ComBat | 원본(CEL/txt) → 발현행렬 |
| `02_deg` | 차등발현 | DEG (limma) | 발현행렬 → DEG 목록 |
| `03_pathway` | 경로분석 | ORA/GSEA/ssGSEA | DEG → 활성 경로 |
| `04_ml` | 머신러닝 | LASSO·ROC | 발현행렬 → FN1·ALDH2 |
| `05_mr` | 멘델리안 | eQTL·MR·IVW | eQTL+GWAS → 인과 유전자 |
| `06_scrna` | 단일세포 | QC·Harmony·pseudotime | scRNA → 세포별 발현 |

## 공통 규칙
- 각 폴더 안: `R/`(기능 함수) + `run_*.R`(실행) + `output/`(결과)
- **원본 코드**(`../R-scripts-for-pipeline-reproducibility`)와 **데이터**(`../data`)는 안 건드림
- 결과는 각 폴더 `output/` 에만 저장

## 실행 순서
```
01_preprocessing → 02_deg → 03_pathway → 04_ml → 05_mr → 06_scrna
```
(앞 단계 출력이 다음 단계 입력)
