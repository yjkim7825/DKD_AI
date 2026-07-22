# 03_pathway — 경로분석 (GSEA)

DEG 결과를 경로(기능) 단위로 해석. 저자 `GSEA.R` 로직 재현.

## 실행
```r
source("run_gsea.R")   # R에서
```

## 입력 (복사 없이 경로 참조)
- DEG: `../02_deg/output/GSE142025_3group/diff_{Early_vs_Control,Late_vs_Control,Late_vs_Early}.txt`
- 유전자셋: `../../data/4-1. Hallmark(50)`, `4-2. KEGG_legacy(186)`

## 과정 (gsea_func.R)
1. `logFC` 로 유전자 순위 (내림차순)
2. `GSEA(logFC, TERM2GENE=gmt, pvalueCutoff=1)`
3. `p.adjust < 0.05` 경로만 필터 → NES>0(켜짐)/NES<0(꺼짐)
4. `gseaplot2` 시각화

## 출력 (output/)
- `GSEA.<DB>.<비교>.txt` — 전체 경로 (NES, p.adjust)
- `GSEA.<DB>.<비교>.{Up,Down}.pdf` — 상위 경로 그림
- `GSEA_summary.txt` — 비교×DB 유의경로 개수

## 구조
```
03_pathway/
├── R/gsea_func.R    # GSEA 함수 (저자 로직)
├── run_gsea.R       # 실행 (3비교 × Hallmark/KEGG)
└── output/          # 결과
```
※ ssGSEA(세포별 경로점수)는 single-cell 단계(06_scrna)에서 `AddModuleScore`로 수행.
