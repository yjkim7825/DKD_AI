# 논문 ↔ 우리 재현 그림 대조표

- **논문 본문 그림(Figure 1~)**: 유정님이 가진 **논문 PDF** 안에 있음 (우리가 파일로 가진 건 없음)
- **논문 보충 그림(S1~S10)**: `Figures S1_S10.docx` 안. 일부는 아래처럼 PNG로 추출해둠
- **우리 재현/추출 그림**: 아래 경로 (전부 `code/` 하위)

---

## Figure 1 — 전체 흐름도

| 논문 | 우리 |
| --- | --- |
| 본문 Figure 1 (컬러 3-Step 흐름도) | (흐름도라 재현 안 함) |
| 보충 **Figure S1** (상세 파이프라인) | `data_output/paper_figures/S1_pipeline.png` ✅추출 |

## Figure 2 — DEG + 경로분석

| 패널 | 논문 | 우리 파일 |
| --- | --- | --- |
| 2A 볼케이노(말기vs초기) | 본문 Fig 2A | `R-scripts-working/results/paper_repro/volcano_Late_vs_Early.png` |
| 2B DEG 히트맵 | 본문 Fig 2B | `.../paper_repro/B_DEG_heatmap.png` |
| 2C Hallmark(초기vs정상) | 본문 Fig 2C | `.../paper_repro/C_Early_vs_Control_NES.png` |
| 2D Hallmark(말기vs정상) | 본문 Fig 2D | `.../paper_repro/D_Late_vs_Control_NES.png` |
| 2E Hallmark(말기vs초기) | 본문 Fig 2E | `.../paper_repro/E_Late_vs_Early_NES.png` |
| 2F GSEA 곡선 | 본문 Fig 2F | `.../paper_repro/F_GSEA_curve.png` |

*추가 볼케이노: `volcano_Early_vs_Control.png`, `volcano_Late_vs_Control.png`*

## MR (Mendelian randomization)

| 항목 | 논문 | 우리 파일 (GCST 기준) |
| --- | --- | --- |
| MR 산점도(scatter) | 본문 MR Fig | `data_output/eqtlgen_mr_gcst/scatter_1.pdf`(FN1), `scatter_2.pdf`(ALDH2) |
| forest plot | 본문/보충 | `.../eqtlgen_mr_gcst/forest_1.pdf`, `forest_2.pdf` |
| funnel plot | 보충 **S2** | `.../eqtlgen_mr_gcst/funnel_1.pdf`, `funnel_2.pdf` |
| leave-one-out | 보충 **S3** | `.../eqtlgen_mr_gcst/leaveoneout_1.pdf`, `leaveoneout_2.pdf` |

*FinnGen 버전은 `data_output/eqtlgen_mr/` 폴더에 동일 구성*

## 임상 상관 (Nephroseq eGFR)

| 논문 | 우리 파일 |
| --- | --- |
| 본문 Fig A/B (eGFR 상관 산점도) | `data_output/paper_figures/Nephroseq_eGFR_corr.png` (내가 만든 표+막대) |

## 분자 도킹

| 논문 | 우리 |
| --- | --- |
| 본문 Fig C/D (Vina −7.0 / −9.8 결합) | ❌ 없음 (CB-Dock2 웹툴) |
| 약물 후보 | Table S13 (자료에 있음) |

## Single-cell RNA-seq

| 패널 | 논문 | 우리 파일 |
| --- | --- | --- |
| t-SNE 군집/주석 | 본문 scRNA-A | `R-scripts-working/results/step7_scrna/UMAP_celltype.pdf` |
| 마커 dotplot | 본문 scRNA-D | `.../step7_scrna/DotPlot_FN1_ALDH2.pdf` (FN1·ALDH2만) |
| FN1·ALDH2 violin | 본문 scRNA-G | `.../step7_scrna/Violin_FN1_ALDH2.pdf` |
| 단계별 발현 박스플롯 | 본문 scRNA-H | `.../paper_repro/H_FN1_ALDH2_boxplot.png` |
| 경로 활성(세포별) | 보충 **S6/S7** | 논문: `data_output/paper_figures/S6_pathway_earlyDKD.png`, `S7_pathway_lateDKD.png` / 우리: `.../step7_scrna/KEGG.ssGSEA.scRNA_heatmap.pdf` |
| CellChat 통신망 | 보충 **S9/S10** | `.../step7_scrna/paper_figures/FigureS9_cellchat_earlyDKD.png`, `FigureS10_cellchat_lateDKD.png` |
| 소변(GSE266146) violin | 보충 **S8** | `.../step7_scrna/Violin_FN1_ALDH2_urinary.pdf` |
| pseudotime 궤적 | 본문 scRNA-F | ❌ 없음 (monocle 미실행) |
| FeaturePlot | — | `.../step7_scrna/FeaturePlot_FN1_ALDH2.pdf` |

## 개념 설명용 (참고)

| 그림 | 파일 |
| --- | --- |
| LASSO λ 곡선(개념) | `data_output/paper_figures/LASSO_lambda_curve_concept.png` |

---

## 수치 비교 (논문 vs 우리)

### DEG 개수

| 비교 | 논문 | 우리 |
| --- | --- | --- |
| 말기 vs 초기 | 2,833 (↑1,557/↓1,276) | 3,314 (↑1,755/↓1,559) |
| 말기 vs 정상 | 3,525 | 4,022 |
| 초기 vs 정상 | 390 | 671 |
| 패턴 | 진행할수록 급증 | 동일 ✓ |

### MR (GCST outcome)

| 유전자 | 논문 Supp9 | 우리(eQTLGen) |
| --- | --- | --- |
| FN1 | OR 2.78, p=0.012 | OR 4.53, p=0.012 (방향·유의 일치, 크기는 스케일차) |
| ALDH2 | OR 0.673, p=0.032 | OR 0.57, p=0.097 (방향 일치, 경계) |

### 경로

| | 논문 | 우리 |
| --- | --- | --- |
| 말기 상위 경로 | EMT·염증·자멸사 | EMT·염증 ✓ |
