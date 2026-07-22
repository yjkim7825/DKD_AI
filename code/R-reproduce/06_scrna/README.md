# 06_scrna — 단일세포 분석 (GSE209781: Control 3 + Late DKD 3)

저자 `single-cell RNA-seq analysis 1~4.R` 재현. **어느 세포가 FN1·ALDH2를 발현/주도하나**를 세포 단위로 규명.

## 파이프라인 (저자 4스크립트 = 우리 4파일)

| 파일 | 저자 대응 | 역할 |
| --- | --- | --- |
| `01_scrna_process.R` | scRNA 1.R | QC → 정규화 → PCA → **Harmony** → 클러스터링(res 0.2) → 세포주석(12종) |
| `02_scrna_proportion.R` | scRNA 2.R | 세포비율 + **FN1/ALDH2 발현정량**(wilcox) + EMT 등 ModuleScore |
| `03_scrna_pseudotime.R` | scRNA 3.R | **Monocle 의사시간** (PCT→dPCT) + FN1·ALDH2 추세 |
| `04_scrna_cellchat.R` | scRNA 4.R | **CellChat** 세포통신 (Control vs Late DKD) |
| `reference/` | — | 저자 원본 4개 보존 |

**데이터 흐름**: 10x raw → [01] → `Late_annotated.RDS` → [02] → `Late_with_geneset_scores.RDS` → [03 monocle] / [04 cellchat]

## 핵심 파라미터 (저자 그대로)
- QC: nFeature 300~5000, percent.mt<10
- Normalize LogNormalize(1e4) · HVG vst 2000 · PCA npcs=30
- **Harmony** group.by=patient · Cluster dims=1:30, **resolution=0.2** (12 클러스터)
- Markers min.pct=0.25, logfc=0.5, only.pos
- Monocle: negbinomial.size, num_cells≥10, ~celltype qval<0.01, DDRTree
- CellChat: Secreted Signaling, computeCommunProb(raw.use=T, trim=0.1), min.cells=10

## 세포타입 12종 (저자 주석)
T · **PCT**(근위세관) · EC(내피) · Mono-Mac · **dPCT**(손상 근위세관) · MES(계막) · LOH-DCT · B · Mast · Neut · Plasma1 · Plasma2

## 논문 라이브러리 버전 (코드에 명시)
Seurat **v5.3.0** · harmony **v1.2.3** · monocle **v2.36.0** · CellChat **v2.2.0** · GSVA **v2.2.0**(ssGSEA 186 KEGG) · Cell Ranger 4

## 그리는 figure 패널 (논문 Fig 7 = 말기)
| 패널 | 내용 | 스크립트 · 출력 |
| --- | --- | --- |
| 7A | t-SNE 세포 아틀라스 | 01 → `01.Fig7A_tSNE_celltype.pdf` |
| 7B | 세포비율(Control vs DKD) | 02 → `02.Fig7B_cell_proportion.pdf` |
| 7C | 근위세관 아형 비율 | 02 → `02.Fig7C_PCT_subtype_proportion.pdf` |
| 7D | 마커 버블 DotPlot | 01 → `01.Fig7D_marker_dotplot.pdf` |
| 7E | ALDOB·CUBN·SLC34A1 바이올린 | 02 → `02.Fig7E_PCT_function_genes.pdf` |
| 7F | monocle 유사시간 궤적 | 03 → `03.trajectory_*.pdf` |
| 7G | FN1·ALDH2 바이올린/동태 | 02·03 → `02.Fig7G_*`, `03.FN1_ALDH2_*` |
| (S6/7) | 경로 활성(모듈스코어·ssGSEA) | 02 → `02.EMT_score_*`, `02.ssGSEA_KEGG_*` |
| 세포소통 | CellChat 비교 | 04 → `04.*` |
- ⚠ Fig 7H(Control·Early·Late 단계별) 는 초기(GSE131882)+말기 둘 다 필요 → 말기 단독 불가.

## 실행 (사용자 R, Seurat 필요)
```r
setwd(".../R-reproduce/06_scrna")
source("run_scrna.R")            # 01→02→03→04 (무거움; 단계별 권장)
```
필요 패키지: Seurat, harmony, monocle(v2), CellChat, GSVA, dplyr, ggplot2, ggpubr, patchwork
- 데이터: `data/1-7. GSE209781_RAW/` 각 GSM tar.gz **압축해제** 필요 (10x 폴더). 이미 처리된 `data/processed/GSE209781_annotated.rds` 있으면 01 대신 그걸 `Late_annotated.RDS`로 써도 됨.

## 주의
- ⚠ 저자 scRNA 1.R 후반엔 다른 프로젝트(RCC/전립선) 잔여 코드가 섞여 있음 — DKD와 무관해 재현에서 **제외**함.
- ⚠ Seurat/monocle/CellChat은 R 전용이라 이 프로젝트 샌드박스에선 실행 불가 → **사용자 PC에서 실행**. (코드는 저자+논문 기준으로 작성, 미실행)
- 상세: `../docs/06_scRNA_논문vs재현.md`, `../docs/06_scRNA_코드차이.md`
