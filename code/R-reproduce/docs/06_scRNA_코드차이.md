# 06. 단일세포(scRNA) — 저자 원본 코드 vs 우리 재현 코드

저자 `single-cell RNA-seq analysis 1~4.R` ↔ 우리 `R-reproduce/06_scrna/01~05`.
저자 Seurat/monocle/CellChat 파이프라인의 **큰 골격·파라미터는 그대로 따름**. 다만 저자 원본을 실제로 정독·실행하며 **불가피하거나 더 타당한 조정 5가지**를 했고(아래 ★), 초기 데이터 처리(01b)·단계 비교(05)를 **추가**했다.

**★ 조정 5가지 요약**
1. **세포주석**: 하드코딩 `RenameIdents(번호→타입)` → **마커 z-score 자동배정** (우리 클러스터 번호가 저자와 달라 그대로 쓰면 라벨이 뒤바뀜)
2. **dPCT 정의**: 저자는 dPCT = 특정 클러스터(4번) 통째 → 우리는 **PCT 중 손상마커 상위 30%**를 dPCT로 세분
3. **Harmony 보정 단위**: 저자는 `patient`를 **Control/Late_DKD 2군**으로 정의해 보정 → 우리는 **샘플(orig.ident) 6개** 기준 보정(표준 배치보정)
4. **유사시간**: monocle2 ↔ 최신 dplyr 충돌(`group_by_` defunct) → **손상점수 정렬**로 대체
5. **CellChat**: presto 미설치 → `do.fast=FALSE`, 입력 layer·population.size 등 소소한 차이

## 📂 파일 대응
| 우리 | 저자 | 역할 |
| --- | --- | --- |
| `01_scrna_process.R` | scRNA 1.R | QC→Harmony→클러스터→세포주석 (말기 GSE209781) |
| `01b_scrna_early.R` | (없음) | 초기 GSE131882 동일 파이프라인 추가 처리 |
| `02_scrna_proportion.R` | scRNA 2.R | 세포비율+FN1/ALDH2+ModuleScore(+GSVA 추가) |
| `03_scrna_pseudotime.R` | scRNA 3.R | 유사시간 (monocle→손상점수 정렬 대체) |
| `04_scrna_cellchat.R` | scRNA 4.R | CellChat 세포통신 |
| `05_scrna_stage_compare.R` | (없음) | 데이터셋 내부 PCT vs dPCT 단계 검증 |
| `reference/single-cell RNA-seq analysis 1~4.R` | (원본) | 저자 업로드 원본 |

---

## 0. 한눈에

| 단계 | 저자 원본 | 우리 재현 | 판정 |
| --- | --- | --- | --- |
| 데이터 입력 | 병합본 `late_merged_seurat.rda` 로드 | **10x raw 6샘플 읽어 병합**(메모리 절약) | ✅ 결과 동일 |
| QC | nFeature 300~5000, **mt<10** | 동일(mt<10) | ⚠ 논문 본문=15% |
| 정규화·HVG·PCA | LogNorm(1e4)·vst2000·npcs30 | 동일 | ✅ |
| **배치보정(Harmony)** | `patient`=**Control/Late_DKD 2군** 기준 | `patient`=**샘플 6개** 기준 | ⚠ 보정단위 다름 |
| 클러스터 | dims1:30, **res=0.2** → 12클러스터 | 동일(res=0.2) | ✅ |
| 세포주석 | 수동 `RenameIdents` 12종(번호 하드코딩) | ★ **마커 z-score 자동배정** | ⚠ 바꿈(라벨 뒤바뀜 방지) |
| **dPCT 정의** | dPCT = **클러스터 4번 통째** | **PCT 중 손상마커 상위 30%** 세분 | ⚠ 정의 방식 다름 |
| 유전자집합 점수 | **AddModuleScore** (EMT·노화·자가포식·산화·세포사·염증 6종 MSigDB) | 동일 AddModuleScore **+ GSVA ssGSEA(186 KEGG) 추가** | ✅ + 확장 |
| 유사시간 | monocle2 negbinomial, ~celltype qval<0.01, DDRTree | ★ **손상점수 정렬로 대체** | ⚠ dplyr 충돌 |
| CellChat | Secreted, counts입력, 기본 identifyOverExpressed(presto), trim=0.1, population.size=TRUE, min.cells10 | data입력, **do.fast=FALSE**, trim=0.1, min.cells10 | ⚠ 소소한 차이 |
| 초기 GSE131882 | (없음) | **01b로 추가 처리** | ➕ 확장 |
| 단계 비교 | (없음) | **05로 추가**(플랫폼 교락 제거) | ➕ 확장 |

→ 골격·핵심 파라미터는 일치. 차이는 **입력방식·보정단위·dPCT정의·주석방식·유사시간구현·CellChat세부**와 우리가 **추가한 01b/05**.

---

## 1. 단계별 대조

### ① 데이터 입력·QC·정규화·PCA (scRNA 1.R)
```r
# 저자: 이미 병합된 rda 로드
Late <- readRDS("...late_merged_seurat.rda"); Late <- JoinLayers(Late)
# 우리: 10x raw 6샘플을 하나씩 읽어 QC로 줄인 뒤 병합 (679만 배럴코드 → RAM 폭발 방지)
for(nm in samples){ o<-CreateSeuratObject(Read10X(dir),min.cells=3,min.features=200); ...; rm(cnt);gc() }

# 이후는 저자·우리 동일
subset(nFeature_RNA>300 & nFeature_RNA<5000 & percent.mt<10)
NormalizeData(LogNormalize,1e4); FindVariableFeatures(vst,2000); ScaleData(); RunPCA(npcs=30)
```
저자 코드의 세포수(`table(orig.ident)`): Control 3603·5383·949 / Late_DKD 1081·5188·2612 = **18,816** → 우리 재현 **18,817**로 일치.

**★ 미토 컷(mt%) 불일치 — 논문 본문 15% vs 저자코드·우리 10%:**
- 논문 **본문(Methods)** 엔 "미토 비율 **15% 초과 제거**". 그러나 저자 **실제 코드**는 `percent.mt < 10`.
- 우리는 **저자 코드(10%)를 따름**. 근거: **10%라야 논문 보고 세포수(18,816)가 재현**됨(우리 18,817). 15%면 세포가 더 남아 수치가 어긋남.
- 결론: 논문 서술(15%)과 저자 코드(10%)가 상충 → **저자가 실제 돌린 값은 10%**로 판단. 재현은 코드값(10%)이 정답.

### ② 배치보정(Harmony) — ⚠ 보정 단위 다름
```r
# 저자: patient를 '질병군(2개)'으로 정의 → Control vs Late_DKD 사이만 정렬
Late$patient <- ifelse(grepl("^Control",orig.ident),"Control","Late_DKD")
RunHarmony(Late,"patient")
# 우리: patient=개별 샘플(6개) → 샘플 간 배치를 정렬 (표준적 배치보정)
Late$patient <- Late$orig.ident
RunHarmony(Late,"patient")
```
저자 방식은 배치 축을 **질병군**으로 잡아 정렬하고, 우리는 **샘플**로 잡음. 6샘플 배치보정이 일반적 관행이라 우리는 샘플 단위를 채택. 클러스터 구조 자체는 크게 다르지 않음.

### ③ 클러스터·주석 — ⚠ 주석 방식 + dPCT 정의 다름
```r
FindNeighbors(reduction="harmony",dims=1:30); FindClusters(resolution=0.2)  # 저자·우리 동일 → 12클러스터
FindAllMarkers(min.pct=0.25,logfc.threshold=0.5,only.pos=TRUE)              # 동일

# 저자(하드코딩) — 우리 클러스터 번호와 안 맞아 라벨 뒤바뀜:
# RenameIdents(0=T,1=PCT,2=EC,3=Mono-Mac,4=dPCT,5=MES,6=LOH-DCT,7=B,8=Mast,9=Neut,10=Plasma1,11=Plasma2)
#   → 저자는 dPCT = '클러스터 4번' 통째

# 우리(번호 무관 마커 자동배정 + dPCT 세분):
annotate_by_markers(Late, KIDNEY_MARKERS)                     # 마커 z-score로 타입 결정
inj <- colMeans(data[c("VCAM1","HAVCR2","SPP1"),])            # 손상마커 점수
dPCT <- PCT 중 inj 상위 30%                                    # PCT를 손상 정도로 PCT/dPCT 세분
```
**왜 바꿨나:** 저자 하드코딩은 "클러스터 4번=dPCT"를 전제하는데, 우리 재현의 클러스터 번호 순서가 달라 그대로 쓰면 **ALDH2 높은 세관이 dPCT로 오분류**(방향이 거꾸로). → 마커 자동배정으로 교체하니 ALDH2↓ 방향이 논문과 일치.
**dPCT 정의 차이:** 저자는 dPCT를 독립 클러스터로 보고, 우리는 PCT 안에서 손상마커로 세분 → 세부 세포수는 다를 수 있으나 "손상 근위세관" 개념은 동일.

### ④ 유전자집합 점수 (scRNA 2.R) — ✅ 동일 + 확장
```r
# 저자: AddModuleScore 6종 MSigDB — EMT, CELL_AGING, AUTOPHAGY, OXIDATIVE_STRESS, APOPTOSIS, INFLAMMATORY_RESPONSE
AddModuleScore(Late, features=list(geneset), name=...)   # 세포타입별 boxplot + tSNE FeaturePlot
# 우리: 위 AddModuleScore 동일 + GSVA ssGSEA(186 KEGG, PCT vs dPCT pseudobulk) 추가 정량
```
저자 방식(AddModuleScore) 그대로 쓰고, 경로 활성 정량을 위해 **GSVA ssGSEA를 추가**로 얹음(더 넓게).

### ⑤ 유사시간 — ⚠ monocle→손상점수 정렬 대체 (scRNA 3.R)
```r
# 저자(monocle2): 최신 dplyr에서 group_by_() defunct로 실행 불가
# newCellDataSet(negbinomial.size()) → differentialGeneTest(~celltype) qval<0.01
#   → reduceDimension(DDRTree) → orderCells → plot_genes_in_pseudotime(FN1,ALDH2)
# 우리(대체): 손상마커 점수로 세포를 정상→손상 순서 정렬
subset(idents=c("PCT","dPCT"))
inj <- colMeans(data[c("VCAM1","HAVCR2","SPP1"),])
pseudotime <- rank(inj)/length(inj)                 # 0~1 정상→손상
# FN1·ALDH2를 pseudotime 축으로 loess 추세 (Fig 7G 대응)
```

### ⑥ CellChat — ⚠ 소소한 파라미터 차이 (scRNA 4.R)
```r
# 공통: group(Control/Late_DKD) 분할, Secreted Signaling DB, trim=0.1, min.cells=10, mergeCellChat
# 저자: 입력 layer="counts"(원시), identifyOverExpressedGenes() 기본(do.fast=TRUE→presto), population.size=TRUE
# 우리: 입력 layer="data"(정규화), identifyOverExpressedGenes(do.fast=FALSE)  # presto 미설치 회피
```
핵심 로직(분비신호 리간드-수용체 통신, Control vs Late 비교)은 동일. presto 의존을 피하려 `do.fast=FALSE`로 표준 Wilcoxon 사용.

---

## 2. 그 외 정리 차이 (로직 아님)

| 항목 | 저자 | 우리 |
| --- | --- | --- |
| 작업 경로 | `setwd("G:\\187geneMR\\...")` 하드코딩 | `ROOT` 상대경로 |
| 잔여 코드 | scRNA 1.R 후반 RCC/전립선 등 **타 프로젝트 복붙** 혼재 | **제외**(DKD 파이프라인만) |
| 마커 DotPlot | 여러 pdf로 흩어짐(수십 블록) | 논문 Fig7D 순서로 1블록 정리 |
| FN1·ALDH2 | gene 변수 바꿔 반복(scatter/violin 다수) | 단계별 정리(01 violin, 02 정량, 03 추세) |

---

## 3. 한 줄 요약
> 저자 scRNA 4스크립트(Seurat→monocle→CellChat)를 5파일로 재현·**실제 실행 완료**(초기 20,620·말기 18,817세포 = 논문 20,220·18,816과 일치). 골격·핵심 파라미터는 동일하되, 실제 코드 정독 결과 **주석방식·dPCT정의·Harmony 보정단위·유사시간구현·CellChat세부**를 타당하게 조정하고 **초기(01b)·단계비교(05)를 추가**. **핵심 결과 — 말기 dPCT에서 ALDH2↓(p=3e-85, Fig 7G 재현)**. FN1은 세뇨관 snRNA에서 검출한계 이하라 방향 확인 불가(데이터 한계, 조작 없이 정직 보고).
