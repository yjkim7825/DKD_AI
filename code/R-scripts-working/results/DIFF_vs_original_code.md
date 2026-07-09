# DIFF_vs_original_code.md — 원본 저자 코드 ↔ 우리 코드 1:1 대조

> 대조 기준 (둘 다 읽기 전용):
> - **원본 저자 코드**: `../R-scripts-for-pipeline-reproducibility/` (구현 authoritative)
> - **논문 원문**: `docs/P25_Multi-omics ... FN1 and ALDH2.pdf`
>
> 분류: **완전일치** / **정렬함**(우리 코드를 원본에 맞춤) / **불가피한 차이**(원본 하드코딩·오프라인 LD·RMA CDF·입력매트릭스 등).
> 경로/데이터 소스는 원본의 `G:\...` 하드코딩을 따르지 않고 우리 `config.R` 방식 유지(규칙).

---

## 0. 스크립트 매핑

| STEP | 우리 스크립트 | 원본 저자 스크립트 | figshare |
|---|---|---|---|
| 1 DEG (GSE142025) | `01_deg_gse142025/02_deg.R` | `data preprocessing 2.R`(정규화) + `differential expression analysis.R`(DEG) | Supp Table 3/4/5 |
| 1 DEG (GSE96804, complementary) | `01_deg_gse142025/03_gse96804_deg.R` | `differential expression analysis.R`(2군 템플릿 재사용) | 논문 p.4 |
| 2 CEL→RMA | `02_rma_microarray/01~03_*.R` | **없음**(원본은 series_matrix부터: `data preprocessing 1.R`) | — |
| 3 병합+ComBat | `03_merge_combat/01_merge_combat.R` | `data preprocessing 3.R` | — |
| 4 LASSO/SVM-RFE/ROC | `04_ml_roc/01_paper_design_lasso_svm_roc.R` | `machine learning modeling 1.R`(LASSO) + `2.R`(SVM-RFE) | Supp 8/9/10/12 |
| 5 GSEA/ssGSEA | `05_gsea/01_gsea.R` | `GSEA.R` + (figshare) `hallmark.gsea.R` | GSEA.result.* |
| 6 MR | `06_mr/01_mr.R` | `Mendelian randomization 1~3.R` | Supp 6/7/8/9/11 |
| 7 scRNA | `07_scrna/01_scrna.R` | `single-cell RNA-seq analysis 1~4.R` + (figshare) `scRNA.Seurat0.7.R` | — |

---

## 1. 세부 대조 (STEP별)

| STEP | 원본 저자 코드 방식 | 우리 코드 방식 | 논문 서술 | 차이 여부 | 정렬조치 / 이유 |
|---|---|---|---|---|---|
| **1 정규화** | `data preprocessing 2.R`: 단일 `GSE142025_matrix.txt` → `avereps` → auto-log2(qx 판정) → `normalizeBetweenArrays` → s1/s2/s3.txt 로 Control/Early/Late 라벨 | 36개 per-sample txt 병합(Symbol 교집합) → `avereps` → 동일 auto-log2 → `normalizeBetweenArrays` → 파일명 접두사 N/B/A 로 라벨 | GSE142025 early/late 비교의 DEG 원천 | **불가피한 차이**(입력) | 정규화·라벨 로직은 원본과 동일. **입력만** 단일 GEO matrix→per-sample 병합(단일 matrix 미보유). DEG 3,314 vs 논문 2,833(70% 겹침) |
| **1 DEG** | `differential expression analysis.R`: limma, `logFCfilter=0.585`, `adj.P.Val=0.05`, `Type=strsplit(colname,"_")[2]`, `model.matrix(~0+factor(Type))`, `makeContrasts(alt-ref)`, `topTable(adjust='fdr', number=2e5)`, sig=`abs(logFC)>0.585 & adj.P<0.05` | 동일: limma, `LOGFC_FILTER=0.585`, `ADJP_FILTER=0.05`, `grp=sub(".*_","",colname)`, `~0+g`, `makeContrasts(alt-ref)`, `topTable(number=Inf, adjust='fdr')`, 동일 sig 필터 | Fig 2A, Supp 3/4/5 | **완전일치** | 원본 2군 템플릿(DKD-Control)을 3대비(LvE/LvC/EvC)로 실행. 임계·방향·함수인자 동일 |
| **1 GSE96804 DEG** | `differential expression analysis.R` 를 GSE96804 에 적용(DKD-Control) | `03_gse96804_deg.R`: 동일 로직으로 GSE96804 DKD vs Control | 논문 p.4 "complementary limma v3.64.1, adjP<0.05,\|logFC\|>0.585" | **정렬함**(추가) | 논문 명시 DEG 였으나 파이프라인에 없어 원본 DEG 로직 그대로 신규 구현. FN1 Up(+2.10)/ALDH2 저발현방향 확인 |
| **2 CEL→RMA** | **원본 없음** — `data preprocessing 1.R` 는 series_matrix + GPL 프로브주석(`str_split(" // ")[2]`)→`avereps` | GSE96804=oligo RMA(HTA-2_0), GSE30529=affy RMA(U133A_2), 104948/54=affy RMA+brainarray ENTREZG CDF | Methods: 8 GEO 데이터 전처리 | **불가피한 차이** | 원본은 series_matrix 시작(우리는 미보유, 제공 데이터=CEL+CDF). 마이크로어레이 표준 RMA 로 구현. 표본수·그룹 논문 Supp1 과 대조(아래) |
| **2 104948/54 Control 정의** | (원본 없음) | `keepCodes=c("DN","LD")` — TN 제외, DN=DKD/LD=Control | Supp1 Control=21 | **정렬함** | LD(21)=논문 21 정확 일치. TN(5) 제외. 4개 데이터 전부 Supp1 표본수 일치(all TRUE) |
| **3 병합+ComBat 로직** | `data preprocessing 3.R`: 폴더 내 모든 txt → `Reduce(intersect)` 교집합 → `header[1]_` 접두사 cbind → `batchType=i` → `ComBat(allTab, batch, par.prior=TRUE)` → preNorm/normalize 저장 | 동일: 교집합→`{GSE}_` 접두사 cbind→batch=데이터셋→`ComBat(par.prior=TRUE)`→preNorm/txt. + 배치 내 무분산 유전자 제거 가드 | Methods: ComBat 배치보정 | **완전일치** | 병합·ComBat 인자 동일. 차이는 폴더 자동수집→세트 명시(경로규칙) + 무분산가드(안정화, 이번 제거 0) |
| **3 훈련/검증 구성** | 원본 ML 경로 `GSE96804_104948`, 입력 `data.train/test.txt`(사구체 vs 세뇨관 암시) | **논문 Figure4 따름**: 훈련=GSE96804 단독 / 검증=ComBat(104948+104954)=`data.valid.paper.txt`. (원본 네이밍 호환 data.train/test.txt 도 참고 생성) | Fig4: train=GSE96804, valid=독립셋(104948+104954) | **정렬함(논문 채택)** | 논문 vs 원본 상충 → **논문 Figure4 채택**(사용자 지시 통일). valid C42/D29/71 = 논문 일치(all TRUE). 원본 네이밍도 병행 유지 |

---

## 2. 데이터셋 표본수 대조 (STEP2 산출 vs 논문 Supp Table 1)

| 데이터 | 논문 Supp1 (Control/DKD) | 우리 RMA (Control/DKD) | 판정 |
|---|---|---|---|
| GSE96804 | 20 / 41 | 20 / 41 | ✅ 일치 |
| GSE30529 | 12 / 10 | 12 / 10 | ✅ 일치 |
| GSE104948 | 21 / 12 | **21 / 12** | ✅ 일치 (2026-07-09 정렬: TN 제외, Control=LD만) |
| GSE104954 | 21 / 17 | **21 / 17** | ✅ 일치 (동일 정렬) |

---

## 3. STEP4 세부 대조 (machine learning modeling 1·2.R)

| 항목 | 원본 저자 코드 | 우리 코드 | 논문 서술 | 차이 | 정렬조치/이유 |
|---|---|---|---|---|---|
| LASSO | `glmnet(family="binomial",alpha=1)` + `cv.glmnet(type.measure="deviance",nfolds=10)`, `coef(s=lambda.min)`, 절편 제거, `set.seed(123)` | 동일(인자·시드 동일) | glmnet v4.1.9, 10-fold CV, λmin | **완전일치** | 결과 6개 = {ALDH2,CDKN1B,FN1,IFI44L,VNN2,XAF1}, 논문 5/6(IFI44L↔TSPYL5) |
| SVM-RFE | `rfe(caretFuncs, method="cv", sizes=c(2:8), preProcess=center/scale)`, `y=as.numeric(as.factor(group))`, `methods="svmRadial"`(오타) | 동일 + `methods`→`method="svmRadial"` 정렬 | (논문 본문 언급 없음) | **정렬함**(오타 수정) | 결과 {CA2,CDKN1B,VNN2}. 원본 회귀형 RFE 그대로 |
| Venn 교집합 | `intersect(LASSO, SVM-RFE)` → interGenes.txt | 동일 | — | **완전일치** | 교집합={CDKN1B,VNN2} → **FN1/ALDH2 미포함**(원본 코드 특성) |
| **최종 유전자 선택** | 코드=SVM-RFE∩LASSO | 논문 방법 병행: **LASSO→ROC(AUC>0.8 both)** → {ALDH2,FN1,CA2} | "AUC>0.8 in both training & test → multivariate" | **논문≠원본코드** | **논문 방법이 FN1·ALDH2 선택**(코드 SVM-RFE 교집합은 아님). 둘 다 산출, 논문을 정본 |
| 단일 ROC | `roc()`+`ci.auc(method="bootstrap")`, train=data.train, valid=data.test, LASSO 유전자 루프 | 동일 함수, train=GSE96804/valid=data.valid.paper, 후보10 전체 | Fig4 D–G, AUC+95%CI | **정렬함(셋 논문)** | 함수 동일. 훈련/검증 셋만 논문 Figure4 채택 |
| 결합모델 | (원본 스크립트 범위 밖 — 논문 Fig4H,I) | GLM/RF/SVM/XGBoost FN1+ALDH2 | Fig4H,I multivariate | 논문 전용 | 논문 Fig4H,I 재현 |

### STEP4 논문값 대조 (LD-only 71 검증 반영)
| 지표 | 논문 | 우리 | 판정 |
|---|---|---|---|
| FN1 AUC train/valid | 0.911 / 0.911 | 0.909 / 0.871 | train 정확, valid 근사(−0.04) |
| ALDH2 AUC train/valid | 0.912 / 0.815 | 0.940 / **0.807** | valid 0.807≈0.815 (LD-only 정렬로 0.784→0.807 개선) |
| LASSO 6개 | {CDKN1B,ALDH2,FN1,XAF1,TSPYL5,VNN2} | {…,IFI44L 대신 TSPYL5 제외} | 5/6 일치 |
| 결합 GLM valid | 0.942 | 0.820 | 갭 존재(RF 0.927/XGB 0.898 근접) |

---

## 4. STEP5 세부 대조 (GSEA.R + hallmark.gsea.R)

| 항목 | 원본 저자 코드 | 우리 코드 | 논문 서술 | 차이 | 정렬조치/이유 |
|---|---|---|---|---|---|
| Hallmark GSEA | `hallmark.gsea.R`: **diff_(유의 DEG)** logFC → `GSEA(TERM2GENE, pvalueCutoff=1)`, msigdbr H | 동일 로직, **입력 diff_ 로 정렬**(이전 all_), gmt 파일(4-1)로 msigdbr 대체 | Fig2C-E | **정렬함** | 입력을 원본과 동일 diff_ 로 변경. EMT +2.78/+2.53(논문 +2.52/+2.59), OXPHOS -1.71/-1.82(-1.62/-2.02) |
| KEGG GSEA | `GSEA.R`: 동일 diff_ logFC → KEGG gmt GSEA | **신규 추가**(이전 파이프라인 누락), gmt(4-2) | Fig2 | **정렬함(추가)** | 원본 KEGG GSEA 를 그대로 구현. ECM_RECEPTOR +1.95(p.adj 9e-4), TCA/OXPHOS 하향 |
| gene set 출처 | msigdbr(H) / `c2.cp.kegg.Hs.symbols.gmt` | 제공 gmt 4-1/4-2 | MSigDB | **불가피한 차이** | 오프라인 → 제공 gmt(동등). GSEA p<1e-10 은 eps 클리핑(정밀도 경고, 결론 무관) |
| ssGSEA | (원본 GSEA.R 은 KEGG GSEA; ssGSEA 스크립트 미제공) | KEGG ssGSEA(GSVA) on GSE96804 + FN1/ALDH2 Spearman | Fig2H,I 단계별 경로점수 | 논문 전용(원본코드 없음) | 논문 Fig2H,I 방법 재현. FN1↔ECM 0.919, ALDH2↔TCA 0.921 |

---

## 5. STEP6 세부 대조 (Mendelian randomization 1~3.R)

| 항목 | 원본 저자 코드 | 우리 코드 | 논문 서술 | 차이 | 정렬조치/이유 |
|---|---|---|---|---|---|
| 노출 | `read_exposure_data(clump=FALSE)` from exposure.F.csv(저자 준비) | 저자 도구변수(Supp8) ∩ eQTL vcf 에서 beta/se/allele/eaf 추출 | eQTL, LD clump r²<0.001 | **불가피한 차이** | 오프라인 LD API 불가 → 저자 instrument(이미 clump) 사용. 결과 Supp9 정확 일치 |
| outcome 읽기 | FinnGen `rsids/beta/sebeta/alt/ref/pval/af_alt`; GCST | 동일 열 매핑 + 필요 SNP만 fread | FinnGen R12 + GCST90435706 | **완전일치** | 열 매핑 동일 |
| **outcome 필터** | `dat[dat$pval.outcome>5e-06,]` | **동일 추가**(NA 가드) | 역인과 방지 | **정렬함** | 누락됐던 필터 추가. nsnp 3/14 = Supp9 유지 |
| MR 방법 | `mr(c("mr_ivw","mr_egger_regression","mr_weighted_median"))` + `generate_odds_ratios` | 동일 | IVW 주분석 + Egger/median | **완전일치** | GCST IVW FN1 b=1.021/nsnp3, ALDH2 b=−0.395/nsnp14 = **Supp9 소수점 일치** |
| 이질성/다면발현 | `mr_heterogeneity`, `mr_pleiotropy_test` | 동일 | Supp11 | **완전일치** | FN1 절편 p=0.99, ALDH2 p=0.70(다면발현 없음) |
| **IVW 필터(MR2.R)** | IVW p<0.05 & 3방법 OR방향 일치 & 다면발현 p>0.05 | **동일 추가** | — | **정렬함(추가)** | FinnGen 통과 0(FN1/ALDH2 비유의=Supp6 정합), GCST 통과 {ALDH2,FN1} |
| 민감도 | `mr_singlesnp`, `mr_leaveoneout`(+plots) | single-SNP/leave-one-out CSV 추가(plot 생략) | leave-one-out | **정렬함(추가)** | 표만 저장(무거운 PDF 생략). 원본 forest(MR3.R)는 cosmetic |

### STEP6 논문 대조 (Supp9 IVW, GWAS Catalog)
| 유전자 | 논문 b / nsnp | 우리 b / nsnp | OR | 판정 |
|---|---|---|---|---|
| FN1 | 1.0212 / 3 | **1.0212 / 3** | 2.777 | ✅ 완전 일치 (위험 인과 OR>1) |
| ALDH2 | −0.3954 / 14 | **−0.3954 / 14** | 0.673 | ✅ 완전 일치 (보호 인과 OR<1) |

---

## 6. STEP7 세부 대조 (scRNA.Seurat0.7.R + single-cell RNA-seq analysis 1~4.R)

| 항목 | 원본 저자 코드 | 우리 코드 | 논문 서술 | 차이 | 정렬조치/이유 |
|---|---|---|---|---|---|
| 데이터 | GSE209781 6샘플(pre-merged rda) | GSE209781 6샘플(Read10X) | scRNA late-stage | **완전일치(데이터)** | 동일 데이터셋. tar→Read10X(원본은 저자 rda) |
| QC | `nFeature>300 & <5000 & percent.mt<10` | 동일 | — | **완전일치** | 18,817 세포(원본 로그 유사) |
| 정규화/HVG/PCA | LogNormalize(1e4)/vst 2000/PCA 30 | 동일 | — | **완전일치** | 인자 동일 |
| Harmony 변수 | `RunHarmony("patient")` (Control/DKD 2군) | `orig.ident`(6샘플) | batch 통합 | **의도적 차이** | 샘플단위 배치보정(표준, 원본 group단위보다 세분). FN1/ALDH2 국소화 결론 무관 |
| 해상도 | `FindClusters(resolution=0.2)` | **0.2 로 정렬**(이전 0.5) | — | **정렬함** | 원본과 동일 0.2 → 15 클러스터 |
| 주석 | 수동 12클러스터→세포타입 | 마커세트 argmax 자동주석 | 9 신장+7 면역 세포 | **방법 차이** | 자동(재현성). 마커세트는 원본 주석 마커와 동일 계열 |
| FN1/ALDH2 | DotPlot/Violin by celltype | 동일 + 세포타입별 mean/%expr 표 | 내피_系膜_FN1 / PCT_ALDH2 | **완전일치(결론)** | FN1=Endothelial, ALDH2=PCT 재현 |

### STEP7 논문 대조 (FN1/ALDH2 세포 국소화)
| 유전자 | 논문 | 우리(res 0.2) | 판정 |
|---|---|---|---|
| FN1 | 내피/系膜/손상PCT 고발현 | **Endothelial** (mean 1.713, 73.8%), Mesangial 2위 | ✅ 정합 |
| ALDH2 | PCT(근위세뇨관) 고발현 | **PCT** (mean 1.533, 75.9%) | ✅ 정합 |

---

## 7. 종합 분류

- **완전일치**: STEP1 DEG, STEP3 병합·ComBat 로직, STEP4 LASSO·Venn·이질성/다면발현, STEP5 gene set 랭킹, STEP6 MR 방법(Supp9 소수점 일치), STEP7 QC·정규화·FN1/ALDH2 결론.
- **정렬함(원본에 맞춤)**: STEP2 Control=LD(21), STEP4 SVM-RFE 오타·단일 ROC 셋, STEP5 diff_입력·KEGG GSEA 추가, STEP6 outcome필터·IVW필터·민감도 추가, STEP7 해상도 0.2.
- **논문 채택(논문≠원본코드)**: STEP3/4 훈련=GSE96804·검증=104948+104954(Figure4), STEP4 최종선택=ROC AUC>0.8(코드 SVM-RFE 교집합 아님).
- **불가피한 차이**: STEP1 GSE142025 입력(per-sample 병합 vs 단일 matrix), STEP2 CEL→RMA(원본은 series_matrix), STEP5 gmt(msigdbr 대체), STEP6 LD clumping(저자 instrument), STEP7 Harmony 변수(sample vs group).
