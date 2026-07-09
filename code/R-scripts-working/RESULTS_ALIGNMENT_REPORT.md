# RESULTS_ALIGNMENT_REPORT.md — 논문/원본코드 정렬 + 중간산출물 저장 (2026-07-09)

P25(*FN1·ALDH2 DKD multi-omics*, Renal Failure 2025) R 파이프라인을 **논문 원문 + 원본 저자 코드**
(`../R-scripts-for-pipeline-reproducibility/`)와 1:1 대조해 정렬하고, 각 STEP 중간산출물을
`results/` STEP별 폴더에 전부 저장한 작업 요약. 상세 대조표는 `results/DIFF_vs_original_code.md`,
결정 근거는 `DECISIONS.md`.

> 절대규칙 준수: 원본 데이터(`../data` 번호폴더)·원본 스크립트 **무수정**. 경로는 `config.R` 참조(하드코딩 0).
> 산출물은 `results/`·`../data/processed/` 로만.

---

## STEP별 변경 요약

### STEP1 — DEG
- 산출물을 `results/step1_deg/` 로 정리(STEP 구분). STEP5 읽기 경로 동기화.
- **논문 명시 GSE96804 complementary limma DEG 추가**(`03_gse96804_deg.R`) — MR 통합 근거.
  FN1 = DKD 고발현(logFC +2.10, Up), ALDH2 = 저발현 방향(컷 경계).
- 대조: GSE142025 Late_vs_Early 3,314(논문 2,833, 70% 겹침 — 입력매트릭스 차이, 불가피).
- 원본 `differential expression analysis.R` 와 임계·방향·함수 **완전일치**.

### STEP2 — CEL→RMA
- 원본엔 RMA 없음(series_matrix 시작) → 우리 RMA 는 추가분(제공 데이터 CEL+CDF).
- **미해결 #1 해소**: GSE104948/54 Control = **LD(21)만**(TN 제외) = 논문 Supp1 정확 일치.
  → `keepCodes=c("DN","LD")` 정렬, 4개 데이터 표본수 전부 Supp1 일치(all TRUE).
- QC 중간산출물 → `results/step2_rma/`(요약/샘플분포/표본수대조).

### STEP3 — 병합 + ComBat
- 병합·ComBat 로직 원본 `data preprocessing 3.R` 와 **완전일치**.
- **논문 Figure4 로 설계 통일**: 훈련=GSE96804 / 검증=ComBat(104948+104954)=`data.valid.paper.txt`(71, C42/D29).
  (원본 ML 네이밍 data.train/test.txt 도 병행 생성.) 논문 대조 all TRUE.
- 중간산출물 → `results/step3_merge/`.

### STEP4 — LASSO/SVM-RFE/ROC/결합
- 입력을 STEP3 `data.valid.paper.txt`(71) 로 정리.
- 원본 ML 코드 반영: **SVM-RFE + Venn 교집합 + ci.auc(bootstrap) CI 추가**(원본 오타 `methods=`→`method=`).
- **논문 vs 원본코드 차이 규명**: 최종 FN1/ALDH2 선택은 **논문 ROC AUC>0.8 방법**(={ALDH2,FN1,CA2})에서 나오며
  코드의 SVM-RFE∩LASSO(={CDKN1B,VNN2})에서는 아님 — 둘 다 산출, 논문 채택.
- LD-only 효과: **ALDH2 검증 0.784→0.807**(논문 0.815 근접), FN1 0.909/0.871. LASSO 5/6 일치.
- 중간산출물 → `results/step4_paper/`.

### STEP5 — GSEA/ssGSEA
- **Hallmark GSEA 입력 all_→diff_** 로 원본 정렬. EMT +2.78/+2.53(논문 +2.52/+2.59), OXPHOS -1.71/-1.82.
- **KEGG GSEA 추가**(원본 `GSEA.R`, 이전 누락). ECM_RECEPTOR +1.95(p.adj 9e-4).
- ssGSEA(논문 Fig2H,I) 유지: FN1↔ECM 0.919, ALDH2↔TCA 0.921.
- 중간산출물 → `results/step5_gsea/`.

### STEP6 — MR
- 원본 MR1·2.R 누락분 반영: **outcome 필터 `pval.outcome>5e-06`**, **IVW 필터(MR2.R)**, **single-SNP/leave-one-out**.
- **저자 Supp9 소수점 일치**: FN1 b=1.0212/nsnp3/OR2.777, ALDH2 b=−0.3954/nsnp14/OR0.673.
  IVW 필터 통과: FinnGen 0(비유의=Supp6 정합), GCST {ALDH2,FN1}. 다면발현 없음(FN1 p=0.99/ALDH2 p=0.70).
- 중간산출물 → `results/step6_mr/`.

### STEP7 — scRNA
- **해상도 0.5→0.2** 원본 정렬(→18,817 세포, 15 클러스터). QC/정규화/PCA 원본 일치.
- Harmony 변수는 orig.ident(샘플단위, 원본 group단위보다 표준) 유지 — 결론 무관.
- FindAllMarkers(원본 인자)·마커점수·대조 CSV 추가.
- **FN1 = Endothelial(1.713/73.8%), ALDH2 = PCT(1.533/75.9%)** — 논문 정합.
- 중간산출물 → `results/step7_scrna/`.

---

## 종합 판정

| 분류 | STEP |
|---|---|
| **완전일치** | STEP1 DEG, STEP3 ComBat, STEP4 LASSO/Venn/이질성, STEP5 랭킹, STEP6 MR(Supp9 소수점), STEP7 결론 |
| **정렬함(원본에 맞춤)** | STEP2 Control=LD21, STEP4 SVM-RFE/CI, STEP5 diff_·KEGG GSEA, STEP6 outcome/IVW필터·민감도, STEP7 res 0.2 |
| **논문 채택(논문≠원본코드)** | STEP3/4 훈련·검증 구성(Fig4), STEP4 최종선택(ROC AUC>0.8) |
| **불가피한 차이** | STEP1 입력매트릭스, STEP2 RMA(vs series_matrix), STEP5 gmt, STEP6 LD clumping, STEP7 Harmony 변수 |

**핵심 결론 전부 재현**: FN1=위험(EMT/ECM/내피·系膜), ALDH2=보호(OXPHOS/대사/PCT),
진단 ROC·MR 인과(Supp9 정확 일치). 논문과 불가피한 차이는 강제 정렬하지 않고 근거와 함께 기록.

## 각 STEP 산출물 위치
`results/step1_deg/ · step2_rma/ · step3_merge/ · step4_paper/ · step5_gsea/ · step6_mr/ · step7_scrna/`
각 폴더에 `compare_paper_vs_ours.*.csv`(논문/원본값 vs 우리값 대조) 포함.
