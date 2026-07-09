# DECISIONS.md — 결정·미해결 기록 (R-scripts-working)

원본 데이터/스크립트는 수정하지 않는다. 여기에는 "왜 이렇게 했는가 / 미해결로 남긴 것"만 기록.

---

## STEP 1 — GSE142025 DEG (2026-07-07, 완료)

### 결과
- `01_GSE142025_to_matrix.R` → 17,184 genes × 36 samples (Control 9 / Early 6 / Late 21, 파일명 접두사 N/B/A).
- `02_GSE142025_DEG.R` → 자동 log2 판정에서 **log2 스킵**(입력이 이미 log 스케일) 후 `normalizeBetweenArrays`.
- **Late_vs_Early: 3,314 DEGs (up 1,755 / down 1,559).**
- 부가: Late_vs_Control 4,022 / Early_vs_Control 671.

### 논문(figshare Supplementary Table 3) 대조
- 논문값: **up 1,557 / down 1,276 / 총 2,833** (정답표 직접 카운트로 확정).
- 유전자 목록 겹침: **1,975개 공유(논문 목록의 70% 회수)**, 우리만 1,339 / 논문만 858.
- 방향(up>down)·비율 일관. **정확 재현은 아님(70% 겹침).**

### 차이 원인 가설 (미해결 — 나중 검토)
1. **입력 매트릭스 차이(가장 유력).** 작업본은 36개 per-sample txt(이미 샘플별 정규화된 log 값)를
   Symbol **교집합 병합**. 논문/원본은 GEO 제공 단일 `GSE142025_matrix.txt`에서 시작.
   시작 값·유전자 세트가 달라 임계값 경계에서 DEG 수가 벌어짐.
2. **정규화·모델 차이.** RNA-seq인데 limma를 count가 아닌 log 값에 바로 적용.
   원본과 동일 로직이나, 본래 RNA-seq는 `voom`(가중치) 경로가 더 적합할 수 있음 → GSEA 단계에서 재검토.
3. 개수를 2,833에 억지로 맞추는 튜닝은 하지 않음(정답 overfit 방지).

### 처리 방침
- **미해결로 남김.** FN1/ALDH2 LASSO/SVM/ROC(핵심)는 GSE96804 등 다른 데이터로 진행되므로 영향 없음.
- GSEA(STEP 5)에서 Late_vs_Early 순위 리스트를 쓸 때 이 차이를 다시 확인.
- 필요 시 정공법: GEO 원본 `GSE142025_matrix.txt` 확보 후 재대조, 또는 voom 경로 비교.

### 논문 정렬 패스 (2026-07-09) — 산출물 재정리 + 누락 DEG 추가
- **출력 정리**: STEP1 산출물을 `results/` 루트 → **`results/step1_deg/`** 로 이동(STEP 구분 저장 규칙). 구 루트 파일 삭제.
  STEP5(`05_gsea/01_gsea.R`)가 이 DEG 를 읽으므로 읽기 경로도 `step1_deg/` 로 동기화.
- **논문값 대조 CSV 추가**: `step1_deg/compare_paper_vs_ours.GSE142025.csv`
  (논문 Late_vs_Early 2,833/up1557/down1276 vs 우리 3,314/up1755/down1559).
- **누락 DEG 추가(논문 명시)**: `03_gse96804_deg.R` — 논문 Methods(p.4) *"Complementary differential expression
  analysis of the GSE96804 dataset using limma (v3.64.1), adjusted p<0.05, |logFC|>0.585"* 를 그대로 구현.
  이 DEG 가 MR 통합(DKD 고발현∩risk / 저발현∩protective)의 근거이며 FN1/ALDH2 방향을 정의함.
  - 결과: **GSE96804 DKD vs Control DEG 654개(up 336/down 318)**.
  - **FN1: logFC +2.10, adjP 8.3e-8 → DKD 고발현(Up)** = 논문 "FN1 = DKD 고발현·MR risk" 방향 정확 일치.
  - **ALDH2: logFC −0.42, adjP 6.6e-6 → 저발현 방향이나 |logFC|<0.585 컷 미달**(이 데이터셋 한정).
    방향(DKD 저발현·protective)은 논문과 일치하나 GSE96804 단독에선 유의컷 경계 아래 — 정직 기록.
    (논문 최종 후보는 저자 명시 10개를 사용하므로 STEP4 에 영향 없음.)
  - 산출: `step1_deg/DEG_{all,diff}_GSE96804_DKD_vs_Control.txt`, `GSE96804.DKD_{high,low}_genes.txt`,
    `GSE96804.FN1_ALDH2_direction.csv`, `vol_GSE96804_DKD_vs_Control.pdf`.
- **논문과 불가피한 차이(강제 정렬 안 함)**: GSE142025 입력이 GEO 단일 matrix 가 아니라 36개 per-sample txt 병합
  → Late_vs_Early 3,314 vs 논문 2,833(70% 겹침). 위 "차이 원인 가설" 그대로 유지.

---

## STEP 2 — CEL → RMA 발현매트릭스 (진행 중)

> ⚠️ 원본 저장소에는 **CEL→RMA 단계가 없음**(원본은 이미 만들어진 발현매트릭스부터 시작).
> 따라서 이 STEP은 [Reviewer] 1:1 대조 대상이 아니라 **마이크로어레이 표준 관례**(oligo/affy RMA)에
> 따라 새로 구현. 검증은 산출물 차원·log 스케일·그룹 수·핵심 유전자(FN1/ALDH2) 존재로 수행.

### GSE96804 (2026-07-07, 완료) — `03_GSE96804_to_matrix.R`
- 플랫폼: Affymetrix **HTA-2_0** (whole-transcript, GPL17586) → `oligo::rma(target="core")` 사용
  (affy 패키지가 아니라 oligo. CLAUDE.md의 "Gene ST"는 근사 표기, 실제 HTA-2_0).
- 주석: 추가 주석 패키지 설치 없이 **`getNetAffx(eset,"transcript")`** 의 geneassignment에서 첫 심볼 추출.
- 그룹: 파일명이 전부 `DN` 이라 `GEOquery::getGEO` 메타데이터로 판정.
  - 대조군 라벨이 "unaffected portion of tumor nephrectomies" 라 일반 패턴(normal/control)에 안 걸림
    → **`grepl("diabet")` 로 DKD를 양성 정의**, 나머지를 Control 로. 결과 **Control 20 / DKD 41** (기지 구성과 일치).
- 산출: `../data/processed/GSE96804.labeled.txt` — **33,720 유전자 × 61 샘플**, log2 범위 1.79–13.08, FN1·ALDH2 존재.
- 겪은 버그(수정됨): `grepl`/`ifelse` 가 벡터 이름을 버려 그룹 매칭 0/61 발생 → `names(metaGrp)` 재부여로 해결.
- 부작용 메모: `getGEO(destdir=OUT_DIR)` 가 `GSE96804_series_matrix.txt.gz` 를 processed/ 에 캐시함(원본 아님, 무해).

### GSE30529 (2026-07-07, 완료) — `04_GSE30529_to_matrix.R`
- 플랫폼: **HG-U133A_2** (3'-IVT, GPL571) → `affy::ReadAffy(cdfname="hgu133a2cdf")` + `affy::rma()`.
- 주석: `hgu133a2.db` PROBEID→SYMBOL (probe당 심볼 1개), `avereps`.
- 그룹: `GEOquery`. source_name = "Tubuli of control kidney"(12) / "Tubuli of DKD kidney"(10) → **Control 12 / DKD 10** (기지 구성 일치).
- 산출: `../data/processed/GSE30529.labeled.txt` — **13,041 유전자 × 22 샘플**, log2 범위 3.18–14.87, FN1·ALDH2 존재.

### GSE104948 / GSE104954 (2026-07-07, 완료) — `05_GSE104948_104954_to_matrix.R`
- 한 GSE 안에 **HG-U133A + HG-U133_Plus_2 두 플랫폼 혼재** → 제공된 **brainarray ENTREZG 커스텀 CDF**로 플랫폼별 RMA.
  - `make.cdf.env` 로 CDF 환경 구축(CDF .gz 는 원본 미변경, scratchpad 로 gunzip 복사 후 사용).
  - 커스텀 CDF env 를 `.GlobalEnv` 에 clean 이름으로 assign → `getCdfInfo` 가 찾도록 `raw@cdfName` 설정.
- 병합: 플랫폼별 Entrez(_at) 발현을 **공통 Entrez 교집합(12,135개)** 으로 cbind → Entrez→SYMBOL(org.Hs.eg.db) → avereps.
- 샘플 선택: 파일명 질환코드에서 **DN=DKD, LD/TN=Control 만** 사용(FSGS/HT/IgA/MCD/MGN/RPGN/SLE/TMD 제외).
- 산출:
  - `GSE104948.labeled.txt` — **12,042 유전자 × 38** (Control 26 / DKD 12), log2 2.88–14.74, FN1·ALDH2 존재.
  - `GSE104954.labeled.txt` — **12,042 유전자 × 43** (Control 26 / DKD 17), log2 2.57–14.47, FN1·ALDH2 존재.
- ⚠️ **알려진 한계(미해결):** 한 GSE 안 두 플랫폼을 RMA 후 그냥 cbind → **GSE 내부 플랫폼 배치효과**가 남음.
  paper/GEO 관례상 GSE 단위 1파일로 처리했고, STEP 3 ComBat 은 데이터셋(파일) 단위 배치만 교정.
  ROC(STEP 4)가 이상하면 여기서 플랫폼을 별도 배치로 분리하거나 ComBat 공변량 추가 검토.

### ✅ 미해결 이슈 #1 해소 (2026-07-09) — Control = LD only
1. **GSE104948/54 Control 수: 논문 Supp1 = 21 = LD(21) 정확 일치로 정렬 완료.**
   - 데이터 카운트: 104948 = DN 12 / **LD 21** / TN 5, 104954 = DN 17 / **LD 21** / TN 5.
   - **LD(living donor)만 = 21 = 논문 Supp Table 1 Control 수와 정확 일치.** TN(tumor nephrectomy 정상부, 5개)은
     논문이 제외한 것으로 확정 → `03_gse104948_104954.R` `keepCodes=c("DN","LD")` 로 정렬(TN 제외).
   - 재실행 결과: **GSE104948 C21/D12(33) , GSE104954 C21/D17(38)** — 논문 Supp1 완전 일치(compare CSV all TRUE).
   - 저위험(검증셋). 다운스트림(STEP3 병합·data.valid.paper.txt)은 이 21-Control 매트릭스로 재생성.
2. **GSE 내부 U133A + Plus2 배치효과 미교정.**
   - 한 GSE 안 두 플랫폼을 RMA 후 cbind 만 함. STEP 3 ComBat 은 데이터셋(파일) 단위 배치만 교정 → GSE 내부
     플랫폼 배치는 남음. STEP 4 ROC 이상 시 플랫폼을 별도 배치로 분리하거나 ComBat 공변량 추가 검토.

### STEP 2 종합
- 산출 4종 완료: GSE96804(33,720×61) / GSE30529(13,041×22) / GSE104948(12,042×38) / GSE104954(12,042×43).
- 네 매트릭스 모두 log2 스케일·FN1/ALDH2 존재 확인. STEP 3(병합+ComBat) 준비됨.

### 논문/원본 정렬 패스 (2026-07-09)
- **원본 대조**: 원본엔 CEL→RMA 단계 자체가 없음. 원본 `data preprocessing 1.R` 는 series_matrix + GPL
  프로브주석(`str_split(" // ")[2]`)→avereps 방식 → **우리 RMA 는 원본 대응 없는 추가분**(제공 데이터=CEL+CDF 라 RMA 가 타당).
  `results/DIFF_vs_original_code.md` STEP2 행에 "불가피한 차이"로 기록.
- **QC 중간산출물 저장**(무거운 RMA 재실행 없이 기존 매트릭스에서 요약): `results/step2_rma/`
  — `RMA_matrix_summary.csv`(차원/그룹/log2/FN1·ALDH2), `RMA_sample_distribution.csv`(샘플별 Q1/median/Q3),
  `compare_paper_vs_ours.samples.csv`(논문 Supp1 표본수 대조).
- **표본수 대조**: 96804(C20/D41)·30529(C12/D10) 논문 정확 일치. 104948/54 는 DKD 일치, Control 26 vs 논문 21
  (LD+TN vs LD만 — 미해결 #1 유지, 강제 정렬 안 함).

---

## STEP 3 — 병합 + ComBat (2026-07-07, 완료) — `06_merge_combat.R`

### 세트 구성 결정 (중요)
- **훈련 `data.train.txt` = ComBat(GSE96804 + GSE104948)** — 둘 다 **사구체(glomeruli)**.
- **검증 `data.test.txt`  = ComBat(GSE104954 + GSE30529)** — 둘 다 **세뇨관(tubule)**.
- 근거: 원본 ML 스크립트(`machine learning modeling 1·2.R`) 작업경로가 **`...\GSE96804_104948`** 이고
  입력이 `data.train.txt`/`data.test.txt`. 조직(사구체/세뇨관) 구분으로도 일관.
- ⚠️ **CLAUDE.md STEP 3 메모("훈련 GSE96804 / 검증 104948·104954")와 다름.** 원본 스크립트가 authoritative라
  그쪽을 따름. (CLAUDE.md 는 대략적 핸드오프 메모.) 만약 실제 논문 구성이 다르면 세트만 바꿔 재실행하면 됨.

### 방식 (원본 `data preprocessing 3.R` 로직)
- 세트 내 데이터셋들의 **공통 유전자 교집합** → 데이터셋명 접두사(`{GSE}_`)로 cbind → `ComBat(batch=데이터셋)`.
- 열이름 = `{GSE}_{GSM}_{Control|DKD}` (뒤 `_Control`/`_DKD` 접미사로 STEP4가 라벨 추출 — 원본과 호환).
- ComBat 전 배치 내 무분산 유전자 제거 가드 추가(이번엔 제거 0). 원본과 차이: 폴더 자동수집 대신 파일목록 명시.

### 결과
- **data.train.txt**: 11,302 유전자 × 99 샘플 (Control 46 / DKD 53), log2 2.14–13.80, FN1·ALDH2 존재.
  - 내역: GSE96804(C20/D41) + GSE104948(C26/D12).
- **data.test.txt**: 11,988 유전자 × 65 샘플 (Control 38 / DKD 27), log2 2.80–14.61, FN1·ALDH2 존재.
  - 내역: GSE104954(C26/D17) + GSE30529(C12/D10).
- 참고용 `data.train.preNorm.txt` / `data.test.preNorm.txt`(ComBat 전)도 함께 저장.

### 위험·한계
- STEP 2 의 GSE 내부 플랫폼 배치효과(104948/54)는 여기서 교정되지 않음(위 미해결 이슈 2 참조).
- 훈련은 사구체 2종, 검증은 세뇨관 2종 → 조직 차이가 곧 train/test 도메인 차이. ROC 가 검증에서 크게
  떨어지면 조직 차이/배치 재검토(STEP 4).

### 논문/원본 정렬 패스 (2026-07-09) — 설계 통일 + 중간산출물
- **병합·ComBat 로직**: 원본 `data preprocessing 3.R` 와 **완전일치**(교집합→접두사 cbind→batch=데이터셋→
  `ComBat(par.prior=TRUE)`→preNorm/txt). 우리 추가분은 배치 내 무분산 유전자 제거 가드(이번 제거 0).
- **훈련/검증 구성 = 논문 Figure4 로 통일**(정본): 훈련=GSE96804 단독 / 검증=ComBat(104948+104954).
  원본 ML 은 data.train/test.txt(사구체 vs 세뇨관)를 쓰지만 **논문 Figure4 를 채택**(사용자 지시). 원본 네이밍
  호환용 data.train/test.txt 도 병행 생성(참고).
- **LD-only 반영으로 검증셋 재생성**: `data.valid.paper.txt` = Control 42 / DKD 29 / **71 샘플**
  (이전 26-Control 시절 81 → 21-Control 로 71). 논문 대조 CSV all TRUE.
- **중간산출물** → `results/step3_merge/`: `merge_combat_summary.csv`(4세트 차원/그룹/ComBat여부/FN1·ALDH2),
  `sample_median_pre_post_ComBat.csv`(샘플별 ComBat 전후 중앙값), `compare_paper_vs_ours.design.csv`(설계 대조).
- 산출 매트릭스: `processed/` 에 data.train.paper(=GSE96804,ComBat생략)/data.valid.paper/data.train/data.test (+ *.preNorm).

---

## STEP 4 — LASSO + SVM-RFE 교집합 → ROC (2026-07-07, 완료)
`07_step4_candidate_lists.R`(후보목록) + `08_step4_ml.R`(ML). 원본 `machine learning modeling 1·2.R` 로직.
두 갈래를 **분리 저장**: `results/step4_A_repro/`, `results/step4_B_explore/`, 비교표 `results/step4_AB_compare.csv`.

### 후보 목록 (재현성 위해 results/ 에 interGenes.List.txt 저장)
- **A(재현)** = Supp4 DEG(Late DKD vs Control) ∩ MR-union(Supp6∪7) = **63개**(train 피처 45), FN1·ALDH2 포함.
- **B(탐색)** = 전체 DEG union(Supp3/4/5), MR 필터 없음 = **4,323개**(train 피처 2,050), FN1·ALDH2 포함.

### 결과
| 갈래 | 후보 | train피처 | LASSO | SVM-RFE | 교집합 | FN1 | ALDH2 |
|---|---|---|---|---|---|---|---|
| A(재현) | 63 | 45 | 12 | 45(전체) | **12** | ✅ | ✅ |
| B(탐색) | 4323 | 2050 | 35 | 8 | 1(DUSP1) | ❌ | ❌ |
- **A 교집합(12)**: ALDH2, CA2, CD83, FN1, GMIP, IL12RB1, NET1, SP140L, TPP1, TRAF3IP3, XYLT1, ZFAND5.
- **핵심 재현 성공**: A 에서 **FN1·ALDH2 가 LASSO∩SVM-RFE 교집합에 생존.**
- ⚠️ 단, A 의 SVM-RFE 가 `optVariables=45`(전체)를 최적으로 반환 → SVM 이 추가 축소를 못 함.
  따라서 A 교집합 12 = 사실상 **LASSO 12 그대로**. 논문은 최종 2개({FN1,ALDH2})로 좁혔으나 우리는 12개.
  (seed/버전/rfe sizes 차이. FN1/ALDH2 생존이라는 핵심 결론엔 영향 없음.)
- B 는 사전 MR 필터를 빼자 교집합이 DUSP1 1개로 붕괴, FN1/ALDH2 탈락 → **MR 사전필터가 FN1/ALDH2 선택에 필수적**임을 보여줌(A 접근 타당성 방증).

### ROC-AUC (train=사구체 병합, test=세뇨관 병합; AUC 는 유전자 고유값이라 A/B 동일)
| 유전자 | AUC_train | AUC_test |
|---|---|---|
| **FN1** | 0.842 | **0.950** |
| **ALDH2** | 0.851 | 0.696 |
- FN1: 훈련·검증 모두 강함(검증 0.95). ALDH2: 훈련 강함(0.85), 검증 약화(0.70).
- ⚠️ **논문 보고 AUC 값과의 수치 대조는 미완**: 논문 ROC 는 본문 그림(Figures)에 있고 오프라인 수치 미확보.
  필요 시 `6. Article related data/Figures S1_S10.docx` 또는 논문 그림에서 값 추출해 대조.

### SVM-RFE 보정 (2026-07-07) — `09_step4_svmrfe_fix.R`
- 문제: v1 SVM-RFE 가 `optVariables=45`(전체) 반환 → 원본이 y 를 `as.numeric()` 로 바꿔 **회귀(RMSE)** RFE
  였고 전체 세트가 최적으로 선택됨.
- 수정: y 를 **factor(Control/DKD) 분류** RFE 로 실행 + `pickSizeTolerance(tol=2)` 로 성능 근접 시 소수 피처 선호.
- 결과: **optsize=4** (Accuracy 크기별 2→0.849, 3→0.888, **4→0.942(최적)**, 45→0.931).
  - SVM-RFE.v2 = {CA2, ALDH2, FN1, ZFAND5}.
  - **교집합 재계산 LASSO(12) ∩ SVM-RFE.v2(4) = {ALDH2, CA2, FN1, ZFAND5} (4개)** — **FN1·ALDH2 생존 유지 ✅.**
  - 산출: `results/step4_A_repro/SVM-RFE.gene.v2.txt`, `interGenes.v2.txt` (v1 은 보존).
- caret rfe 의 `method="svmRadial"`: 원본 오타 `methods=` → 올바른 `method=` 로 수정.

### 논문 ROC-AUC 대조 시도 결과 (docx/보충표)
- `Figures S1_S10.docx` 는 **그림 범례 텍스트만** 포함, **AUC 수치 없음**(수치는 그림 이미지 안).
- 발견: 보충 ROC 그림 패널은 **CDKN1B, TSPYL5, VNN2, XAF1** 4개 유전자용 → **FN1/ALDH2 가 아님**.
  이 4개는 우리 A후보(63)·LASSO·SVM 어디에도 없음(train/test 엔 존재). 즉 **논문 ML 진단 패널은 우리 MR-필터
  경로와 다른 유전자 세트**이며, FN1/ALDH2 는 multi-omics 하이라이트로 본문 그림에서 개별 ROC 제시된 것으로 보임.
  참고: 우리 데이터에서 이 4개 AUC(test)는 낮음(CDKN1B 0.49, TSPYL5 0.57, VNN2 0.69, XAF1 0.54).
- 어떤 보충표에도 AUC 열 없음 → **논문 FN1/ALDH2 AUC 수치는 오프라인 figshare 로 추출 불가**(본문 그림 필요).

### 위험·한계
- ⚠️ **ALDH2 검증 AUC(0.70) 가 훈련(0.85) 대비 낮음** — 세뇨관(test) 조직 차이/GSE 내부 플랫폼 배치효과
  (미해결 #2)/Control 정의(LD vs LD+TN, 미해결 #1) 영향 가능. **(사용자 지시로 지금은 원인 규명 보류, 한계로만 기록.)**
- 논문 ML 패널(CDKN1B/TSPYL5/VNN2/XAF1) 과 우리 패널({ALDH2,CA2,FN1,ZFAND5}) 이 다름 — 후보 유전자
  선정 경로(DEG 정의·MR DB·pathway 유니버스) 차이로 추정. FN1/ALDH2 재현이라는 목표엔 부합.

---

## STEP 3-4 재실행 — 논문 본문(Figure 4) 확정 설계 (2026-07-07) — `10_step34_paper_design.R`

### 설계 (이전 STEP3/4 와 다름, 이쪽이 논문 정본)
- **훈련 = GSE96804 단독**(사구체 microarray). **검증 = ComBat(GSE104948 + GSE104954)** = 81샘플(Control 52/DKD 29).
  **GSE30529 는 ML 검증에서 제외.** (검증 매트릭스 `processed/data.valid.paper.txt` 저장.)
- **후보 = 저자 확정 10개** (Supp Table 8/9/10/12 "candidate genes" 컬럼에서 직접 확보 — 재구성 불필요):
  **ALDH2, FN1, VNN2, CREB5, XAF1, CA2, CDKN1B, IFI44L, SYTL2, TSPYL5.**
  = (MR risk ∩ DKD 고발현 DEG) ∪ (MR protective ∩ DKD 저발현 DEG) 후 MR 견고성(이질성/다면발현/Steiger) 통과분.
  방향 확인: GSE96804 DEG 부호와 MR b 부호가 6개 모두 일치(FN1↑risk, ALDH2/CDKN1B/TSPYL5/VNN2↓prot, XAF1↑risk).
  후보 목록 `results/step4_paper/interGenes.List.txt` 저장.

### 3) LASSO(10) → 6개
- 우리 결과: **{ALDH2, CDKN1B, FN1, IFI44L, VNN2, XAF1}** (6개).
- 논문 6개: {CDKN1B, ALDH2, FN1, XAF1, TSPYL5, VNN2}. → **5/6 일치** (우리 IFI44L ↔ 논문 TSPYL5 1개만 교체).
  단일 AUC 로 보면 TSPYL5(tr0.924)가 IFI44L(tr0.627)보다 우수 → seed/glmnet 버전/λ 차이로 경계에서 뒤바뀜.
  **핵심(FN1·ALDH2 선택)은 재현.**

### 4) 단일유전자 ROC-AUC (train=GSE96804 / valid=104948+104954) vs 논문
| 유전자 | 우리 train | 논문 train | 우리 valid | 논문 valid |
|---|---|---|---|---|
| **FN1**  | 0.909 | 0.911 | **0.915** | 0.911 | ← 거의 정확 재현
| **ALDH2**| 0.940 | 0.912 | 0.784 | 0.815 | ← train 약간↑, valid 약간↓ (근사)
- (참고) CA2 tr0.966/va0.912, VNN2 tr0.974/va0.710, CDKN1B tr0.967/va0.695, TSPYL5 tr0.924/va0.731.
- 전 STEP4(train=사구체병합/test=세뇨관)의 FN1 0.842/0.950·ALDH2 0.851/0.696 대비, **논문 설계(train=GSE96804 단독)로 바꾸니 FN1 이 논문값(0.911/0.911)에 훨씬 근접.**

### 5) FN1+ALDH2 결합 모델 (train=GSE96804, valid=104948+104954) vs Figure 4H,I
| 모델 | AUC_train | AUC_valid | 논문(Fig4) |
|---|---|---|---|
| GLM | 0.980 | **0.826** | valid 0.942 |
| RF | 1.000 | 0.914 | — |
| SVM | 0.973 | 0.785 | — |
| XGBoost | 0.997 | 0.873 | — |
- ⚠️ **GLM 검증 0.826 < 논문 0.942.** 단일 FN1(valid 0.915)보다도 낮음 → 훈련 계수가 검증(사구체+세뇨관 혼합)로
  잘 전이 안 됨(ALDH2 valid 0.784). 원인 후보: 검증 전처리/ComBat, 조직 이질성, 모델 학습 세부(교차검증/스케일).
  RF(0.914)·XGBoost(0.873)는 상대적으로 나음. **결합모델 정합은 부분적 — 추가 조정 여지(한계로 기록).**

### 남은 한계
- ALDH2·결합 GLM 의 검증 성능 갭 → (미해결 #1 Control 정의, #2 GSE 내부 플랫폼 배치) 및 검증 전처리와 연관 가능.
  사용자 지시로 원인 심층 규명은 보류, 한계로만 기록.

### 논문/원본 정렬 패스 (2026-07-09) — 원본 ML 코드 1:1 반영 + LD-only 재실행
- **입력 정리**: STEP4 가 검증셋을 자체 재생성하던 것을 **STEP3 산출물 `data.valid.paper.txt`(71) 를 읽도록** 변경(단일 소스).
  훈련=GSE96804(61). LD-only 정렬로 검증 81→71.
- **원본 ML 코드 반영(누락분 추가)**: 기존 paper 스크립트에 없던 **SVM-RFE + Venn 교집합 + ci.auc(bootstrap) CI** 를
  원본 `machine learning modeling 1·2.R` 그대로 추가.
  - LASSO: 원본과 **완전일치**(glmnet+cv.glmnet, deviance, nfolds10, lambda.min, seed123) → 6개 {ALDH2,CDKN1B,FN1,IFI44L,VNN2,XAF1}.
  - SVM-RFE: 원본과 동일(caretFuncs/cv/sizes2-8/center·scale), 원본 오타 `methods=`→`method=` 정렬 → {CA2,CDKN1B,VNN2}.
  - **교집합(SVM-RFE∩LASSO) = {CDKN1B,VNN2} → FN1/ALDH2 미포함**(원본 코드 방식의 특성).
- **논문 vs 원본코드 차이(중요)**: 논문 본문은 *"LASSO→ROC AUC>0.8(train&valid 모두)→multivariate"* 방법 →
  우리 계산 **AUC>0.8 통과 = {ALDH2, FN1, CA2}** 로 **FN1·ALDH2 를 정확히 선택**. 즉 **FN1/ALDH2 는 논문의
  ROC-필터 방법에서 선택되며, 코드의 SVM-RFE 교집합에서는 아님**. 둘 다 산출하고 논문 방법을 정본으로 채택.
- **LD-only 재실행 효과**: ALDH2 검증 AUC **0.784→0.807**(논문 0.815 에 근접), FN1 검증 0.915→0.871,
  결합 GLM valid 0.826→0.820(RF 0.927/XGB 0.898). 단일유전자 CI 도 재계산.
- **중간산출물**(`results/step4_paper/`): interGenes.List.txt(후보10), LASSO.gene.txt, SVM-RFE.gene.txt,
  interGenes.txt(교집합), single_gene_ROC_withCI.csv(10유전자 AUC+95%CI), single_gene_ROC_AUC.csv,
  genes_AUC_over_0.8_both.txt, combined_model_AUC.csv, compare_paper_vs_ours.AUC.csv, ROC PDF 6종, ROC_figures_index.csv.

### STEP 4 최종 확정 (2026-07-07)
- **STEP 4 는 논문 본문 설계 재현으로 확정**: FN1 단일 ROC 재현(0.909/0.915 ≈ 논문 0.911/0.911),
  ALDH2 근사(0.940/0.784 ≈ 0.912/0.815), LASSO 5/6 일치.
- **아래 2건은 "파지 않는 알려진 한계"로 확정 기록**(추가 조정 안 함):
  1) FN1+ALDH2 결합 GLM 검증 0.826 < 논문 0.942 (RF 0.914/XGB 0.873 는 근접).
  2) LASSO 6개 중 IFI44L ↔ 논문 TSPYL5 1개 교체(λ 선택/표준화/glmnet 버전 경계 차이).

## STEP 7 — 단일세포 RNA-seq (2026-07-08, 완료) — `14_step7_scrna.R`
원본 `scRNA.Seurat0.7.R` 로직(Seurat 표준). 출력 `results/step7_scrna/`.

### 데이터·방법
- **GSE209781**(10x, 6샘플: NM01-03=Control, DKD01-03=DKD) — 원본 스크립트가 쓴 바로 그 데이터셋.
  (GSE131882=dgecounts rds, GSE266146=BD Rhapsody zip 은 포맷 상이 → 이번엔 제외, 필요 시 추가.)
- Read10X → 병합(총 132k 세포) → QC(nFeature 300~5000, mt<10%) → **18,817 세포** →
  LogNormalize → HVG(2000) → ScaleData → PCA(30) → **Harmony 통합(orig.ident)** →
  FindClusters(res 0.5, 21 클러스터) → UMAP → **마커기반 세포주석**(클러스터별 표준 신장 마커세트 평균발현 최대).
- 세포유형 11종: PCT, LOH, Endothelial, Mesangial, Mono_Mac, T_cell, NK, B_cell, Plasma, Mast, Neutrophil.
- 겪은 버그(수정): (1) Seurat v5 AverageExpression 이 숫자 클러스터에 'g' 접두사 → 제거. (2) 이름붙은 벡터를
  meta 대입 시 바코드 매칭 오류 → unname(). 무거운 전처리는 체크포인트(rds)로 1회만.

### 핵심 결과 — FN1/ALDH2 세포유형별 발현 (mean, %expressing)
| 유전자 | 최고발현 세포유형 | mean | %expr | 2위 |
|---|---|---|---|---|
| **FN1** | **Endothelial(내피)** | 1.721 | 73.5% | Mesangial(系膜) 0.929 / 53.4% |
| **ALDH2** | **PCT(근위세뇨관)** | 1.648 | 78.1% | LOH 0.609 / 40.1% |
- 나머지 세포유형은 FN1<0.12 / ALDH2<0.51 로 낮음.
- **해석**: FN1 = 사구체 내피·系膜세포(ECM/섬유화 구조세포), ALDH2 = 근위세뇨관(대사 활성세포).
  **STEP 5(FN1=EMT/ECM, ALDH2=산화적인산화/대사)와 완전 정합**, 논문 scRNA 그림(내피_系膜_FN1 / PCT_ALDH2)과 일치.

### 위험·한계
- 마커기반 자동주석(argmax) — Podocyte/DCT/CD 는 별도 군집으로 안 잡힘(저발현/병합). FN1/ALDH2 국소화 결론엔 무관.
- 논문은 손상PCT(dPCT)에서도 FN1 상향 언급 — 본 주석은 PCT 단일 군집으로 병합. 세분화는 추가 작업 필요.
- GSE131882/266146 미처리(포맷 상이). 필요 시 별도 진행.

### 논문/원본 정렬 패스 (2026-07-09) — 원본 scRNA.Seurat0.7.R 반영
- **해상도 정렬**: `FindClusters(resolution 0.5→0.2)` 로 원본과 동일. → 18,817 세포, **15 클러스터**.
- **Harmony 변수**: 원본은 `RunHarmony("patient")`(Control/DKD 2군), 우리는 `orig.ident`(6샘플) 유지 —
  샘플단위 배치보정(표준, 원본 group단위보다 세분). FN1/ALDH2 국소화 결론 무관 → DIFF "의도적 차이" 기록.
- **중간산출물 추가**: `cluster_marker_scores.csv`(주석 근거 점수행렬), `cluster_markers_top20.csv`(FindAllMarkers
  min.pct0.25/logfc0.5/only.pos = 원본 인자), `compare_paper_vs_ours.scRNA.csv`.
- **결과(res 0.2)**: **FN1 = Endothelial(mean 1.713, 73.8%)**, **ALDH2 = PCT(mean 1.533, 75.9%)** →
  논문(내피/系膜_FN1, PCT_ALDH2) 정합. res 0.5 결과(FN1 1.721/ALDH2 1.648)와 거의 동일 — 결론 견고.

---

## STEP 6 — 2-표본 MR (2026-07-08, 완료) — `13_step6_mr.R`
원본 `Mendelian randomization 1~3.R` 로직. 출력 `results/step6_mr/`.

### 데이터·방법
- 노출 = FN1(ENSG00000115414)·ALDH2(ENSG00000111275) cis-eQTL (eqtl-a GWAS-VCF, FORMAT ES:SE:LP:AF:SS:ID).
- **FN1 vcf(2-1, 비압축)는 원본 미변경 — scratchpad 복사본을 Rsamtools 로 bgzip+tabix 재압축 성공** 후 사용.
- ⚠️ 오프라인이라 LD clumping API 불가 → **저자 도구변수(Supp Table 8, 이미 clump됨)를 노출 SNP 집합으로 사용**,
  exposure 통계(beta/se/allele/eaf)는 vcf 에서 추출. vcf 매칭: FN1 3/3, ALDH2 15/15 전부 존재(파싱 검증됨).
- 결과 = FinnGen R12 DKD(3-1, 비압축 TSV 2.1GB) 먼저 → GWAS Catalog GCST90435706(3-2, gz) 추가.
- harmonise_data → mr(IVW, MR-Egger, weighted median) → 이질성·다면발현 검정.

### 결과 — 저자 Supp9 와 사실상 정확 재현
| 유전자 | 결과 | 우리 IVW | 저자 Supp9 IVW |
|---|---|---|---|
| **FN1** | GCST | OR 2.777, b=1.021, p=0.0122, nsnp 3 | b=1.021, p=0.0122, nsnp 3 |
| **ALDH2** | GCST | OR 0.673, b=−0.395, p=0.0321, nsnp 14 | b=−0.395, p=0.0321, nsnp 14 |
- **b·se·p·nsnp 전부 일치** → MR 완전 재현.
- **FinnGen outcome**: FN1 OR 1.025 p=0.785, ALDH2 OR 0.973 p=0.459 → **둘 다 비유의**.
  이는 저자 Supp6(FinnGen MR)에 FN1/ALDH2 가 없던 것과 정합(이들은 GWAS Catalog 에서만 유의).
- 다면발현(MR-Egger 절편) 통과: FN1 p=0.99, ALDH2 p=0.70 (수평 다면발현 없음).
- **결론 재현**: **FN1 = DKD 위험인과(OR>1), ALDH2 = 보호인과(OR<1)** — 논문 핵심 MR 주장 재현.

### 위험·한계
- LD clumping 을 저자 도구변수(Supp8)로 대체(오프라인 제약). 독립 clumping 은 미수행 → 저자 instrument 에 의존.
  단 결과가 Supp9 와 정확히 일치하므로 재현 목적상 타당.
- FinnGen 은 전체 파일 fread(2.1GB, 필요열만) — 메모리 사용 큼(1회성, 정상 완료).

### 논문/원본 정렬 패스 (2026-07-09) — 원본 MR1·2.R 누락분 반영
- **outcome 필터 추가**: 원본 `Mendelian randomization 1.R` 의 `dat[dat$pval.outcome>5e-06,]`(역인과 방지)를 추가(NA 가드).
  적용 후에도 **GCST FN1 nsnp3 / ALDH2 nsnp14 = Supp9 유지**.
- **IVW 필터 추가(원본 MR2.R)**: IVW p<0.05 & 3방법 OR방향 일치 & 다면발현 p>0.05 →
  **FinnGen 통과 0**(FN1/ALDH2 비유의 = Supp6 정합), **GCST 통과 {ALDH2, FN1}** → `IVW.filter.GCST.csv`.
- **민감도 추가(원본 MR1.R)**: `mr_singlesnp` + `mr_leaveoneout` 표 저장(무거운 PDF plot 은 생략).
- **Supp9 대조 CSV**: `compare_paper_vs_ours.MR_IVW.csv` — FN1 b=1.0212/nsnp3, ALDH2 b=−0.3954/nsnp14 **소수점 일치**.
- 다면발현: FN1 절편 p=0.99, ALDH2 p=0.70(수평 다면발현 없음).
- 중간산출물(`results/step6_mr/`): exposure, table.SNP/MRresult/heterogeneity/pleiotropy/singleSNP/leaveoneout × FinnGen·GCST,
  IVW.filter.GCST.csv, compare_paper_vs_ours.MR_IVW.csv.

---

### STEP 4 ROC 그림 (2026-07-07) — `12_step4_roc_plots.R`
- 논문 Figure 4 스타일 pROC 그림 6개 → `results/step4_paper/`:
  단일 `train/valid_ROC.{FN1,ALDH2}.pdf`(AUC+95%CI bootstrap 라벨, 빨강+대각선),
  결합 `combined_ROC_{train,valid}.pdf`(GLM/RF/SVM/XGBoost 4곡선 범례).
- CSV 대조 전부 일치(그림 AUC ≈ single_gene_ROC_AUC.csv / combined_model_AUC.csv, tol 0.02).
  단일유전자 그림 AUC 는 `ci.auc(bootstrap)` 점추정이라 CSV 대비 ±0.002 수준 미세차(무시 가능).
  95%CI: FN1 tr 0.821–0.974 / va 0.828–0.977, ALDH2 tr 0.878–0.985 / va 0.668–0.887.

---

## STEP 5 — Hallmark GSEA + KEGG ssGSEA (2026-07-07, 완료) — `11_step5_gsea.R`
원본 `hallmark.gsea.R` 로직(clusterProfiler::GSEA, logFC 랭킹). 제공 gmt 사용(원본 msigdbr 대체).
출력 `results/step5_gsea/`.

### (1) Hallmark GSEA — GSE142025 DEG 3비교 (유의 p.adj<0.05 개수: Early_vs_Ctrl 26 / Late_vs_Ctrl 29 / Late_vs_Early 32)
FN1/ALDH2 관련 경로 방향·유의성:
- **EMT (EPITHELIAL_MESENCHYMAL_TRANSITION)**: Late_vs_Control **NES +2.52 (p.adj 1e-9)**,
  Late_vs_Early **NES +2.59 (p.adj 5e-10)** → DKD 진행에서 강한 상향. **FN1(EMT 유전자)와 정합.**
- 염증(INFLAMMATORY_RESPONSE): Late 비교에서 강한 상향(NES +2.2~+2.7, p.adj<1e-8).
- **산화적 인산화(OXIDATIVE_PHOSPHORYLATION)**: Late_vs_Control **NES -1.62 (p.adj 4.6e-4)**,
  Late_vs_Early **NES -2.02 (p.adj 7.7e-9)** → DKD 진행에서 하향. **ALDH2(미토콘드리아 효소) 보호역할과 정합.**
- ⚠️ Hallmark **REACTIVE_OXYGEN_SPECIES 경로 자체는 비유의**(p.adj 0.22~0.99). "산화" 신호는 ROS 경로가 아니라
  산화적 인산화(미토콘드리아 대사)로 나타남 — 정직하게 기록.

### (2) KEGG ssGSEA — GSE96804(train), 186 경로 x 61 샘플
- DKD vs Control(Wilcoxon): 상위 유의 경로 대부분 **대사경로 하향**(adipocytokine/insulin signaling, pentose phosphate,
  glycerolipid, 여러 아미노산·당대사; deltaDKD<0) → DKD 에서 대사 억제.
- **FN1 ↔ KEGG_ECM_RECEPTOR_INTERACTION rho=0.919** (FN1 최상위 상관) → **FN1=fibronectin, ECM/섬유화 축과 정합.**
- **ALDH2 ↔ 대사경로 강상관**: TRYPTOPHAN(0.936)/ARGININE_PROLINE(0.933)/TCA CYCLE(0.921)/여러 아미노산 대사 rho~0.91–0.94
  → **ALDH2=미토콘드리아 대사 효소 역할과 정합.**

### 요약
- **FN1 축 = EMT/ECM(섬유화)**, **ALDH2 축 = 미토콘드리아 산화적 인산화/대사** — 둘 다 유의하게 재현.
- 원본 대비 차이: gmt 파일 사용(msigdbr 아님), GSEA p<1e-10 은 eps 기본값이라 하한 클리핑(정밀도 경고, 결론 무관).

### 논문/원본 정렬 패스 (2026-07-09) — 원본 GSEA.R/hallmark.gsea.R 1:1 반영
- **Hallmark GSEA 입력 정렬**: 원본 `hallmark.gsea.R` 는 **diff_(유의 DEG)** logFC 를 GSEA 함 → 우리도 all_→**diff_ 로 변경**.
  결과 EMT +2.78/+2.53(논문 +2.52/+2.59), OXPHOS -1.71/-1.82(-1.62/-2.02) — 부호·크기 근접 유지.
- **KEGG GSEA 추가(누락분)**: 원본 `GSEA.R`(KEGG GSEA on diff_ DEG)을 그대로 구현 →
  Late_vs_Early ECM_RECEPTOR +1.95(p.adj 9e-4), TCA -1.58/OXPHOS -1.48 하향. (이전엔 ssGSEA 만 있었음.)
- **ssGSEA 유지**: 원본 GSEA.R 은 KEGG **GSEA** 이고 ssGSEA 스크립트는 미제공 → ssGSEA(GSVA)+FN1/ALDH2 상관은
  **논문 Fig2H,I 방법**으로 유지(원본 코드 대응 없음). FN1↔ECM 0.919, ALDH2↔TCA 0.921 재현.
- **중간산출물**(`results/step5_gsea/`): Hallmark.GSEA.*(3), KEGG.GSEA.*(3), Hallmark_focus/KEGG.GSEA_focus,
  compare_paper_vs_ours.Hallmark_NES.csv, KEGG.ssGSEA.*(scores/DKDvsControl/FN1_ALDH2_corr).
