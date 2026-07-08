# CLAUDE.md — 작업 인수인계 (R-scripts-working)

## 0. 한 줄 요약
P25 논문(*FN1·ALDH2 DKD multi-omics*, Renal Failure 2025)의 R 파이프라인을 **로컬에서 실제로
돌려보는 중**입니다. **데이터는 이제 전부 다운로드 완료**(`../data/` 참고). 완벽 재현보다
"한 단계씩 실제로 돌아가게" 우선. **STEP 1(GSE142025 DEG) 완료 → 현재 STEP 2(CEL→RMA) 진행.**

### 진행 현황
- ✅ STEP 1 완료: GSE142025 DEG 3종 실행. Late_vs_Early = 3,314개(논문 2,833과 방향·비율 일치,
  유전자 70% 겹침). 차이 원인(입력 매트릭스/정규화, RNA-seq는 voom 고려)은 **미해결·나중 검토**로
  `DECISIONS.md` 에 기록하고 넘어감. 이 차이는 GSEA(STEP 5)에서 다시 보며 FN1/ALDH2 ROC엔 영향 없음.
- ▶ STEP 2 진행: CEL→RMA 발현매트릭스.

## 1. 프로젝트 배경
- 목표: bulk 전사체 → DEG → LASSO/SVM 바이오마커 → ROC (핵심), 이후 GSEA·MR·scRNA 확장.
- **원본 스크립트**(수정 금지, 참고용): `../R-scripts-for-pipeline-reproducibility/`
  - 경로가 `G:\187geneMR\...` 로 하드코딩돼 있고, "이미 만들어진 발현매트릭스"부터 시작함.
- **이 폴더**(작업 복제본, 여기서 수정): `R-scripts-working/`
- 참고 문서: `../../pipeline_port/docs/`
  - `01_필요데이터_목록.md` (확정 데이터 목록), `03_R코드_흐름분석.md` (각 단계 이유),
    `02_Python_재구성_가이드.md`. 그리고 `../../pipeline_port/` 에 Python 이식본도 있음.

## 2. 절대 규칙 (반드시 지킬 것)
1. **원본 데이터(`../data/` 의 번호 폴더) 절대 수정/삭제 금지.** 모든 산출물은
   `../data/processed/` (중간 매트릭스) 와 `./results/` (분석결과) 로만 쓴다.
2. **원본 스크립트 폴더(`../R-scripts-for-pipeline-reproducibility/`) 수정 금지.**
3. 데이터 파일을 고쳐야 하면(예: vcf 재압축) **원본 말고 복사본**을 만들어 수정한다.
4. 경로 하드코딩 금지 — 항상 `00_config_local.R` 를 `source()` 해서 경로를 참조.
5. R 실행은 이 폴더에서: `setwd(".../R-scripts-working")` 또는 `Rscript 파일.R`.

## 3. 운영 방식 — 멀티 에이전트 + Y/N 게이트
아래 3가지 관점(가능하면 서브에이전트로 분담)으로 진행:
- **[Porter]** 원본 R 로직을 이 폴더의 실행 가능한 스크립트로 구현/수정.
- **[Reviewer]** 원본 R(`../R-scripts-for-pipeline-reproducibility/`)과 1:1 대조 검토
  (대비 방향, 필터 임계, 정규화 방식 등 일치 확인).
- **[Tester]** `Rscript` 로 실제 실행하고 산출물/개수/에러를 검증.

**한 단계가 끝날 때마다 반드시:**
1) 아래 형식으로 **요약 보고**를 먼저 출력
   - 무엇을 했는가 / 원본 대비 달라진 점 / 실행·검증 결과(수치) / 위험·한계
2) 그 다음 **"다음 단계로 진행할까요? (Y/N)"** 로 묻고 **사용자 답을 기다린다.**
   사용자가 Y 하기 전에는 **다음 단계 파일을 만들거나 실행하지 않는다.**
- **한 번에 한 단계만.** 에러가 나면 다음으로 넘어가지 말고 그 단계에서 고친다.

## 4. 데이터 현황 (`../data/`, 번호 폴더 = 원본 읽기전용) — **전부 다운로드 완료**
| 경로 | 내용 | 상태 |
|---|---|---|
| `1-1. GSE96804_RAW` | 훈련 bulk, Affy CEL(61) + .chp | ✅ RAW |
| `1-2. GSE104948_RAW` | 검증 bulk(사구체), CEL(196) + **brainarray ENTREZG CDF** | ✅ RAW |
| `1-3. GSE104954_RAW` | 검증 bulk(세뇨관), CEL(195) + CDF | ✅ RAW |
| `1-4. GSE142025_RAW` | 초기/후기 DKD, RNA-seq txt(36) | ✅ (STEP1 처리 완료 → `../data/processed/GSE142025.labeled.txt`) |
| `1-5. GSE30529_RAW` | DEG 검증, Affy CEL(22, U133A_2) | ✅ RAW |
| `1-6. GSE131882_RAW` / `1-7. GSE209781_RAW` / `1-8. GSE266146_RAW` | scRNA(신장/신장/소변) | ✅ RAW |
| `2-1. eqtl-...115414` (FN1) | eQTL vcf | ✅ 있음 (⚠️ **비압축 `.vcf`, tbi는 `.vcf.gz`용 → STEP6에서 복사본 bgzip+tabix**) |
| `2-2. eqtl-...111275` (ALDH2) | eQTL vcf.gz + tbi | ✅ |
| `3-1. finngen_R12_DM_NEPHROPATHY_EXMORE` | MR outcome(FinnGen, 2.1GB) | ✅ |
| `3-2. GCST90435706` | MR outcome(GWAS Catalog, 903MB) | ✅ |
| `4-1. h.all.v2026.1.Hs.symbols.gmt` | Hallmark gene sets(50) | ✅ |
| `4-2. c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt` | KEGG gene sets(≈186) | ✅ |
| `5-1. Nephroseq_FN1_...all_analyses.csv` / `5-2. ..._GFR_correlation.csv` | FN1 발현/eGFR 상관 | ✅ (⚠️ **ALDH2 Nephroseq 없음** — 필요 시 별도 쿼리) |
| `6. Article related data` | figshare 번들: 보충표 S1~S13, `hallmark.gsea.R`, `scRNA.Seurat0.7.R`, GSEA 결과 | ✅ |

- 그룹 라벨: GSE142025 = 파일명 접두사 N=Control(9)/B=Early(6)/A=Late(21).
  GSE104948/54 = 파일명에 질환 인코딩(DN=DKD, LD/TN=정상대조). GSE96804/30529 = GEO 메타데이터 필요(GEOquery::getGEO).

## 5. 이미 만들어진 것 (이 폴더)
- `00_config_local.R` — 경로/파라미터. **모든 스크립트가 이걸 source.**
- `01_GSE142025_to_matrix.R` — GSE142025 txt 36개 → `processed/GSE142025.labeled.txt`
- `02_GSE142025_DEG.R` — limma DEG 3종(Late vs Early / Late vs Control / Early vs Control)

## 6. 해야 할 일 (순서 — 위에서부터)
- ✅ **STEP 1. GSE142025 DEG** — 완료. (`01_`,`02_` 실행됨. 차이 `DECISIONS.md` 기록)
- ▶ **STEP 2. CEL → 발현매트릭스 (RMA)** ← 지금 여기
  - GSE96804(oligo, Gene ST), GSE30529(affy, U133A_2), GSE104948/104954(affy + 제공된 brainarray ENTREZG CDF).
  - 그룹 라벨: 파일명(104948/54 = DN=DKD, LD/TN=control) or `GEOquery::getGEO`(96804/30529). 출력 `processed/*.labeled.txt`.
  - 필요 패키지(이 STEP에서 설치): `oligo`, `affy`, `GEOquery`, GSE96804용 pd.* 주석 패키지.
- **STEP 3. 병합 + ComBat** → 훈련(GSE96804)/검증(104948·104954) 세트 (원본 `data preprocessing 3.R` 참고).
- **STEP 4. LASSO + SVM-RFE 교집합 → ROC** (원본 `machine learning modeling 1·2.R`) : FN1/ALDH2 재현 (핵심).
- **STEP 5. GSEA/Hallmark** (figshare `hallmark.gsea.R` 참고). gmt 준비됨: `../data/4-1..gmt`(Hallmark), `../data/4-2..gmt`(KEGG). STEP1 DEG 차이도 여기서 재검토.
- **STEP 6. MR** : FN1 vcf **복사본** bgzip+tabix 먼저 → TwoSampleMR (원본 `Mendelian randomization 1~3.R`).
  outcome 둘 다 준비됨: `../data/3-1..`(FinnGen), `../data/3-2. GCST90435706`(GWAS Catalog). FinnGen 먼저 → GCST 추가.
- **STEP 7. scRNA** (figshare `scRNA.Seurat0.7.R` 참고) : GSE131882/209781/266146. (무겁고 느림 — 맨 마지막)

## 7. 환경
- Windows + R 설치됨. 패키지는 각 STEP 진입 시 필요한 것만 설치(무거운 Seurat/monocle/TwoSampleMR는 해당 STEP에서).
- `TwoSampleMR` 은 GitHub 설치일 수 있음: `remotes::install_github("MRCIEU/TwoSampleMR")`.

---
**시작 지시**: 위 규칙대로 [Porter/Reviewer/Tester] 관점으로 **STEP 1**부터. 실행·검증 후
요약 보고 → "다음 단계 진행? (Y/N)" 로 멈춰서 사용자 확인을 받을 것.
