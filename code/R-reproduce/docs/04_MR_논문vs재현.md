# 04. MR — 논문 vs 우리 재현

## 📂 파일 경로
- **코드**: `R-reproduce/04_mr/` (01_mr_run · 02_mr_filter · 03_mr_forest)
- **결과**: `04_mr/output/` (table.MRresult.*, forest_MR_*.png, compare_paper_vs_ours.MR_IVW.csv)
- **과정 노트북**: `00_explore/06_mr_preprocess.ipynb`(전처리), `05_mr_steps.ipynb`(IVW 계산)

---

## ★ 현재 재현 (04_mr, VCF 노출) — 논문 Supp9 완전 일치

우리 `04_mr` 구현은 저자 도구변수(Supp8, eqtl-a VCF)를 그대로 써서 **논문과 소수점까지 동일**.

| 유전자 | 우리(GCST IVW) | 논문 Supp9 | 판정 |
| --- | --- | --- | --- |
| **FN1** | OR 2.777, p=0.0122, nsnp 3 | OR 2.78, p=0.0122, nsnp 3 | ✅ 위험 (완전일치) |
| **ALDH2** | OR 0.673, p=0.0321, nsnp 14 | OR 0.673, p=0.0321, nsnp 14 | ✅ 보호 (완전일치) |

- **서사도 동일**: FinnGen(핀란드)선 비유의 → GCST(영국 38만)서 유의
- 3방법(IVW·Egger·median) 방향 일치, 민감도(이질성·다면발현) 이상 없음
- `compare_paper_vs_ours.MR_IVW.csv`: ours_b = paper_b (15자리 동일)

## ★ Figure 3D 완전 재현 — 인과 유전자 10개 (`forest_MR_FigureD.png`)

Supp8 도구변수 10개 전체를 두 코호트(FinnGen·GCST)에 MR → **IVW 필터 통과 유전자 = 논문 D의 10개와 정확히 일치**.
논문 D는 각 유전자를 **유의하게 나온 코호트** 결과로 보여준 것 (우리 그림은 Cohort 컬럼으로 이를 드러냄).

| 유전자 | nsnp | OR (95% CI) | p | 통과 코호트 | 방향 |
| --- | --- | --- | --- | --- | --- |
| ALDH2 | 14 | 0.673 (0.469–0.967) | 0.032 | GCST | 보호 |
| CA2 | 13 | 0.915 (0.848–0.987) | 0.022 | FinnGen | 보호 |
| CDKN1B | 3 | 0.788 (0.642–0.967) | 0.023 | FinnGen | 보호 |
| CREB5 | 7 | 1.067 (1.002–1.136) | 0.044 | FinnGen | 위험 |
| **FN1** | 3 | **2.777 (1.249–6.172)** | 0.012 | GCST | 위험 |
| IFI44L | 7 | 1.174 (1.030–1.339) | 0.016 | FinnGen | 위험 |
| SYTL2 | 6 | 1.090 (1.009–1.179) | 0.029 | FinnGen | 위험 |
| TSPYL5 | 4 | 0.897 (0.819–0.983) | 0.019 | FinnGen | 보호 |
| VNN2 | 7 | 0.769 (0.597–0.992) | 0.043 | GCST | 보호 |
| XAF1 | 5 | 1.219 (1.015–1.464) | 0.034 | FinnGen | 위험 |

- **10개 전부 OR·95%CI·p·nsnp까지 논문 D와 동일** ✅
- 필터 통과: **FinnGen 7개**(CA2·CDKN1B·CREB5·IFI44L·SYTL2·TSPYL5·XAF1) + **GCST 3개**(ALDH2·FN1·VNN2)
- → 저자 도구변수(Supp8) + 같은 방법으로 **인과 유전자 선별을 완전 재현**. (Venn A의 코호트별 유의 유전자 선별 로직도 이로써 확인)

---

# [부록] eQTLGen 원본 도구변수 버전 (탐색)

> 아래는 저자 완제품(VCF) 대신 **eQTLGen 원본에서 직접 도구변수를 만든** 별도 탐색 (노출 소스 차이로 OR 크기만 다름, 방향·유의는 동일). 위 04_mr(VCF)이 주 재현.

## 1. 목적

논문의 MR 핵심 서사 — **"FinnGen에서는 비유의, GWAS Catalog에서 인과 신호 발현, FN1=위험·ALDH2=보호"** — 를,
저자가 쓴 완제품(eqtl-a VCF)이 아니라 **eQTLGen 원본(2-3 cis-eQTL + 2-4 빈도)** 에서 직접 만든 도구변수로 재현할 수 있는지 검증.

## 2. 방법 (요약)

1. **노출 생성**: eQTLGen cis-eQTL(2-3, 3.7GB)을 청크 스트리밍으로 FN1(ENSG00000115414)·ALDH2(ENSG00000111275)만 필터 → 2-4 빈도로 Z-score를 beta/SE로 변환.
2. **LD 클럼핑**: 로컬 PLINK + 1000G EUR 참조패널, `p<5e-8, r²<0.001, kb=10000` (유전자별 독립) → FN1 3개 / ALDH2 3개 도구변수.
3. **2-표본 MR** (논문 `Mendelian randomization 1.R` 과 100% 동일 흐름):
   `read_exposure(clump=FALSE)` → outcome SNP 추출 → harmonise → `pval.outcome>5e-6` 필터 → IVW·MR-Egger·Weighted median → 이질성·다면발현·leave-one-out.
4. **두 outcome 각각 실행**: 3-1 FinnGen(핀란드 64,663) / 3-2 GCST90435706(영국 388,955).

FN1 도구변수 3개 = **rs10932612, rs17525860, rs615857** (= 논문 Supp9 FN1 IV와 동일).

## 3. 결과

### 3-1. 두 outcome IVW 결과

| Outcome | 유전자 | nsnp | OR | p | 판정 |
| --- | --- | --- | --- | --- | --- |
| FinnGen | FN1 | 3 | 1.04 | 0.792 | 비유의 |
| FinnGen | ALDH2 | 3 | 1.005 | 0.934 | 비유의 |
| **GCST** | **FN1** | 3 | **4.53** | **0.0121** | **유의 (위험, OR>1)** |
| GCST | ALDH2 | 3 | 0.572 | 0.097 | 보호 방향, 경계 비유의 |

GCST FN1은 세 방법 모두 일관: IVW OR 4.53(p=0.012), Weighted median OR 4.58(p=0.017), MR-Egger OR 4.34(p=0.63, 3-SNP라 저검정력).

### 3-2. 논문 Supp9 · 이전 세션(VCF) 대조 (outcome=GCST)

| 유전자 | 우리 (eQTLGen→GCST) | 논문 Supp9 | 이전 세션 (VCF→GCST) |
| --- | --- | --- | --- |
| FN1 | OR 4.53, p=0.0121, nsnp 3 | OR 2.78, p=0.0122, nsnp 3 | OR 2.78, p=0.0122, nsnp 3 |
| ALDH2 | OR 0.572, p=0.097, nsnp 3 | OR 0.673, p=0.0321, nsnp 14 | OR 0.673, p=0.0321, nsnp 14 |

- 이전 세션(VCF 노출)은 논문 Supp9와 **완전 일치** — 저자 IV와 동일했기 때문(재확인됨).
- 우리(eQTLGen 노출)는 **방향·유의성 재현**, 세부 수치는 노출 스케일·IV 수 차이로 다름(아래 해석).

### 3-3. 민감도 분석 (outcome=GCST, 모두 이상 없음)

| 유전자 | 이질성 Q_pval(IVW) | 다면발현 Egger 절편 p | 판정 |
| --- | --- | --- | --- |
| FN1 | 0.49 | 0.99 | 이질성·다면발현 없음 |
| ALDH2 | 0.49 | 0.45 | 이질성·다면발현 없음 |

(FinnGen도 동일하게 이질성·다면발현 없음: FN1 0.30/0.37, ALDH2 0.34/0.47)

## 4. 해석

**FN1 — 방향·유의성 정확 재현 (효과크기만 스케일 차이)**

- 같은 3개 SNP, p=0.0121 ≈ Supp9 0.0122 로 **유의성 일치**, 방향도 위험(OR>1)으로 일치.
- OR이 4.53 vs 2.78로 큰 이유는 **노출 단위 차이**: 우리 노출 beta는 eQTLGen Z-변환(SD 단위)이라 저자 VCF ES보다 약 1.4~1.5배 작음 → MR 비율(결과β/노출β)이 그만큼 커짐.
- ⚠️ 따라서 **OR 크기는 논문과 직접 비교 불가**(더 강한 효과가 아님). 재현된 것은 **방향 + 유의성**.

**ALDH2 — 인과 방향 재현, 유의성은 검정력 미달**

- 방향은 논문과 같은 **보호적(OR 0.57 < 1)**, 하지만 p=0.097로 경계 비유의(논문 0.032 유의).
- 원인은 순수 **검정력**: 엄격 클럼핑(r²<0.001)으로 eQTLGen cis 독립신호가 **3개**로 축약된 반면, 논문은 **14개** IV 사용. IV가 적으면 IVW 표준오차가 커져 유의성을 잃음. 효과 **방향 자체는 논문과 동일**.

## 5. 결론

| 유전자 | FinnGen (우리) | GCST (우리) | 논문 패턴 | 판정 |
| --- | --- | --- | --- | --- |
| FN1 | null (p=0.79) | 유의 위험 (OR 4.53, p=0.012) | FinnGen 비유의 → GCST 유의 위험 | **완전 일치 ✓** |
| ALDH2 | null (p=0.93) | 보호 방향, 경계 (OR 0.57, p=0.097) | 방향 일치, 유의성은 IV 수(3 vs 14) 차이 | **방향 일치 ✓** |

**종합**: 논문의 핵심 서사(FinnGen null → GCST 유의, FN1=위험·ALDH2=보호)를 **eQTLGen 원본에서 직접 만든 도구변수로 재현**함. FN1은 outcome별 패턴·유의성까지 정확히 일치, ALDH2는 인과 방향 일치(GCST 유의성은 엄격 클럼핑에 따른 IV 수 한계로 경계).

## 6. 한계 / 주의

- FN1 OR 크기는 노출 스케일 때문에 논문과 직접 비교 불가 — 보고 시 "방향·유의성 재현, 효과크기는 단위 차이" 명시 필요.
- ALDH2 GCST 유의성 미달은 오류가 아니라 IV 수(검정력) 문제.
- 3-SNP MR은 MR-Egger(다면발현 검정)의 신뢰도가 낮음 — 3-IV MR의 정상적 한계.

## 7. 산출물 위치

- 노출: `code/data_output/exposure.eQTLGen.FN1_ALDH2.clumped.csv`
- FinnGen MR: `code/data_output/eqtlgen_mr/` (표·플롯)
- GCST MR: `code/data_output/eqtlgen_mr_gcst/` (표·플롯)
- 참조패널: `code/data_output/ref/EUR.*`

## 8. 남은 선택지 (ALDH2 처리)

- **(A) 현재대로 "방향 일치, 유의성은 IV 수 한계"로 정직하게 보고** ← 추천(가장 방어 가능)
- (B) 클럼핑 완화(r²·kb 상향)로 IV를 늘려 유의성 회복 시도 (단, 도구변수 독립성 약해짐)
- (C) 논문의 14개 IV를 그대로 가져와 대조
