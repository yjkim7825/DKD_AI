# eQTLGen 원본 도구변수로 DKD MR 재현 — 결과 요약

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
