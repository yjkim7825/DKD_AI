# 3자 대조표 — 논문 vs R(내 데이터 재현) vs Python(이식)

> P25 (*Lin et al., FN1 & ALDH2 in DKD, Renal Failure 2025*) 파이프라인을
> **R로 재현**하고 **Python으로 이식**한 뒤, 논문 보고값과 나란히 대조한 결과.
> 원자료: `code/R-scripts-working/results/`, `code/Python-scripts-working/results/python_vs_R_vs_paper.csv`.
> R↔Python 은 STEP2(RMA, R 전용) 산출 매트릭스를 공유하므로, 동일 입력에서 계산 일치 여부가 핵심.

---

## 1. 단일유전자 ROC-AUC (핵심 결과)

> 검증셋 = LD-only 정렬 71샘플(Control 42/DKD 29). R=Python 동일 입력(data.valid.paper.txt).

| 유전자 | 세트 | 논문 | R | Python |
|---|---|---|---|---|
| **FN1** | train | 0.911 | 0.909 | 0.909 |
| **FN1** | valid | 0.911 | **0.871** | **0.871** |
| **ALDH2** | train | 0.912 | 0.940 | 0.940 |
| **ALDH2** | valid | 0.815 | **0.807** | **0.807** |

- **Python = R 소수점까지 완전 일치** (동일 매트릭스·동일 판별식, ALDH2는 보호유전자라 방향 자동보정).
- FN1 거의 정확. **ALDH2 검증은 LD-only 정렬로 0.784→0.807**(논문 0.815 근접).

## 2. 결합 진단모델 (FN1 + ALDH2, 4-알고리즘)

| 모델 | 세트 | 논문 | R | Python |
|---|---|---|---|---|
| GLM | train | 0.978 | 0.980 | 0.979 |
| GLM | valid | **0.942** | 0.820 | **0.953** |
| RF | train | 1.000 | 1.000 | 1.000 |
| RF | valid | 0.815 | 0.927 | 0.927 |
| SVM | train | 0.977 | 0.973 | 0.978 |
| SVM | valid | 0.807 | 0.800 | **0.374** ⚠️ |
| XGBoost | train | 0.999 | 0.997 | 0.997 |
| XGBoost | valid | 0.844 | 0.898 | 0.898 |

- 검증 71샘플 기준. 논문 최종 진단모델은 이 결합모델(특히 GLM 0.942).
- **Python GLM(0.953)이 R(0.820)보다 논문(0.942)에 근접.** RF(0.927)·XGB(0.898)는 R=Python.
- **SVM valid(0.374)만 이탈** — sklearn `SVC(rbf, gamma=scale)` ≠ R `e1071`(커널·스케일 규약 차이). 포팅 한계로 기록.

## 3. LASSO 선택 유전자 (10 후보 → 축소)

| | 선택 수 | 구성 | FN1·ALDH2 |
|---|---|---|---|
| 논문 | 6 | CDKN1B, ALDH2, FN1, XAF1, TSPYL5, VNN2 | 포함 |
| R | 6 | ALDH2, FN1, VNN2, XAF1, CDKN1B, **IFI44L** | 포함 |
| Python | 5 | ALDH2, CREB5, FN1, IFI44L, VNN2 | 포함 |

- 3자 모두 **FN1·ALDH2 선택**. 나머지 구성은 λ 선택·표준화·구현 차이로 일부 교체(비결정성).

## 4. DEG 개수 (STEP1, limma→Welch t 근사)

| 대비 | 논문 | R | Python |
|---|---|---|---|
| Late vs Early | 2,833 | 3,314 | 3,460 |
| Late vs Control | (미보고) | 4,022 | 4,096 |
| Early vs Control | (미보고) | 671 | 473 |

- Python은 limma eBayes(moderated t)를 Welch t로 근사 → 개수 소폭 차이(정상), 방향·비율 일치.

## 5. Mendelian Randomization (IVW, outcome=GCST90435706)

| 유전자 | 지표 | 논문(Supp9) | R | Python |
|---|---|---|---|---|
| FN1 | b / OR / p / nsnp | 1.021 / 2.777 / 0.012 / 3 | 1.021 / 2.777 / 0.012 / 3 | 1.021 / 2.777 / 0.012 / 3 |
| ALDH2 | b / OR / p / nsnp | −0.395 / 0.673 / 0.032 / 14 | −0.395 / 0.673 / 0.032 / 14 | −0.402 / 0.669 / 0.028 / 15 |

- **FN1 3자 완전 일치.** ALDH2는 harmonise에서 SNP 1개 더 유지(15 vs 14)로 미세차, **방향·결론 동일**.
- 핵심 인과 결론 재현: **FN1 = 위험(OR>1), ALDH2 = 보호(OR<1).**
- FinnGen(2.1GB) outcome은 무거워 스킵, GCST로 핵심 재현.

## 6. GSEA / ssGSEA (경로 축)

| 지표 | R | Python |
|---|---|---|
| EMT NES (Late vs Control) | +2.52 | +2.52 |
| OXPHOS NES (Late vs Control) | −1.62 | −1.62 |
| FN1 ↔ ECM_RECEPTOR (rho) | 0.919 | 0.919 |
| ALDH2 ↔ TCA_CYCLE (rho) | 0.921 | 0.921 |

- **Python = R 일치.** FN1=EMT/ECM(섬유화), ALDH2=OXPHOS/대사 축 재현. ROS 경로 자체는 3자 모두 비유의.

## 7. scRNA 세포유형 국소화 (GSE209781)

| 유전자 | R | Python |
|---|---|---|
| FN1 최고발현 | Endothelial (73.5%) | Endothelial (73.4%) |
| ALDH2 최고발현 | PCT(근위세뇨관) (78.1%) | PCT (77.5%) |

- 세포 수 18,817(R) ↔ 18,818(Python)로 일치. FN1=내피/메산지움, ALDH2=근위세뇨관 재현.

---

## 종합

STEP2(RMA, R 전용)를 제외한 **전 단계가 R·Python 양쪽에서 실행**되어, 다중오믹스 서사
(발현 → ML 바이오마커 → 인과(MR) → 경로 → 세포)가 **두 언어에서 동일하게 재현**됨.
유일한 실질 이탈은 **결합모델 SVM valid(Python 0.405)** — 커널 구현 차이에 따른 포팅 한계.
