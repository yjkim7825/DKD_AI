# 2. Python 재구성용 폴더 구조 + 터미널 Claude 운영 가이드

전체 파이프라인(bulk + MR + scRNA)을 Python으로 재구성할 때의 **폴더 구조**와,
**터미널 Claude(Claude Code)** 를 "멀티 에이전트처럼 역할 분담 + 단계별 보고 + yes/no 게이트"
방식으로 운영하는 **명령/프롬프트 템플릿**.

---

## A. 제안 폴더 구조

```
dkd_biomarker/
├── README.md
├── pyproject.toml            # 또는 requirements.txt (의존성)
├── config/
│   └── config.yaml           # 경로·파라미터(하드코딩 금지, 여기만 수정)
├── data/
│   ├── raw/                  # GEO/GWAS 원본 (git 제외)
│   ├── interim/              # 전처리 중간 산출
│   └── processed/            # 정규화·병합 최종 입력
├── src/dkd/
│   ├── io/                   # GEO/GPL 파서, eQTL/GWAS 로더
│   ├── bulk/                 # 전처리·정규화·ComBat·DEG·GSEA
│   ├── ml/                   # LASSO·SVM-RFE·ROC
│   ├── mr/                   # Mendelian randomization (또는 R 위임)
│   ├── scrna/                # scanpy 통합·주석·pseudotime·통신
│   └── utils/                # 공통(limma eBayes 등)
├── pipelines/                # 단계 실행 스크립트(01_~ 07_)
│   ├── 01_preprocess.py
│   ├── 02_deg.py
│   ├── 03_ml.py
│   ├── 04_roc.py
│   ├── 05_gsea.py
│   ├── 06_mr.py
│   └── 07_scrna.py
├── results/                  # 표·그림(단계별 하위폴더)
├── tests/                    # smoke/unit test (합성데이터)
└── docs/                     # 데이터 목록·흐름·의사결정 로그
```

> 지금 만든 `pipeline_port/python` 이 이 구조의 **bulk+ml 부분 초안**입니다.
> 확장 시 `src/dkd/` 로 모듈화하고 `pipelines/` 로 단계 실행을 분리하면 관리가 쉬움.

**설계 원칙 3가지**
1. **경로는 config 한 곳** (모든 하드코딩 제거) — 원본 R의 `setwd("G:\\...")` 반복 문제 해결.
2. **단계 = 파일 1개**, 입력→출력이 파일로 명확 → 중간부터 재실행 가능.
3. **tests/ 에 합성데이터 스모크테스트** → 실제 데이터 없이도 로직 검증(현 `smoke_test.py` 방식).

---

## B. 터미널 Claude(Claude Code) 운영 방식

### B-0. 공통 규칙 프롬프트 (세션 맨 처음 1번 붙여넣기)
```
너는 DKD 바이오마커 파이프라인을 R→Python으로 재구성하는 작업을 한다.
반드시 아래 운영 규칙을 지켜라:
1) 역할 분담(멀티 에이전트처럼):
   - [Porter] 원본 R 로직을 Python으로 이식(구현)
   - [Reviewer] Porter 결과를 원본 R과 대조 검토(수치/로직 일치, 엣지케이스)
   - [Tester] 합성데이터로 스모크테스트 작성·실행
   각 단계에서 세 역할의 관점을 각각 명시해서 진행하라.
2) 한 단계가 끝나면 반드시 아래 형식으로 '요약 보고'를 먼저 출력:
   - 무엇을 했는가 / 원본 대비 달라진 점 / 검증 결과 / 위험·한계
3) 그 다음 "다음 단계로 진행할까요? (yes/no)" 로 물어보고 내 답을 기다려라.
   내가 yes 하기 전에는 절대 다음 단계 파일을 만들지 마라.
4) 경로는 config 파일만 참조하고 절대 하드코딩하지 마라.
```

### B-1. 서브에이전트로 역할 분담 (Claude Code 방식)
Claude Code에서는 서브에이전트를 띄워 병렬/전문화할 수 있다. 예:
```
1단계(전처리)를 진행하자.
- Porter 서브에이전트: R의 data preprocessing 1~3.R 로직을 src/dkd/bulk/preprocess.py 로 이식
- Reviewer 서브에이전트: 이식본이 원본과 동일한지(quantile 정규화, ComBat 파라미터) 대조 리뷰
- Tester 서브에이전트: tests/test_preprocess.py 에 합성데이터 검증 추가
끝나면 세 결과를 통합 요약 보고하고 yes/no로 물어봐.
```

### B-2. 단계별 진행 명령 (순서대로, 각 단계 후 yes/no)
| 단계 | 붙여넣을 지시(요약) | 원본 R |
|---|---|---|
| 1 | "전처리 모듈 이식: probe→gene, log2/quantile, ComBat" | preprocessing 1~3 |
| 2 | "DEG 이식: limma 2군 moderated t-test (eBayes 그대로)" | differential expression |
| 3 | "GSEA 이식: gseapy.prerank(KEGG)" | GSEA |
| 4 | "특징선택 이식: LASSO(LogRegCV L1)+SVM-RFE(RFECV)+교집합" | ML modeling 1 |
| 5 | "ROC 이식: sklearn AUC + bootstrap CI, 훈련/검증" | ML modeling 2 |
| 6 | "MR: TwoSampleMR 상당(파이썬) 또는 R 위임 래핑" | MR 1~3 |
| 7 | "scRNA 이식: scanpy 통합·주석·pseudotime·통신(liana/cellchat)" | single-cell 1~4 |

각 단계 지시 끝에 항상 덧붙일 문장:
```
끝나면 요약 보고 후 "다음 단계 진행? (yes/no)"로 멈춰라.
```

### B-3. 검토(Reviewer) 강화용 지시 예시
```
방금 만든 DEG 이식본을 Reviewer 관점에서 원본 differential expression analysis.R와
1:1 대조해라. 특히 (a) 대비(contrast) 방향 DKD-Control, (b) adj.P가 BH인지,
(c) logFC 임계 0.585, (d) eBayes 분산 shrinkage 반영 여부를 표로 점검하고,
같은 입력에 대해 R과 Python 결과가 일치하는지 확인할 방법을 제시하라.
```

---

## C. 실무 팁
- **한 번에 한 단계만.** 7단계를 한꺼번에 시키면 검토가 얕아짐 → yes/no 게이트가 품질을 지킴.
- **R↔Python 대조**는 "같은 GEO 1개"로: R로 뽑은 DEG 표와 Python 표를 gene 기준 join → logFC/adjP 상관 확인(권장 r>0.99).
- **MR·scRNA는 R 유지도 정당한 선택.** 무리한 이식보다 `rpy2`/`subprocess`로 R 단계를 감싸고, 결과만 Python으로 이어받는 하이브리드가 안전.
- 의사결정은 `docs/decision_log.md`에 남겨(어느 단계에서 무엇을 왜 바꿨는지) → 논문 Methods 작성 때 그대로 재활용.
