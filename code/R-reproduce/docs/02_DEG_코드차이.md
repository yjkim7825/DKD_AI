# 02. DEG — 저자 원본 코드 vs 우리 코드

저자 `differential expression analysis.R` ↔ 우리 `02_deg/R/deg_limma.R` + `run_deg.R` 비교.
**핵심 로직(limma)은 100% 동일**, 우리는 함수화 + 3그룹 확장만 함.

## 📂 파일 경로 (입력 → 결과·사진)

**코드**: `02_deg/R/deg_limma.R`(함수), `run_deg.R`(실행), `run_deg_3group_plots.R`(3그룹 그림)

**입력**: `01_preprocessing/bulk/output/02_normalized/GSE{142025,30529}.normalized.txt`

**결과·사진** — `02_deg/output/`
| 폴더 | 사진(pdf) | 표(txt) |
| --- | --- | --- |
| `GSE142025_3group/` | `heatmap_*.pdf`, `vol_*.pdf` (3비교 각각) | `diff_*`, `all_*`, `up_/down_*`, `TOP20_*` |
| `GSE142025_3group_combined/` | `heatmap_3group.pdf`, `vol_3stage.pdf`, `venn_stage.pdf` | `TOP_3stage.txt` |
| `GSE30529_2group/` | `heatmap_*.pdf`, `vol_*.pdf` (검증) | `diff_*` 등 |

- 3비교 = `Early_vs_Control`, `Late_vs_Control`, `Late_vs_Early`
- 과정 노트북: `00_explore/03_deg_steps.ipynb`

---

## 0. 한눈에

| 항목 | 저자 원본 | 우리 |
| --- | --- | --- |
| 파일 | `differential expression analysis.R` 1개 | `deg_limma.R`(함수) + `run_deg.R`(실행) 분리 |
| 비교 | **2그룹 1개** (DKD vs Control) | **함수 1개로 여러 비교** 재사용 |
| 대상 | GSE142025_**twoGroups** | GSE142025 **3그룹**(3비교) + GSE30529 검증 |
| 방법 | model.matrix→lmFit→makeContrasts→eBayes→topTable | **완전 동일** |
| 기준 | `logFCfilter=0.585`, `adj.P<0.05` | **동일** (`LOGFC=0.585`, `ADJP=0.05`) |
| 출력 | all/diff/up/down/heatmap/vol | **동일 + TOP20 테이블 추가** |

→ 알고리즘·기준·출력이 같고, **우리는 "함수 한 개 = 두 그룹 비교"로 만들어 3그룹·검증셋에 돌려쓴 것**만 다름.

---

## 1. 단계별 코드 대조

### ① 입력 읽기 + 중복 평균

**저자**
```r
rt = read.table(inputFile, header=T, sep="\t", check.names=F)
rt = as.matrix(rt); rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data = avereps(data)          # 같은 유전자 여러 행 → 평균
```

**우리** (`deg_limma.R` `read_expr`)
```r
rt <- read.table(path, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
mat <- as.matrix(rt)
grp <- sapply(strsplit(colnames(mat), "_"), function(x) tail(x, 1))
```
- 차이: 저자는 `avereps`를 **DEG 단계에서** 다시 함. 우리는 **전처리 단계(prep1)에서 이미 avereps** 해서 넘어옴 → 여기선 생략(값 동일).
- 그룹 라벨: 저자 `x[2]`(두 번째 조각), 우리 `tail(x,1)`(마지막 조각) — 열이름 형식에 맞춘 것뿐, 결과 같음.

### ② 그룹 필터 + 정렬

**저자**
```r
Type = sapply(strsplit(colnames(data), "_"), function(x) x[2])
keep_samples = Type %in% c("Control", "DKD")   # 두 그룹만
data = data[, keep_samples]; Type = Type[keep_samples]
data = data[, order(Type)]                      # 그룹 정렬
```

**우리** (`deg_two` 안)
```r
keep <- grp %in% c(caseGrp, ctrlGrp)            # 두 그룹만 (인자로 받음)
d <- mat[, keep]; g <- grp[keep]
d <- d[, order(g)]; g <- sort(g)
```
- **완전 동일 로직**. 저자는 "Control/DKD" 하드코딩, 우리는 **`caseGrp`/`ctrlGrp` 인자**로 받아 아무 두 그룹이나 비교 가능(그래서 Late/Early/Control 3비교에 재사용).

### ③ limma 차등발현 (★ 핵심 — 완전 동일)

**저자**
```r
design <- model.matrix(~0 + factor(Type))
colnames(design) <- levels(factor(Type))
fit <- lmFit(data, design)
cont.matrix <- makeContrasts(DKD - Control, levels=design)
fit2 <- eBayes(contrasts.fit(fit, cont.matrix))
allDiff = topTable(fit2, adjust='fdr', number=200000)
```

**우리**
```r
design <- model.matrix(~0 + factor(g)); colnames(design) <- levels(factor(g))
fit <- lmFit(d, design)
cont <- makeContrasts(contrasts = paste0(caseGrp, " - ", ctrlGrp), levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
```
- **똑같음**. 대비식만 저자는 `DKD - Control` 고정, 우리는 `paste0(caseGrp," - ",ctrlGrp)` 로 동적 생성.

### ④ 유의 DEG 필터

**저자**
```r
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
```
**우리**
```r
sig <- allDiff[abs(allDiff$logFC) > LOGFC & allDiff$adj.P.Val < ADJP, ]
```
- **동일** (`0.585` & `0.05`).

### ⑤ 저장 (all/diff/up/down)
- 저자·우리 모두 `all_`, `diff_`, `up_`, `down_` 파일 저장 — 형식 동일.
- **우리만 추가**: `TOP20_*.txt` (|logFC| 큰 순 상위20 — "가장 많이 변한 유전자" 빠르게 보려고).

### ⑥ 히트맵 + 볼케이노
- **완전 동일**: top50↑+50↓ 히트맵(`scale="row"`, blue2-white-red2, cluster_cols=FALSE), 볼케이노(Up/Down/Not 3색).
- 저자 코드 그대로 옮김.

---

## 2. 우리가 "더" 한 것 (3그룹 통합) — `run_deg_3group_plots.R`

저자는 2그룹만. 우리는 3그룹이라 **추가 시각화**를 만듦:

| 추가물 | 내용 |
| --- | --- |
| 3그룹 히트맵 | Control(초록)/Early(주황)/Late(빨강) 컬러바 한 번에, `gaps_col`로 그룹 구분 |
| 3분류 볼케이노 | ①초기에만(N→E) 주황 / ②둘다(계속) 초록 / ③말기에만(E→L) 빨강 |
| 벤 다이어그램 | N→E ∩ E→L 겹침 |
| TOP_3stage.txt | 구간별 최다 변화 유전자 |

→ **저자에 없는 "언제 변하는 유전자인가"(초기/계속/말기) 3분류**가 우리 확장 포인트.

---

## 3. 결과 대조 (개수)

| 비교 | 우리 | 논문 | 차이 |
| --- | --- | --- | --- |
| Late vs Early | 3,314 | 2,833 | +481 |
| Late vs Control | 4,022 | 3,525 | +497 |
| Early vs Control | 671 | 390 | +281 |
| GSE30529(검증) | 1,178 | — | 검증셋 |

- **패턴 일치**: Early<Late (진행할수록 DEG 급증), Late vs Control 최다.
- **개수는 우리가 ~300~500 많음** → 전처리(정규화·avereps 시점) + 임계값 경계 차이. 방향·핵심 유전자(FN1↑, IEG↓)는 동일.

**3분류 (우리 고유)**: ①초기에만 432 / ②둘다-계속 239 / ③말기에만 3,075
- 초기 대표(↓): FOSB, NR4A1/2/3, DUSP1 (IEG — 초기에 꺼짐)
- 계속 변함(↓): FOS, EGR1, ATF3 (IEG)
- 말기에만: 3,075개로 폭증 → **말기에 대규모 발현 붕괴**

---

## 4. 한 줄 요약
> 저자 limma 로직을 **그대로** 함수화(`deg_two`)해서 2그룹→3그룹·검증셋에 재사용. 알고리즘·기준·출력 동일, **3분류 시각화와 TOP20 테이블만 우리가 추가**. 개수는 전처리 차이로 논문보다 소폭 많지만 패턴·핵심 유전자는 일치.
