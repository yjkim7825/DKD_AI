# ============================================================================
# ml_func.R — 05_ml 공용 함수 (데이터 로딩 · ROC)
#   ※ 함수 정의만 — 실행은 01_ml_feature.R·02_ml_model.R / ml_author.R
# ============================================================================
suppressMessages({ library(pROC) })

# 발현행렬(행=유전자, 열=샘플) 읽어 후보 유전자만, 샘플×유전자로 전치 + 라벨 추출
#   반환: list(X = 샘플×유전자 행렬, y = factor(Control/DKD), genes = 유전자명)
load_expr <- function(exprFile, interGenes) {
  rt <- read.table(exprFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  # 후보목록 순서로 정렬 (학습·검증 파일 간 행 순서가 달라도 동일하게 맞춤)
  keep <- interGenes[interGenes %in% rownames(rt)]
  rt <- rt[keep, , drop = FALSE]
  X  <- t(as.matrix(rt))                                   # 샘플 × 유전자
  y  <- gsub(".*_(Control|DKD)$", "\\1", rownames(X))      # 열이름 끝 라벨
  list(X = X, y = factor(y, levels = c("Control", "DKD")), genes = rownames(rt))
}

# 학습셋 평균·표준편차로 학습·검증 둘 다 표준화 (배치는 ComBat이 이미 제거 → 학습기준 표준화가 정석·무누수)
scale_train_test <- function(Xtr, Xte) {
  mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[sdv == 0] <- 1e-8
  list(tr = scale(Xtr, center = mu, scale = sdv),
       te = scale(Xte, center = mu, scale = sdv))
}

# ── 논문 방식: 학습+검증을 합쳐 sva::ComBat 배치보정 후 후보 유전자만 분리 ──────
#   전처리(prep3)는 검증 내부(104948↔104954)만 보정 → 학습 vs 검증은 미보정.
#   따라서 여기 batch = "학습 vs 검증" (검증은 이미 내부 harmonize됨 → 재분할 안 함).
#   use_mod=TRUE → mod=model.matrix(~그룹)로 Control/DKD 생물학적 차이 명시 보존.
#     (단, 검증 라벨을 보정에 쓰게 되므로 성능 부풀림 논란 있음 → 기본 FALSE)
load_and_combat <- function(trainFile, validFile, interGenes, use_mod = FALSE) {
  suppressMessages(library(sva))
  rt <- read.table(trainFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  rv <- read.table(validFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  common <- intersect(rownames(rt), rownames(rv))
  comb <- cbind(rt[common, , drop = FALSE], rv[common, , drop = FALSE])
  samples <- colnames(comb)
  isTr  <- samples %in% colnames(rt)
  batch <- ifelse(isTr, "train", "valid")                    # 학습 vs 검증 2배치
  grp   <- factor(gsub(".*_(Control|DKD)$", "\\1", samples), levels = c("Control", "DKD"))
  mod   <- if (use_mod) model.matrix(~ grp) else NULL
  cat("[ComBat] batch=학습/검증 (", sum(isTr), "/", sum(!isTr), ") · mod=",
      if (use_mod) "~그룹(생물보존)" else "NULL", "\n", sep = "")
  combat <- sva::ComBat(dat = as.matrix(comb), batch = batch, mod = mod, par.prior = TRUE)
  keep <- interGenes[interGenes %in% rownames(combat)]
  X <- t(combat[keep, , drop = FALSE])                       # 샘플 × 후보유전자 (보정된 값)
  y <- factor(gsub(".*_(Control|DKD)$", "\\1", rownames(X)), levels = c("Control", "DKD"))
  list(Xtr = X[isTr, , drop = FALSE], ytr = y[isTr],
       Xte = X[!isTr, , drop = FALSE], yte = y[!isTr], genes = keep,
       full = combat, batch = batch, isTr = isTr)             # full = PCA용
}

# 유전자 하나의 학습/검증 ROC-AUC (방향 무관: <0.5면 1-AUC)
gene_auc <- function(y, expr) {
  a <- as.numeric(pROC::auc(pROC::roc(y, as.numeric(expr), quiet = TRUE)))
  max(a, 1 - a)
}
