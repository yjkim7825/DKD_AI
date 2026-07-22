# ============================================================================
# prep3_merge_combat.R  — 전처리 3: 여러 발현행렬 병합 + ComBat 배치보정
#   · 공통 유전자(intersect)만 남겨 열 방향으로 병합
#   · batch = 데이터셋 번호 → ComBat 으로 배치효과 제거
#   · 1개만 넣으면 배치가 없으므로 ComBat 생략(그대로 반환)
#   ※ 함수 정의만 — 실행은 run_preprocess.R 에서
# ============================================================================
suppressMessages(library(sva))

# mats = 발현행렬 리스트 (각 행렬: 행=유전자, 열=샘플)
prep3_merge_combat <- function(mats) {
  genes  <- Reduce(intersect, lapply(mats, rownames))          # 모든 데이터 공통 유전자
  allTab <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE]))
  batch  <- rep(seq_along(mats), sapply(mats, ncol))           # 샘플별 배치(데이터셋) 번호
  if (length(mats) < 2) return(allTab)                         # 1개면 배치 없음 → 생략
  keep <- apply(allTab, 1, function(r) all(tapply(r, batch, function(x) sd(x) > 0)))  # 배치내 분산0 제거
  ComBat(dat = allTab[keep, , drop = FALSE], batch = batch, par.prior = TRUE)
}
