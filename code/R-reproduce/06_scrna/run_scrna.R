# ============================================================================
# run_scrna.R — 06_scrna 원샷 실행 (데이터 압축해제 → 01→02→03→04)
#   실행: setwd(".../R-reproduce/06_scrna"); source("run_scrna.R")
#   ※ 무거움. 03(monocle)·04(CellChat)은 패키지 없으면 건너뛰고 경고만.
# ============================================================================
ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MLDIR <- file.path(ROOT, "R-reproduce/06_scrna")
RAW   <- file.path(ROOT, "data/1-7. GSE209781_RAW")

## ── 0) 데이터 압축해제 (tar.gz → 10x 폴더). 이미 있으면 건너뜀 ───────────────
need <- c("NM01","NM02","NM03","DKD01","DKD02","DKD03")
if (!all(dir.exists(file.path(RAW, need)))) {
  cat("[압축해제] GSE209781 tar.gz 6개 푸는 중...\n")
  for (f in list.files(RAW, pattern = "\\.tar\\.gz$", full.names = TRUE)) untar(f, exdir = RAW)
}
cat("[데이터] 10x 폴더: ", paste(basename(list.dirs(RAW, recursive = FALSE)), collapse = ", "), "\n")

## ── 패키지 있으면 실행, 없으면 경고 후 계속 ─────────────────────────────────
step <- function(title, file, pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  cat("\n##########", title, "##########\n")
  if (length(miss)) { cat("  [건너뜀] 패키지 없음:", paste(miss, collapse = ", "),
                          "\n  설치:  install.packages(c(", paste0('\"', miss, '\"', collapse = ","), "))\n"); return(invisible()) }
  source(file.path(MLDIR, file))
}

step("01. QC → Harmony → 클러스터링 → 세포주석", "01_scrna_process.R", c("Seurat","harmony"))
step("02. 세포비율 + FN1/ALDH2 + 유전자세트", "02_scrna_proportion.R", c("Seurat","ggpubr"))
step("03. Monocle 의사시간 (PCT→dPCT)", "03_scrna_pseudotime.R", c("monocle"))
step("04. CellChat 세포통신 (Control vs Late DKD)", "04_scrna_cellchat.R", c("CellChat"))

cat("\n★ 06_scrna 완료 — 결과는 output/ 폴더\n")
cat("  (03 monocle·04 CellChat 은 해당 패키지 설치 후 다시 source 하면 채워짐)\n")
