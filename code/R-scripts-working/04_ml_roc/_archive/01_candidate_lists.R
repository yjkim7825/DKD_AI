# 07_step4_candidate_lists.R ------------------------------------------------
# STEP 4 준비: 후보 유전자 목록 2종 생성 (결과 분리 저장, 재현성 확보)
#   A(재현) = Supp4 DEG(DKD vs Control) ∩ MR-union(Supp6 FinnGen ∪ Supp7 GWAScat) = 63개
#   B(탐색) = 전체 DEG union(Supp3/4/5), 사전 MR 필터 없음
# 저자 보충표(figshare)에서 직접 파생 → interGenes.List.txt (scan() 호환: 심볼 1줄씩)
# ---------------------------------------------------------------------------

library(here)
source(here::here("config.R"))

SUPP <- file.path(DATA_ROOT, "6. Article related data")
sf <- function(n) file.path(SUPP, n)
A_DIR <- file.path(RES_DIR, "step4_A_repro")
B_DIR <- file.path(RES_DIR, "step4_B_explore")
dir.create(A_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(B_DIR, showWarnings = FALSE, recursive = TRUE)

# MR exposure 유전자(유의) — 표 자체가 유의 유전자만 담고 있음
read_mr <- function(n) {
  x <- read.csv(sf(n), header = FALSE, stringsAsFactors = FALSE, skip = 1)
  colnames(x) <- as.character(x[1, ]); x <- x[-1, ]
  unique(trimws(x$exposure))
}
mr6 <- read_mr("Supplementary Table 6. DKD-associated genes identified through two-sample MR analysis using FinnGen database summary statistics..csv")
mr7 <- read_mr("Supplementary Table 7. DKD-associated genes identified through two-sample MR analysis using GWAS Catalog summary statistics..csv")
mrUnion <- unique(c(mr6, mr7))

# DEG 유전자 (id 열)
read_deg <- function(n) {
  x <- read.csv(sf(n), header = TRUE, stringsAsFactors = FALSE, skip = 1)
  unique(trimws(x$id))
}
s3 <- read_deg("Supplementary Table 3.  Differentially expressed genes between late and early DKD patients.csv")
s4 <- read_deg("Supplementary Table 4. Differentially expressed genes between late DKD patients and controls .csv")
s5 <- read_deg("Supplementary Table 5. Differentially expressed genes between early DKD patients and controls.csv")
degUnion <- unique(c(s3, s4, s5))

# 훈련셋에 존재하는 유전자(피처 가용성 확인용)
train <- read.table(file.path(OUT_DIR, "data.train.txt"), header = TRUE, sep = "\t",
                    check.names = FALSE, row.names = 1)
tg <- rownames(train)

# A 목록: Supp4 ∩ MR-union
A_list <- sort(intersect(s4, mrUnion))
# B 목록: DEG union(3/4/5)
B_list <- sort(degUnion)

save_list <- function(genes, dir, provenance) {
  writeLines(genes, file.path(dir, "interGenes.List.txt"))
  writeLines(provenance, file.path(dir, "candidate_list.provenance.txt"))
}
save_list(A_list, A_DIR, c(
  "STEP4-A 재현용 후보 목록",
  sprintf("정의: Supp Table 4 DEG(Late DKD vs Control) ∩ MR-union(Supp6 FinnGen ∪ Supp7 GWAScat)"),
  sprintf("총 %d genes | data.train 존재 %d | FN1=%s ALDH2=%s",
          length(A_list), length(intersect(A_list, tg)), "FN1" %in% A_list, "ALDH2" %in% A_list)))
save_list(B_list, B_DIR, c(
  "STEP4-B 탐색용 후보 목록 (사전 MR 필터 없음)",
  sprintf("정의: 전체 DEG union(Supp3 Late-Early ∪ Supp4 Late-Ctrl ∪ Supp5 Early-Ctrl)"),
  sprintf("총 %d genes | data.train 존재 %d | FN1=%s ALDH2=%s",
          length(B_list), length(intersect(B_list, tg)), "FN1" %in% B_list, "ALDH2" %in% B_list)))

cat("== A(재현) ==\n")
cat("  후보:", length(A_list), " | train 존재:", length(intersect(A_list, tg)),
    " | FN1:", "FN1" %in% A_list, " ALDH2:", "ALDH2" %in% A_list, "\n")
cat("== B(탐색) ==\n")
cat("  후보:", length(B_list), " | train 존재:", length(intersect(B_list, tg)),
    " | FN1:", "FN1" %in% B_list, " ALDH2:", "ALDH2" %in% B_list, "\n")
cat("저장 ->", A_DIR, "/ ", B_DIR, "\n")
