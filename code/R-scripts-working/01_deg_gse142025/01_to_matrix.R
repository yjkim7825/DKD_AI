# 01_GSE142025_to_matrix.R -------------------------------------------------
# RNA-seq(GSE142025) 샘플별 txt(Symbol + 값) 를 하나의 발현행렬로 병합 + 그룹 라벨링.
# 입력 : DIR_GSE142025 의 GSM*_{샘플}.txt.gz  (2열: Symbol, 값)
# 출력 : OUT_DIR/GSE142025.labeled.txt  (열 = '{샘플}_{그룹}')
# 그룹 : 파일명 접두사  N->Control(9), B->Early(6), A->Late(21)
# ※ 원본은 읽기만, 결과는 processed/ 로만 저장.

library(here)
source(here::here("config.R"))

files <- list.files(DIR_GSE142025, pattern = "\\.txt\\.gz$", full.names = TRUE)
stopifnot(length(files) > 0)

read_one <- function(f) {
  s <- sub(".*_(.+)\\.txt\\.gz$", "\\1", basename(f))   # 예: 'A11A'
  d <- read.delim(gzfile(f), header = TRUE, check.names = FALSE)
  colnames(d) <- c("Symbol", s)
  d
}
lst <- lapply(files, read_one)

# Symbol 기준 병합(교집합 유전자)
mat <- Reduce(function(x, y) merge(x, y, by = "Symbol", all = FALSE), lst)
rownames(mat) <- mat$Symbol
mat$Symbol <- NULL

# 그룹 라벨(접두사)
prefix <- substr(colnames(mat), 1, 1)
grp <- ifelse(prefix == "N", GROUP_CONTROL,
        ifelse(prefix == "B", GROUP_EARLY,
        ifelse(prefix == "A", GROUP_LATE, NA)))
if (any(is.na(grp))) stop("그룹 미지정 샘플 존재: ", paste(colnames(mat)[is.na(grp)], collapse=", "))
colnames(mat) <- paste0(colnames(mat), "_", grp)

out <- cbind(geneNames = rownames(mat), mat)
write.table(out, file.path(OUT_DIR, "GSE142025.labeled.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("[GSE142025] ", nrow(mat), " genes x ", ncol(mat), " samples\n", sep="")
print(table(grp))
cat("저장: ", file.path(OUT_DIR, "GSE142025.labeled.txt"), "\n")
