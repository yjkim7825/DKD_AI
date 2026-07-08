# 13_step6_mr.R -------------------------------------------------------------
# STEP 6 : 2-표본 MR (원본 Mendelian randomization 1~3.R 로직)
#   노출(exposure) = FN1·ALDH2 cis-eQTL (eqtl-a GWAS-VCF)
#   결과(outcome)  = FinnGen R12 DKD (3-1) 먼저 → GWAS Catalog GCST90435706 (3-2) 추가
#   주분석 = IVW (+ MR-Egger, weighted median), 이질성/다면발현 검정.
#   ※ FN1 vcf(2-1) 는 원본 미변경 — scratchpad 복사본을 bgzip+tabix 재압축해 사용.
#   ※ 오프라인이라 LD clumping API 불가 → 저자 도구변수(Supp Table 8, 이미 clump됨)를
#     노출 SNP 집합으로 사용하되 exposure 통계(beta/se/allele/eaf)는 vcf 에서 추출(검증 포함).
#   출력: results/step6_mr/
# ---------------------------------------------------------------------------

suppressMessages({ library(data.table); library(TwoSampleMR); library(Rsamtools) })
library(here)
source(here::here("config.R"))
OUT <- file.path(RES_DIR, "step6_mr"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SCR <- file.path(SCRATCH_DIR, "mr")
dir.create(SCR, showWarnings = FALSE, recursive = TRUE)

FN1_VCF   <- file.path(DATA_ROOT, "2-1. eqtl-a-ENSG00000115414.vcf")
ALDH2_VCF <- file.path(DATA_ROOT, "2-2. eqtl-a-ENSG00000111275.vcf.gz")
FINNGEN   <- file.path(DATA_ROOT, "3-1. finngen_R12_DM_NEPHROPATHY_EXMORE")
GCST      <- file.path(DATA_ROOT, "3-2. GCST90435706", "GCST90435706.tsv.gz")
SUPP8     <- file.path(DATA_ROOT, "6. Article related data",
  "Supplementary Table 8. SNP characteristics and MR instrumental variables for candidate genes.csv")

## ---- 0) FN1 vcf 복사본 bgzip + tabix (원본 미변경) ----
fn1_copy <- file.path(SCR, "FN1.eqtl.vcf")
file.copy(FN1_VCF, fn1_copy, overwrite = TRUE)
fn1_bgz <- tryCatch({
  bg <- Rsamtools::bgzip(fn1_copy, dest = paste0(fn1_copy, ".gz"), overwrite = TRUE)
  idx <- tryCatch(Rsamtools::indexTabix(bg, format = "vcf"),
                  error = function(e) { message("[bgzip] tabix 인덱싱 실패(정렬 이슈 가능): ", conditionMessage(e)); NA })
  message("[bgzip] FN1 복사본 bgzip 완료: ", bg, if(!is.na(idx)) paste0(" | tabix: ", idx) else " | tabix 생략")
  bg
}, error = function(e) { message("[bgzip] 실패: ", conditionMessage(e)); NA })

## ---- 1) eQTL vcf 파싱 -> exposure 후보 (전체) ----
parse_vcf <- function(path, gene) {
  v <- fread(path, skip = "#CHROM", sep = "\t", header = TRUE)
  setnames(v, 1, "CHROM")
  smp <- names(v)[ncol(v)]
  fmt <- strsplit(v$FORMAT[1], ":")[[1]]
  p <- tstrsplit(v[[smp]], ":", fixed = TRUE); names(p) <- fmt
  data.frame(SNP = p$ID, chr.exposure = v$CHROM, pos.exposure = v$POS,
             effect_allele.exposure = v$ALT, other_allele.exposure = v$REF,
             beta.exposure = as.numeric(p$ES), se.exposure = as.numeric(p$SE),
             pval.exposure = 10^(-as.numeric(p$LP)), eaf.exposure = as.numeric(p$AF),
             samplesize.exposure = as.numeric(p$SS), exposure = gene,
             id.exposure = gene, stringsAsFactors = FALSE)
}
fn1 <- parse_vcf(if (!is.na(fn1_bgz)) fn1_bgz else FN1_VCF, "FN1")
aldh2 <- parse_vcf(ALDH2_VCF, "ALDH2")
vcfAll <- rbind(fn1, aldh2)

## ---- 2) 저자 도구변수(Supp8) ∩ vcf → 노출 SNP 집합 ----
s8 <- fread(SUPP8, skip = 1, header = TRUE)
authSNP <- split(s8$SNP, s8$exposure)[c("FN1","ALDH2")]
expo <- do.call(rbind, lapply(names(authSNP), function(g) {
  sub <- vcfAll[vcfAll$exposure == g & vcfAll$SNP %in% authSNP[[g]], ]
  message("[exposure ", g, "] 저자 SNP ", length(authSNP[[g]]), " 중 vcf 매칭 ", nrow(sub))
  sub
}))
# 노출 포맷
exposure_dat <- format_data(expo, type = "exposure",
  snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure", pval_col = "pval.exposure", phenotype_col = "exposure",
  samplesize_col = "samplesize.exposure", chr_col = "chr.exposure", pos_col = "pos.exposure")
write.csv(exposure_dat, file.path(OUT, "exposure.FN1_ALDH2.csv"), row.names = FALSE)

## ---- 3) outcome 읽기 (필요 SNP만) ----
snps <- unique(exposure_dat$SNP)
read_outcome <- function(file, gz, cols, map, name) {
  dt <- fread(file, select = cols, header = TRUE, showProgress = FALSE)
  setnames(dt, names(map), unname(map))     # 표준 열명으로
  dt <- dt[SNP %in% snps]
  message("[outcome ", name, "] 매칭 SNP ", nrow(dt), " / ", length(snps))
  od <- format_data(as.data.frame(dt), type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval")
  od$outcome <- name
  od
}
finn <- read_outcome(FINNGEN, FALSE,
  cols = c("rsids","beta","sebeta","alt","ref","pval","af_alt"),
  map = c(rsids="SNP", beta="beta", sebeta="se", alt="effect_allele", ref="other_allele", pval="pval", af_alt="eaf"),
  name = "FinnGen_R12_DKD")
gcst <- read_outcome(GCST, TRUE,
  cols = c("variant_id","beta","standard_error","effect_allele","other_allele","p_value","effect_allele_frequency"),
  map = c(variant_id="SNP", beta="beta", standard_error="se", effect_allele="effect_allele",
          other_allele="other_allele", p_value="pval", effect_allele_frequency="eaf"),
  name = "GWAScatalog_GCST90435706")

## ---- 4) harmonise + MR ----
run_mr <- function(outcome_dat, tag) {
  message("\n===== MR: ", tag, " =====")
  dat <- harmonise_data(exposure_dat, outcome_dat)
  outTab <- dat[dat$mr_keep == TRUE, ]
  write.csv(outTab, file.path(OUT, paste0("table.SNP.", tag, ".csv")), row.names = FALSE)
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  mrTab <- generate_odds_ratios(res)
  write.csv(mrTab, file.path(OUT, paste0("table.MRresult.", tag, ".csv")), row.names = FALSE)
  het <- tryCatch(mr_heterogeneity(dat), error = function(e) NULL)
  ple <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  if(!is.null(het)) write.csv(het, file.path(OUT, paste0("table.heterogeneity.", tag, ".csv")), row.names = FALSE)
  if(!is.null(ple)) write.csv(ple, file.path(OUT, paste0("table.pleiotropy.", tag, ".csv")), row.names = FALSE)
  cat("\n[", tag, "] IVW 결과:\n")
  ivw <- mrTab[mrTab$method == "Inverse variance weighted", c("exposure","nsnp","or","or_lci95","or_uci95","pval")]
  print(ivw, row.names = FALSE)
  invisible(mrTab)
}
mr_finn <- run_mr(finn, "FinnGen")
mr_gcst <- run_mr(gcst, "GCST")

## ---- 5) 저자 MR(Supp9 IVW)와 대조 ----
s9 <- fread(file.path(DATA_ROOT, "6. Article related data",
  "Supplementary Table 9. Mendelian randomization results of candidate genes IVW, MR-Egger, and weighted median analyses.csv"),
  skip = 1, header = TRUE)
s9$b <- as.numeric(s9$b)
s9ivw <- s9[s9$method == "Inverse variance weighted" & s9$exposure %in% c("FN1","ALDH2"),
            c("exposure","nsnp","b","se","pval")]
cat("\n== 저자 Supp9 IVW (GWAS Catalog outcome) ==\n"); print(s9ivw, row.names = FALSE)
message("\n[STEP6] 완료 -> ", OUT)
