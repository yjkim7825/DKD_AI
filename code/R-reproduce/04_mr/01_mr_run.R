# ============================================================================
# 01_mr_run.R  (= 저자 Mendelian randomization 1.R)
#   노출(eQTL) + 결과(GWAS) → harmonise → mr(IVW·Egger·median) → 민감도
#   [전 단계] 02_deg 후보 유전자(FN1·ALDH2)를 인과 검증
#   출력(→ 02·03 스크립트 입력): output/table.MRresult.*.csv, table.pleiotropy.*.csv 등
# ============================================================================
ROOT   <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MRDIR  <- file.path(ROOT, "R-reproduce/04_mr")
DATA   <- file.path(ROOT, "data")
DEGDIR <- file.path(ROOT, "R-reproduce/02_deg/output/GSE142025_3group")
OUT    <- file.path(MRDIR, "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(MRDIR, "R/mr_func.R"))

FN1_VCF   <- file.path(DATA, "2-1. eqtl-a-ENSG00000115414.vcf")
ALDH2_VCF <- file.path(DATA, "2-2. eqtl-a-ENSG00000111275.vcf.gz")
FINNGEN   <- file.path(DATA, "3-1. finngen_R12_DM_NEPHROPATHY_EXMORE.gz")
GCST      <- file.path(DATA, "3-2. GCST90435706", "GCST90435706.tsv.gz")
SUPP8     <- file.path(DATA, "6. Article related data",
  "Supplementary Table 8. SNP characteristics and MR instrumental variables for candidate genes.csv")
SUPP9     <- file.path(DATA, "6. Article related data",
  "Supplementary Table 9. Mendelian randomization results of candidate genes IVW, MR-Egger, and weighted median analyses.csv")

## [전 단계 확인] MR 대상이 DEG인지
deg <- read.table(file.path(DEGDIR, "diff_Late_vs_Control.txt"),
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
for (g in c("FN1","ALDH2"))
  cat("[DEG확인]", g, if (g %in% rownames(deg)) "→ 말기vs정상 DEG ✔" else "→ DEG 아님", "\n")

## 1) 노출: 저자 Supp8 = 후보 유전자 10개 도구변수 (이미 p<5e-8+LD clump, beta/se/eaf 있음)
#   (VCF 파싱값과 동일함이 검증됨 → Supp8 직접 사용, 10개 유전자 모두 → forest D 재현)
s8 <- data.table::fread(SUPP8, skip = 1, header = TRUE)
cat("[노출] Supp8 후보 유전자 ", length(unique(s8$exposure)), "개: ",
    paste(unique(s8$exposure), collapse = ", "), "\n", sep = "")
exposure_dat <- TwoSampleMR::format_data(as.data.frame(s8), type = "exposure",
  snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure", pval_col = "pval.exposure", phenotype_col = "exposure",
  samplesize_col = "samplesize.exposure", chr_col = "chr.exposure", pos_col = "pos.exposure")
write.csv(exposure_dat, file.path(OUT, "exposure.candidates10.csv"), row.names = FALSE)
snps <- unique(exposure_dat$SNP)

## 2) 결과(outcome) 2종
finn <- read_outcome(FINNGEN,
  cols = c("rsids","beta","sebeta","alt","ref","pval","af_alt"),
  map  = c(rsids="SNP", beta="beta", sebeta="se", alt="effect_allele", ref="other_allele", pval="pval", af_alt="eaf"),
  name = "FinnGen_R12_DKD", snps = snps)
gcst <- read_outcome(GCST,
  cols = c("variant_id","beta","standard_error","effect_allele","other_allele","p_value","effect_allele_frequency"),
  map  = c(variant_id="SNP", beta="beta", standard_error="se", effect_allele="effect_allele",
           other_allele="other_allele", p_value="pval", effect_allele_frequency="eaf"),
  name = "GWAScatalog_GCST90435706", snps = snps)

## 3) MR 실행 (결과별) → CSV 저장 (02·03이 읽음)
run_mr(exposure_dat, finn, "FinnGen", OUT)
mr_gcst <- run_mr(exposure_dat, gcst, "GCST", OUT)

## 4) 논문 Supp9(IVW) 대조 — 후보 유전자 10개 전체
s9 <- data.table::fread(SUPP9, skip = 1, header = TRUE); s9$b <- as.numeric(s9$b)
cand <- unique(exposure_dat$exposure)
s9ivw <- s9[s9$method == "Inverse variance weighted" & s9$exposure %in% cand, c("exposure","nsnp","b","se","pval")]
data.table::setnames(s9ivw, c("nsnp","b","se","pval"), c("paper_nsnp","paper_b","paper_se","paper_pval"))
ourIVW <- mr_gcst$mrTab[mr_gcst$mrTab$method == "Inverse variance weighted", c("exposure","nsnp","b","se","or","pval")]
data.table::setnames(ourIVW, c("nsnp","b","se","pval"), c("ours_nsnp","ours_b","ours_se","ours_pval"))
cmp <- merge(as.data.frame(s9ivw), as.data.frame(ourIVW), by = "exposure", all = TRUE)
write.csv(cmp, file.path(OUT, "compare_paper_vs_ours.MR_IVW.csv"), row.names = FALSE)
cat("\n== 논문 Supp9 vs 우리(GCST IVW) ==\n"); print(cmp, row.names = FALSE)
cat("\n★ [01] MR 실행 완료 → output/ (다음: 02_mr_filter.R)\n")
