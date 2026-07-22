# ============================================================================
# mr_func.R — MR1 공용 함수 (노출 파싱 · 결과 읽기 · MR 실행)
#   저자 Mendelian randomization 1.R 의 노출/harmonise/mr/민감도 로직
#   ※ 필터(MR2)·forest(MR3) 는 02/03 스크립트에 별도 (저자처럼 파일 분리)
#   ※ 함수 정의만 — 실행은 01_mr_run.R
# ============================================================================
suppressMessages({ library(data.table); library(TwoSampleMR) })

# ── eQTL VCF 파싱 → exposure 후보 (SNP별 beta/se/eaf) ────────────────────────
#   eqtl-a GWAS-VCF: FORMAT(ES:SE:LP:AF:SS) 에서 노출 통계 추출
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

# ── outcome(GWAS) 읽기 (필요 SNP만) ─────────────────────────────────────────
read_outcome <- function(file, cols, map, name, snps) {
  dt <- fread(file, select = cols, header = TRUE, showProgress = FALSE)
  setnames(dt, names(map), unname(map))
  dt <- dt[SNP %in% snps]
  message("[outcome ", name, "] 매칭 SNP ", nrow(dt), " / ", length(snps))
  od <- format_data(as.data.frame(dt), type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    eaf_col = "eaf", pval_col = "pval")
  od$outcome <- name; od
}

# ── MR 실행 (저자 MR1.R): harmonise → mr → 민감도 ───────────────────────────
run_mr <- function(exposure_dat, outcome_dat, tag, outDir) {
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  message("\n===== MR: ", tag, " =====")
  dat <- harmonise_data(exposure_dat, outcome_dat)
  # 저자 MR1.R: outcome 유의 SNP(p<=5e-06) 제거 (역인과/다면발현 방지)
  n0 <- nrow(dat)
  dat <- dat[!is.na(dat$pval.outcome) & dat$pval.outcome > 5e-06, ]
  message("[", tag, "] pval.outcome>5e-06 필터: ", n0, " -> ", nrow(dat), " SNP")
  write.csv(dat[dat$mr_keep == TRUE, ], file.path(outDir, paste0("table.SNP.", tag, ".csv")), row.names = FALSE)

  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  mrTab <- generate_odds_ratios(res)
  write.csv(mrTab, file.path(outDir, paste0("table.MRresult.", tag, ".csv")), row.names = FALSE)

  # 민감도 (저자 MR1.R)
  het <- tryCatch(mr_heterogeneity(dat),   error = function(e) NULL)
  ple <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
  ss  <- tryCatch(mr_singlesnp(dat),       error = function(e) NULL)
  loo <- tryCatch(mr_leaveoneout(dat),     error = function(e) NULL)
  # 역인과(방향) 검정 — MR Steiger (논문 Supp12). 유전자별 방향 판정.
  #   DKD는 이진 형질 → r.outcome을 get_r_from_lor로 먼저 계산해야 Steiger 작동.
  st <- tryCatch({
    d2 <- dat
    # 결과(DKD)는 로그오즈 → r 계산 (유병률 대략치 prevalence=0.1 가정). 노출은 정량으로 근사.
    if (!"r.outcome" %in% names(d2) || all(is.na(d2$r.outcome))) {
      ncase <- 5000; ncontrol <- 300000              # 대략치(방향 판정엔 스케일 무관)
      d2$r.outcome <- get_r_from_lor(lor = d2$beta.outcome, af = d2$eaf.exposure,
                                     ncase = ncase, ncontrol = ncontrol,
                                     prevalence = 0.1, model = "logit")
    }
    directionality_test(d2)
  }, error = function(e) { message("[", tag, "] Steiger 실패: ", conditionMessage(e)); NULL })
  if (!is.null(het)) write.csv(het, file.path(outDir, paste0("table.heterogeneity.", tag, ".csv")), row.names = FALSE)
  if (!is.null(ple)) write.csv(ple, file.path(outDir, paste0("table.pleiotropy.", tag, ".csv")), row.names = FALSE)
  if (!is.null(ss))  write.csv(ss,  file.path(outDir, paste0("table.singleSNP.", tag, ".csv")), row.names = FALSE)
  if (!is.null(loo)) write.csv(loo, file.path(outDir, paste0("table.leaveoneout.", tag, ".csv")), row.names = FALSE)
  if (!is.null(st))  write.csv(st,  file.path(outDir, paste0("table.steiger.", tag, ".csv")), row.names = FALSE)

  ivw <- mrTab[mrTab$method == "Inverse variance weighted", c("exposure","nsnp","or","or_lci95","or_uci95","pval")]
  cat("\n[", tag, "] IVW:\n"); print(ivw, row.names = FALSE)
  list(mrTab = mrTab, ple = ple, dat = dat)
}
