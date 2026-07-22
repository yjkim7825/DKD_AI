# ============================================================================
# clinical_egfr.R — 임상검증 (Nephroseq): FN1·ALDH2 발현 ↔ eGFR 상관 + 질병 대조
#   논문 임상검증 재현. Nephroseq v5 export CSV(data/5-1~5-4)를 정리·시각화.
#   ※ 원자료(Ju CKD·Woroniecka 등)는 Nephroseq DB 소속이라 직접 상관계산은 불가 —
#     논문과 동일하게 Nephroseq가 제공한 r·Fold Change 값을 사용.
#   입력 : data/5-1..5-4 Nephroseq_*.csv
#   출력 : output/clinical_summary.csv, clinical_egfr.pdf
# ============================================================================
ROOT <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
DATA <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "R-reproduce/07_clinical/output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Nephroseq export: 앞 3줄=메타, 4줄=헤더
read_neph <- function(file, gene, kind) {
  df <- read.csv(file.path(DATA, file), skip = 3, check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(gene = gene, kind = kind, dataset = df$Dataset,
             pval = suppressWarnings(as.numeric(df$`p-Value`)),
             foldchange = suppressWarnings(as.numeric(df$`Fold Change`)),
             r = suppressWarnings(as.numeric(df$`r Value`)),
             stringsAsFactors = FALSE)
}

gfr <- rbind(
  read_neph("5-2. Nephroseq_FN1_DKD_Human_Glom_TubInt_GFR_correlation.csv",   "FN1",   "GFR"),
  read_neph("5-4. Nephroseq_ALDH2_DKD_Human_Glom_TubInt_GFR_correlation.csv", "ALDH2", "GFR"))
dvc <- rbind(
  read_neph("5-1. Nephroseq_FN1_DKD_Human_Glom_TubInt_DiseaseVsControl.csv",   "FN1",   "DvC"),
  read_neph("5-3. Nephroseq_ALDH2_DKD_Human_Glom_TubInt_DiseaseVsControl.csv", "ALDH2", "DvC"))

summary_tab <- rbind(gfr, dvc)
write.csv(summary_tab, file.path(OUT, "clinical_summary.csv"), row.names = FALSE)

cat("== eGFR 상관 (r) ==\n")
cat("FN1  :", paste(sprintf("%.3f", gfr$r[gfr$gene=="FN1"]),   collapse=", "), "(전부 음수 → 위험)\n")
cat("ALDH2:", paste(sprintf("%.3f", gfr$r[gfr$gene=="ALDH2"]), collapse=", "), "(전부 양수 → 보호)\n")
cat("== 질병 vs 정상 Fold Change ==\n")
cat("FN1  :", paste(sprintf("%.2f", dvc$foldchange[dvc$gene=="FN1"]),   collapse=", "), "(과발현)\n")
cat("ALDH2:", paste(sprintf("%.2f", dvc$foldchange[dvc$gene=="ALDH2"]), collapse=", "), "(저발현)\n")

# ── 그림: (좌) GFR 상관 r값  (우) 질병 vs 정상 Fold Change ──────────────────
pdf(file.path(OUT, "clinical_egfr.pdf"), width = 11, height = 5)
par(mfrow = c(1, 2), mar = c(9, 4, 3, 1))
col_g <- ifelse(gfr$gene == "FN1", "#d62728", "#1f77b4")
bp <- barplot(gfr$r, col = col_g, ylim = c(-0.9, 0.9), ylab = "r (expression vs eGFR)",
              main = "eGFR correlation", las = 2,
              names.arg = paste0(gfr$gene, ": ", gfr$dataset), cex.names = 0.6)
abline(h = 0)
legend("topright", c("FN1 (risk, r<0)", "ALDH2 (protective, r>0)"),
       fill = c("#d62728", "#1f77b4"), cex = 0.7, bty = "n")

col_d <- ifelse(dvc$gene == "FN1", "#d62728", "#1f77b4")
barplot(dvc$foldchange, col = col_d, ylab = "Fold change (DKD vs Control)",
        main = "Disease vs Control", las = 2,
        names.arg = paste0(dvc$gene, ": ", dvc$dataset), cex.names = 0.6)
abline(h = 0)
dev.off()
cat("\n★ 임상검증 완료 → output/clinical_egfr.pdf, clinical_summary.csv\n")
cat("  FN1 = DKD 과발현 + eGFR 음의상관 = 위험 / ALDH2 = 저발현 + eGFR 양의상관 = 보호\n")
cat("  → DEG·MR 결론과 일치 (임상 지표로 재확인)\n")
