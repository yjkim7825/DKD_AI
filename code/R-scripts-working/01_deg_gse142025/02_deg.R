# 02_GSE142025_DEG.R -------------------------------------------------------
# GSE142025 발현행렬로 limma 차등발현(DEG). 논문 Figure 2 / Table S3~S5 대응.
#   Late  vs Early    -> Table S3
#   Late  vs Control  -> Table S4
#   Early vs Control  -> Table S5
# 입력 : OUT_DIR/GSE142025.labeled.txt (01 스크립트 산출물)
# 출력 : RES_DIR/DEG_*.txt, diff_*.txt, vol_*.pdf
# ※ 결과는 results/ 로만 저장(원본 불변).

library(here)
source(here::here("config.R"))
suppressMessages({ library(limma); library(ggplot2) })

# STEP1 산출물 전용 폴더 (STEP 구분 저장)
S1_DIR <- file.path(RES_DIR, "step1_deg"); dir.create(S1_DIR, showWarnings = FALSE, recursive = TRUE)

rt <- read.table(file.path(OUT_DIR, "GSE142025.labeled.txt"),
                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
data <- as.matrix(rt)
data <- avereps(data)

# 자동 log2 (이미 로그면 스킵) + quantile 정규화
qx <- as.numeric(quantile(data, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
if ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)) {
  data[data < 0] <- 0; data <- log2(data + 1)
  message("[DEG] log2 변환 적용")
} else message("[DEG] 값이 이미 로그 스케일로 판단 -> log2 스킵")
data <- normalizeBetweenArrays(data)

grp <- sub(".*_", "", colnames(data))

run_contrast <- function(alt, ref, tag) {
  keep <- grp %in% c(alt, ref)
  d <- data[, keep]; g <- factor(grp[keep], levels = c(ref, alt))
  design <- model.matrix(~0 + g); colnames(design) <- levels(g)
  fit <- eBayes(contrasts.fit(lmFit(d, design),
                              makeContrasts(contrasts = paste0(alt, "-", ref), levels = design)))
  tab <- topTable(fit, number = Inf, adjust = "fdr")
  write.table(cbind(id = rownames(tab), tab),
              file.path(S1_DIR, paste0("DEG_all_", tag, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  sig <- tab[abs(tab$logFC) > LOGFC_FILTER & tab$adj.P.Val < ADJP_FILTER, ]
  write.table(cbind(id = rownames(sig), sig),
              file.path(S1_DIR, paste0("DEG_diff_", tag, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # 화산도
  tab$Sig <- ifelse(tab$adj.P.Val < ADJP_FILTER & abs(tab$logFC) > LOGFC_FILTER,
                    ifelse(tab$logFC > 0, "Up", "Down"), "Not")
  p <- ggplot(tab, aes(logFC, -log10(adj.P.Val), col = Sig)) +
    geom_point(size = .7) +
    scale_color_manual(values = c(Down="blue2", Not="grey", Up="red2")) +
    labs(title = gsub("_", " ", tag)) +
    theme(plot.title = element_text(hjust = .5, face = "bold"))
  ggsave(file.path(S1_DIR, paste0("vol_", tag, ".pdf")), p, width = 5.5, height = 4.5)
  cat(sprintf("[%s] DEG %d개 (up %d, down %d)\n",
              tag, nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))
  c(n = nrow(sig), up = sum(sig$logFC > 0), down = sum(sig$logFC < 0))
}

cat("그룹 구성: "); print(table(grp))
r_LvE <- run_contrast(GROUP_LATE,  GROUP_EARLY,   "Late_vs_Early")     # Table S3
r_LvC <- run_contrast(GROUP_LATE,  GROUP_CONTROL, "Late_vs_Control")   # Table S4
r_EvC <- run_contrast(GROUP_EARLY, GROUP_CONTROL, "Early_vs_Control")  # Table S5

# ---- 논문값 vs 우리값 대조 CSV (STEP1 GSE142025) ----
# 논문(Fig 2A / Supp S3): Late vs Early = up 1,557 / down 1,276 / total 2,833 (|logFC|>0.585, adjP<0.05)
cmp <- data.frame(
  metric      = c("GSE142025 Late_vs_Early DEG(total)", "  up", "  down",
                  "GSE142025 Late_vs_Control DEG(total)", "GSE142025 Early_vs_Control DEG(total)"),
  paper       = c(2833, 1557, 1276, NA, NA),
  ours        = c(r_LvE["n"], r_LvE["up"], r_LvE["down"], r_LvC["n"], r_EvC["n"]),
  stringsAsFactors = FALSE)
cmp$diff <- cmp$ours - cmp$paper
write.csv(cmp, file.path(S1_DIR, "compare_paper_vs_ours.GSE142025.csv"), row.names = FALSE)
cat("\n[대조] 논문 Late_vs_Early 2,833(up1557/down1276) vs 우리 ",
    r_LvE["n"], "(up", r_LvE["up"], "/down", r_LvE["down"], ")\n", sep="")
cat("[저장] ", file.path(S1_DIR, "compare_paper_vs_ours.GSE142025.csv"), "\n")
