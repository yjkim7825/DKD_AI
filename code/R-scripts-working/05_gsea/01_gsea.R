# 01_gsea.R -----------------------------------------------------------------
# STEP 5 : GSEA(Hallmark+KEGG) + ssGSEA.  원본 GSEA.R / hallmark.gsea.R 1:1 대조.
#   (1) Hallmark GSEA — 원본 hallmark.gsea.R 와 동일: GSE142025 DEG 3비교의 **diff_(유의 DEG)**
#       logFC 랭킹 → GSEA(TERM2GENE, pvalueCutoff=1). 제공 gmt(4-1) 사용(원본 msigdbr 대체, 오프라인).
#   (2) KEGG GSEA — 원본 GSEA.R 와 동일: 동일 diff_ DEG 를 KEGG gmt(4-2)로 GSEA. (기존 파이프라인 누락분 추가)
#   (3) KEGG ssGSEA — 논문 Fig 2H,I(단계별 경로점수) 방법: GSE96804 발현에 ssGSEA →
#       DKD vs Control Wilcoxon + FN1/ALDH2 발현 상관(멀티오믹스 연결). 원본 GSEA.R(=GSEA) 와는 다른 분석.
#   출력: results/step5_gsea/. 원본 데이터·스크립트 무수정.
# ---------------------------------------------------------------------------

suppressMessages({
  library(clusterProfiler); library(enrichplot); library(GSVA)
  library(ggplot2); library(limma)
})
library(here)
source(here::here("config.R"))
OUT <- file.path(RES_DIR, "step5_gsea"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

HALLMARK_GMT <- file.path(DATA_ROOT, "4-1. h.all.v2026.1.Hs.symbols.gmt")
KEGG_GMT     <- file.path(DATA_ROOT, "4-2. c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt")

## read.gmt 컬럼명 버전차 방어 (term/gene 로 표준화)
read_gmt_std <- function(path) {
  g <- read.gmt(path)
  names(g)[1:2] <- c("term", "gene")
  g
}

## ================= (1) Hallmark GSEA =================
hm <- read_gmt_std(HALLMARK_GMT)   # term2gene: (term, gene)
message("[Hallmark] gene sets: ", length(unique(hm$term)))

# 원본 GSEA.R/hallmark.gsea.R 와 동일: diff_(유의 DEG) 의 logFC 랭킹으로 GSEA.
read_diff_lfc <- function(degFile) {
  d <- read.table(file.path(RES_DIR, "step1_deg", degFile), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  lfc <- d$logFC; names(lfc) <- rownames(d)
  sort(lfc[!is.na(lfc) & !duplicated(names(lfc))], decreasing = TRUE)
}
run_gsea <- function(degFile, tag, term2gene, prefix) {
  lfc <- read_diff_lfc(degFile)
  kk <- GSEA(lfc, TERM2GENE = term2gene, pvalueCutoff = 1, seed = TRUE, verbose = FALSE)
  tab <- as.data.frame(kk)
  write.table(tab, file.path(OUT, paste0(prefix, ".GSEA.", tag, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  message("[", prefix, " ", tag, "] ranked ", length(lfc), " genes | total ", nrow(tab),
          " | sig(p.adj<0.05) ", sum(tab$p.adjust < 0.05))
  tab
}
# 원본과 동일하게 diff_ (유의 DEG) 입력
cmp <- list(Early_vs_Control = "DEG_diff_Early_vs_Control.txt",
            Late_vs_Control  = "DEG_diff_Late_vs_Control.txt",
            Late_vs_Early    = "DEG_diff_Late_vs_Early.txt")
hres <- Map(function(f, n) run_gsea(f, n, hm, "Hallmark"), cmp, names(cmp))

# FN1/ALDH2 관련 Hallmark 경로(EMT/ROS/염증 등) 요약
focus_pat <- "EPITHELIAL_MESENCHYMAL|REACTIVE_OXYGEN|OXIDATIVE|INFLAMMATORY|TGF_BETA|APOPTOSIS"
cat("\n== Hallmark FN1/ALDH2 관련 경로 (NES/p.adjust) ==\n")
focusTab <- do.call(rbind, lapply(names(hres), function(n) {
  t <- hres[[n]]; t <- t[grepl(focus_pat, t$ID), c("ID","NES","pvalue","p.adjust")]
  if (nrow(t)) cbind(Comparison = n, t) else NULL
}))
print(focusTab, row.names = FALSE)
write.csv(focusTab, file.path(OUT, "Hallmark_focus_pathways.csv"), row.names = FALSE)

# 논문 대조 CSV (Hallmark EMT/OXPHOS NES) — 논문 Fig2: EMT Late_vs_Ctrl +2.52 / Late_vs_Early +2.59, OXPHOS -1.62 / -2.02
getNES <- function(comp, pat) { t <- hres[[comp]]; v <- t$NES[grepl(pat, t$ID)]; if (length(v)) v[1] else NA }
hall_cmp <- data.frame(
  pathway = c("EMT", "EMT", "OXIDATIVE_PHOSPHORYLATION", "OXIDATIVE_PHOSPHORYLATION"),
  comparison = c("Late_vs_Control", "Late_vs_Early", "Late_vs_Control", "Late_vs_Early"),
  paper_NES = c(2.52, 2.59, -1.62, -2.02),
  ours_NES = c(getNES("Late_vs_Control", "EPITHELIAL_MESENCHYMAL"),
               getNES("Late_vs_Early",   "EPITHELIAL_MESENCHYMAL"),
               getNES("Late_vs_Control", "OXIDATIVE_PHOSPHORYLATION"),
               getNES("Late_vs_Early",   "OXIDATIVE_PHOSPHORYLATION")))
hall_cmp$ours_NES <- round(hall_cmp$ours_NES, 2)
write.csv(hall_cmp, file.path(OUT, "compare_paper_vs_ours.Hallmark_NES.csv"), row.names = FALSE)
cat("\n== Hallmark NES 논문 대조 ==\n"); print(hall_cmp, row.names = FALSE)

## ================= (2) KEGG GSEA (원본 GSEA.R) =================
kegg <- read_gmt_std(KEGG_GMT)
message("\n[KEGG GSEA] gene sets: ", length(unique(kegg$term)))
kres <- Map(function(f, n) run_gsea(f, n, kegg, "KEGG"), cmp, names(cmp))
# FN1(ECM)·ALDH2(TCA/대사) 관련 KEGG 경로 요약
kfocus_pat <- "ECM_RECEPTOR|CITRATE_CYCLE|TCA|OXIDATIVE_PHOSPHORYLATION|FOCAL_ADHESION|TRYPTOPHAN|ARGININE"
kfocusTab <- do.call(rbind, lapply(names(kres), function(n) {
  t <- kres[[n]]; t <- t[grepl(kfocus_pat, t$ID), c("ID","NES","pvalue","p.adjust")]
  if (nrow(t)) cbind(Comparison = n, t) else NULL
}))
if (!is.null(kfocusTab)) { cat("\n== KEGG GSEA FN1/ALDH2 관련 경로 ==\n"); print(kfocusTab, row.names = FALSE)
  write.csv(kfocusTab, file.path(OUT, "KEGG.GSEA_focus_pathways.csv"), row.names = FALSE) }

## ================= (3) KEGG ssGSEA (논문 Fig 2H,I) =================
keggList <- split(kegg$gene, kegg$term)
message("\n[KEGG ssGSEA] gene sets: ", length(keggList))

expr <- as.matrix(read.table(file.path(OUT_DIR, "GSE96804.labeled.txt"),
                             header = TRUE, sep = "\t", check.names = FALSE, row.names = 1))
grp <- ifelse(grepl("_Control$", colnames(expr)), "Control", "DKD")

par_ss <- ssgseaParam(expr, keggList)
ss <- gsva(par_ss)                       # 경로 x 샘플 점수
write.table(rbind(pathway = colnames(t(ss)), t(ss)), file.path(OUT, "KEGG.ssGSEA.scores.txt"),
            sep = "\t", quote = FALSE, col.names = FALSE)
message("[KEGG ssGSEA] 경로 x 샘플 = ", nrow(ss), " x ", ncol(ss))

# DKD vs Control Wilcoxon per pathway
wilcox_tab <- data.frame(pathway = rownames(ss),
  p = apply(ss, 1, function(v) wilcox.test(v ~ grp)$p.value),
  meanDKD = rowMeans(ss[, grp == "DKD", drop = FALSE]),
  meanCtrl = rowMeans(ss[, grp == "Control", drop = FALSE]))
wilcox_tab$padj <- p.adjust(wilcox_tab$p, "BH")
wilcox_tab$deltaDKD <- wilcox_tab$meanDKD - wilcox_tab$meanCtrl
wilcox_tab <- wilcox_tab[order(wilcox_tab$padj), ]
write.csv(wilcox_tab, file.path(OUT, "KEGG.ssGSEA.DKDvsControl.csv"), row.names = FALSE)
cat("\n== KEGG ssGSEA DKD vs Control (top 12 by padj) ==\n")
print(head(wilcox_tab[, c("pathway","deltaDKD","padj")], 12), row.names = FALSE)

# FN1/ALDH2 발현과 경로점수 상관
cor_focus <- function(gene) {
  if (!(gene %in% rownames(expr))) return(NULL)
  g <- as.numeric(expr[gene, ])
  cc <- apply(ss, 1, function(v) cor(v, g, method = "spearman"))
  data.frame(gene = gene, pathway = names(cc), rho = as.numeric(cc))
}
corTab <- rbind(cor_focus("FN1"), cor_focus("ALDH2"))
corTab <- corTab[order(-abs(corTab$rho)), ]
write.csv(corTab, file.path(OUT, "KEGG.ssGSEA.FN1_ALDH2_corr.csv"), row.names = FALSE)
cat("\n== FN1/ALDH2 vs KEGG ssGSEA 상관 상위 (|rho|) ==\n")
print(head(corTab, 12), row.names = FALSE)

message("\n[STEP5] 완료 -> ", OUT)
