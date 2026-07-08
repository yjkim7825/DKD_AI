# 11_step5_gsea.R -----------------------------------------------------------
# STEP 5 : Hallmark GSEA + KEGG ssGSEA (원본 hallmark.gsea.R 로직)
#   (1) Hallmark GSEA — GSE142025 DEG 3비교(Early/Late vs Control, Late vs Early)의 logFC 랭킹.
#       제공 gmt(../data/4-1. h.all...) 사용(원본은 msigdbr; 오프라인이라 gmt 파일로 대체).
#       FN1(EMT·산화)·ALDH2(산화) 관련 Hallmark 경로 유의성 확인.
#   (2) KEGG ssGSEA — GSE96804(train) 발현에 KEGG(../data/4-2) ssGSEA 로 샘플별 경로점수.
#       DKD vs Control Wilcoxon + FN1/ALDH2 발현과 상관. FN1/ALDH2 관련 KEGG 경로 확인.
#   출력: results/step5_gsea/
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

run_gsea <- function(degFile, tag) {
  d <- read.table(file.path(RES_DIR, degFile), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  lfc <- d$logFC; names(lfc) <- rownames(d)
  lfc <- sort(lfc[!is.na(lfc) & !duplicated(names(lfc))], decreasing = TRUE)
  kk <- GSEA(lfc, TERM2GENE = hm, pvalueCutoff = 1, seed = TRUE, verbose = FALSE)
  tab <- as.data.frame(kk)
  write.table(tab, file.path(OUT, paste0("Hallmark.GSEA.", tag, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  message("[Hallmark ", tag, "] total ", nrow(tab), " | sig(p.adj<0.05) ", sum(tab$p.adjust < 0.05))
  tab
}
cmp <- list(Early_vs_Control = "DEG_all_Early_vs_Control.txt",
            Late_vs_Control  = "DEG_all_Late_vs_Control.txt",
            Late_vs_Early    = "DEG_all_Late_vs_Early.txt")
hres <- Map(function(f, n) run_gsea(f, n), cmp, names(cmp))

# FN1/ALDH2 관련 Hallmark 경로(EMT/ROS/염증 등) 요약
focus_pat <- "EPITHELIAL_MESENCHYMAL|REACTIVE_OXYGEN|OXIDATIVE|INFLAMMATORY|TGF_BETA|APOPTOSIS"
cat("\n== Hallmark FN1/ALDH2 관련 경로 (NES/p.adjust) ==\n")
focusTab <- do.call(rbind, lapply(names(hres), function(n) {
  t <- hres[[n]]; t <- t[grepl(focus_pat, t$ID), c("ID","NES","pvalue","p.adjust")]
  if (nrow(t)) cbind(Comparison = n, t) else NULL
}))
print(focusTab, row.names = FALSE)
write.csv(focusTab, file.path(OUT, "Hallmark_focus_pathways.csv"), row.names = FALSE)

## ================= (2) KEGG ssGSEA =================
kegg <- read_gmt_std(KEGG_GMT)
keggList <- split(kegg$gene, kegg$term)
message("\n[KEGG] gene sets: ", length(keggList))

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
