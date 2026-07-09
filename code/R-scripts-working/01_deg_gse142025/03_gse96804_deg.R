# 03_gse96804_deg.R --------------------------------------------------------
# 논문이 명시한 "complementary" DEG: GSE96804(DKD vs Control) limma DEG.
#   논문 원문(Methods, p.4): "Complementary differential expression analysis of the
#   GSE96804 dataset using limma (v3.64.1) identified DEGs at adjusted p < 0.05 with
#   |logFC| > 0.585." → MR 통합(고발현 logFC>0.585 ∩ MR risk / 저발현 logFC<-0.585 ∩ MR protective)의 근거.
# 이 DEG 로부터 FN1(DKD 고발현·risk)·ALDH2(DKD 저발현·protective) 방향이 정의됨.
# 입력 : OUT_DIR/GSE96804.labeled.txt (STEP2 RMA 산출물, 열 = '{GSM}_{Control|DKD}')
# 출력 : RES_DIR/step1_deg/  (DEG 전체표/유의표 + FN1/ALDH2 방향 확인 + 대조 CSV)
# ※ 원본 불변, 결과는 results/ 로만.

library(here)
source(here::here("config.R"))
suppressMessages({ library(limma); library(ggplot2) })

S1_DIR <- file.path(RES_DIR, "step1_deg"); dir.create(S1_DIR, showWarnings = FALSE, recursive = TRUE)

rt <- read.table(file.path(OUT_DIR, "GSE96804.labeled.txt"),
                 header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
data <- avereps(as.matrix(rt))

# 이미 RMA(log2) 산출물 → log2 스킵. 어레이 간 정규화만.
data <- normalizeBetweenArrays(data)
grp  <- factor(sub(".*_", "", colnames(data)), levels = c(GROUP_CONTROL, GROUP_DKD))
cat("[GSE96804] 그룹 구성: "); print(table(grp))

design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
fit <- eBayes(contrasts.fit(lmFit(data, design),
              makeContrasts(contrasts = paste0(GROUP_DKD, "-", GROUP_CONTROL), levels = design)))
tab <- topTable(fit, number = Inf, adjust = "fdr")

write.table(cbind(id = rownames(tab), tab),
            file.path(S1_DIR, "DEG_all_GSE96804_DKD_vs_Control.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
sig <- tab[abs(tab$logFC) > LOGFC_FILTER & tab$adj.P.Val < ADJP_FILTER, ]
write.table(cbind(id = rownames(sig), sig),
            file.path(S1_DIR, "DEG_diff_GSE96804_DKD_vs_Control.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# MR 통합 근거: DKD 고발현(up) / 저발현(down)
up   <- rownames(sig[sig$logFC > 0, ]);  down <- rownames(sig[sig$logFC < 0, ])
writeLines(up,   file.path(S1_DIR, "GSE96804.DKD_high_genes.txt"))
writeLines(down, file.path(S1_DIR, "GSE96804.DKD_low_genes.txt"))

# 화산도
tab$Sig <- ifelse(tab$adj.P.Val < ADJP_FILTER & abs(tab$logFC) > LOGFC_FILTER,
                  ifelse(tab$logFC > 0, "Up", "Down"), "Not")
p <- ggplot(tab, aes(logFC, -log10(adj.P.Val), col = Sig)) + geom_point(size = .6) +
  scale_color_manual(values = c(Down="blue2", Not="grey", Up="red2")) +
  labs(title = "GSE96804 DKD vs Control (limma)") +
  theme(plot.title = element_text(hjust = .5, face = "bold"))
ggsave(file.path(S1_DIR, "vol_GSE96804_DKD_vs_Control.pdf"), p, width = 5.5, height = 4.5)

# FN1/ALDH2 방향 확인 (논문: FN1 = DKD 고발현·risk, ALDH2 = DKD 저발현·protective)
foc <- intersect(c("FN1", "ALDH2"), rownames(tab))
focus <- tab[foc, c("logFC", "adj.P.Val", "Sig"), drop = FALSE]
write.csv(cbind(gene = rownames(focus), focus),
          file.path(S1_DIR, "GSE96804.FN1_ALDH2_direction.csv"), row.names = FALSE)

cat(sprintf("[GSE96804 DKD vs Control] DEG %d개 (up %d / down %d)\n",
            nrow(sig), length(up), length(down)))
cat("[FN1/ALDH2 방향]\n"); print(focus)
cat("[저장] ", S1_DIR, " (DEG_all/diff_GSE96804*, DKD_high/low_genes, vol, FN1_ALDH2_direction)\n", sep="")
