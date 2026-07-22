# ============================================================================
# 01_ml_feature.R  — 방식 A (논문 Figure 4), 1/2단계: 유전자 선별
#   배치보정(sva ComBat) → LASSO(10겹 CV) → 학습·검증 ROC-AUC>0.8 필터
#   → 최종 진단 유전자(FN1·ALDH2) 확정
#   입력 : data/processed/data.{train,valid}.paper.txt, interGenes.List.txt(후보 10개)
#   출력 : output/A_PCA_batch.pdf(Fig4A) · A_cvfit.pdf(Fig4C) · A_lasso.gene.txt
#          · A_roc_filter.csv · A_ROC_<gene>.pdf · A_final_genes.txt(→ 02가 읽음)
# ============================================================================
set.seed(123)
suppressMessages({ library(glmnet); library(pROC); library(sva) })

ROOT  <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MLDIR <- file.path(ROOT, "R-reproduce/05_ml")
DATA  <- file.path(ROOT, "data/processed")
OUT   <- file.path(MLDIR, "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(MLDIR, "R/ml_func.R"))

interGenes <- scan(file.path(MLDIR, "interGenes.List.txt"), what = "character", quiet = TRUE)

## ── A0) 배치보정 (논문 Methods: sva::ComBat, batch=학습/검증) ────────────────
d   <- load_and_combat(file.path(DATA, "data.train.paper.txt"),
                       file.path(DATA, "data.valid.paper.txt"), interGenes)
sc  <- scale_train_test(d$Xtr, d$Xte); Xtr <- sc$tr; Xte <- sc$te
ytr <- d$ytr; yte <- d$yte; genes <- d$genes
cat(sprintf("[데이터] 학습 %d명 / 검증 %d명 / 유전자 %d개 (ComBat 보정됨)\n",
            nrow(Xtr), nrow(Xte), length(genes)))

# Fig 4A: 배치보정 후 PCA (색=데이터셋, 모양=Train/Valid)
pca <- prcomp(t(d$full), scale. = TRUE)
grp <- ifelse(d$isTr, "Train", "Valid"); bat <- d$batch
pdf(file.path(OUT, "A_PCA_batch.pdf"), width = 6, height = 5.5)
plot(pca$x[, 1], pca$x[, 2], col = as.integer(factor(bat)),
     pch = ifelse(grp == "Train", 16, 17), xlab = "PC1", ylab = "PC2",
     main = "PCA after ComBat batch correction")
legend("topright", legend = levels(factor(bat)), col = seq_along(levels(factor(bat))),
       pch = 16, cex = 0.8)
dev.off()

## ── A1) LASSO (10겹 CV, λmin) — 계수 0 아닌 유전자 선택 ──────────────────────
cvfit <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1,
                   type.measure = "deviance", nfolds = 10)
pdf(file.path(OUT, "A_cvfit.pdf"), width = 6, height = 5.5); plot(cvfit); dev.off()
co <- coef(cvfit, s = "lambda.min"); sel <- rownames(co)[which(as.numeric(co) != 0)]
lassoGene <- setdiff(sel, "(Intercept)")
write.table(lassoGene, file.path(OUT, "A_lasso.gene.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("[A1 LASSO] 선택", length(lassoGene), "개:", paste(lassoGene, collapse = ", "), "\n")

## ── A2) ROC 필터: 학습·검증 둘 다 AUC>0.8 → 최종 유전자 ─────────────────────
roc_tab <- data.frame()
for (g in lassoGene) {
  a_tr <- gene_auc(ytr, Xtr[, g]); a_te <- gene_auc(yte, Xte[, g])
  roc_tab <- rbind(roc_tab, data.frame(gene = g, AUC_train = round(a_tr, 3),
                                       AUC_test = round(a_te, 3),
                                       pass = (a_tr > 0.8 & a_te > 0.8)))
  pdf(file.path(OUT, paste0("A_ROC_", g, ".pdf")), width = 3.5, height = 3.5)
  plot(roc(ytr, as.numeric(Xtr[, g]), quiet = TRUE), col = "#4575b4", legacy.axes = TRUE,
       main = g, print.auc = TRUE, print.auc.y = 0.45)
  plot(roc(yte, as.numeric(Xte[, g]), quiet = TRUE), col = "#f16913", add = TRUE,
       print.auc = TRUE, print.auc.y = 0.35, print.auc.col = "#f16913")
  legend("bottomright", c("Train", "Valid"), col = c("#4575b4", "#f16913"), lwd = 2, cex = 0.7)
  dev.off()
}
write.csv(roc_tab, file.path(OUT, "A_roc_filter.csv"), row.names = FALSE)
finalGenes <- roc_tab$gene[roc_tab$pass]
write.table(finalGenes, file.path(OUT, "A_final_genes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)   # → 02가 읽음
cat("[A2 ROC>0.8] 최종 유전자:", paste(finalGenes, collapse = ", "), "\n")
print(roc_tab, row.names = FALSE)
cat("\n★ [01 선별] 완료 → A_final_genes.txt (다음: 02_ml_model.R)\n")
