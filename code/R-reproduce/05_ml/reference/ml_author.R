# ============================================================================
# ml_author.R  — 방식 B (저자 업로드 코드 방식)
#   저자 "machine learning modeling 1.R" + "2.R" 재현:
#     LASSO ∩ SVM-RFE 교집합으로 유전자 선택 → 유전자별 학습/검증 ROC
#   입력 : data/processed/data.train.paper.txt, data.valid.paper.txt,
#          05_ml/interGenes.List.txt
#   출력 : output/B_LASSO.gene.txt, B_SVM-RFE.gene.txt, B_interGenes.txt,
#          B_venn.pdf, B_SVM-RFE.pdf, B_intersect_roc.csv, B_ROC_<gene>.pdf
#   ※ 참고용 재현: FN1·ALDH2는 잡지만 SVM-RFE가 덜 좁혀 최종 목록이 넓게 남음
#      → 결론적으로 논문 방식(01)이 더 깔끔해 그쪽을 채택.
# ============================================================================
set.seed(123)
suppressMessages({
  library(glmnet); library(e1071); library(kernlab); library(caret)
  library(VennDiagram); library(grid); library(pROC)
})

ROOT   <- "C:/Users/dbdl0/Downloads/ETRI_LAB/DKD_AI/code"
MLDIR  <- file.path(ROOT, "R-reproduce/05_ml")
DATA   <- file.path(ROOT, "data/processed")
OUT    <- file.path(MLDIR, "output"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(MLDIR, "R/ml_func.R"))

interGenes <- scan(file.path(MLDIR, "interGenes.List.txt"), what = "character", quiet = TRUE)
tr <- load_expr(file.path(DATA, "data.train.paper.txt"), interGenes)
te <- load_expr(file.path(DATA, "data.valid.paper.txt"), interGenes)
sc <- scale_train_test(tr$X, te$X); Xtr <- sc$tr; Xte <- sc$te
ytr <- tr$y; yte <- te$y; genes <- tr$genes

## ── B1) LASSO (저자 1.R) ────────────────────────────────────────────────────
cvfit <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1,
                   type.measure = "deviance", nfolds = 10)
co <- coef(cvfit, s = "lambda.min"); sel <- rownames(co)[which(as.numeric(co) != 0)]
lassoGene <- setdiff(sel, "(Intercept)")
write.table(lassoGene, file.path(OUT, "B_LASSO.gene.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("[B1 LASSO]", length(lassoGene), "개:", paste(lassoGene, collapse = ", "), "\n")

## ── B2) SVM-RFE (저자 1.R: caret rfe + svmRadial) ───────────────────────────
group <- factor(ytr, levels = c("Control", "DKD"))
Profile <- rfe(x = Xtr, y = as.numeric(group),
               sizes = c(2, 3, 4, 5, 6, 7, 8),
               rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
               methods = "svmRadial", preProcess = c("center", "scale"))
pdf(file.path(OUT, "B_SVM-RFE.pdf"), width = 6, height = 5.5); par(las = 1)
plot(Profile$results$Variables, Profile$results$RMSE, type = "b",
     xlab = "Variables", ylab = "RMSE (CV)", col = "darkgreen")
wmin <- which.min(Profile$results$RMSE)
points(Profile$results$Variables[wmin], Profile$results$RMSE[wmin], col = "blue", pch = 16)
dev.off()
svmGene <- Profile$optVariables
write.table(svmGene, file.path(OUT, "B_SVM-RFE.gene.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("[B2 SVM-RFE]", length(svmGene), "개:", paste(svmGene, collapse = ", "), "\n")

## ── B3) 교집합 (저자 1.R: Venn) ─────────────────────────────────────────────
geneList <- list(LASSO = unique(lassoGene), `SVM-RFE` = unique(svmGene))
vp <- venn.diagram(x = geneList, filename = NULL,
                   fill = c("cornflowerblue", "darkorchid1"), alpha = 0.5,
                   cat.col = c("cornflowerblue", "darkorchid1"), margin = 0.1)
pdf(file.path(OUT, "B_venn.pdf"), width = 5, height = 5); grid.draw(vp); dev.off()
interG <- Reduce(intersect, geneList)
write.table(interG, file.path(OUT, "B_interGenes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("[B3 교집합]", paste(interG, collapse = ", "), "\n")

## ── B4) 교집합 유전자별 ROC (저자 2.R: 학습·검증 둘 다) ──────────────────────
#   저자 ML2.R처럼 유전자마다 training_ROC / test_ROC 각각 그림 (AUC + 95%CI)
roc_one <- function(y, expr, g, tag) {
  r <- roc(y, as.numeric(expr), quiet = TRUE); ci <- as.numeric(ci.auc(r, method = "bootstrap"))
  pdf(file.path(OUT, paste0("B_ROC_", g, "_", tag, ".pdf")), width = 3.5, height = 3.5)
  plot(r, print.auc = TRUE, col = "red", legacy.axes = TRUE,
       main = paste0(g, if (tag == "train") " of Training" else " of Validation"))
  text(0.34, 0.38, sprintf("95%% CI: %.3f-%.3f", ci[1], ci[3]), col = "red")
  dev.off()
}
roc_tab <- data.frame()
for (g in interG) {
  a_tr <- gene_auc(ytr, Xtr[, g]); a_te <- gene_auc(yte, Xte[, g])
  roc_tab <- rbind(roc_tab, data.frame(gene = g,
                     AUC_train = round(a_tr, 3), AUC_test = round(a_te, 3)))
  roc_one(ytr, Xtr[, g], g, "train")      # 저자: training_ROC.<gene>
  roc_one(yte, Xte[, g], g, "valid")      # 저자: test_ROC.<gene>
}
write.csv(roc_tab, file.path(OUT, "B_intersect_roc.csv"), row.names = FALSE)
cat("\n[B4 교집합 유전자별 ROC]\n"); print(roc_tab, row.names = FALSE)
cat("\n★ [방식 B] 완료 — 저자 코드(LASSO∩SVM-RFE) 재현. 참고용.\n")
