# 12_step4_roc_plots.R ------------------------------------------------------
# STEP 4(논문설계) ROC 그림 저장 — 논문 Figure 4 스타일 (pROC). 데이터 재수신 없음.
#   훈련 = GSE96804.labeled.txt, 검증 = data.valid.paper.txt (기존 processed)
#   (1) 단일유전자 ROC 4개: {train,valid} x {FN1,ALDH2}, AUC+95%CI(ci.auc bootstrap) 라벨, 빨강+대각선
#   (2) 결합모델 ROC 2개: combined_ROC_{train,valid}.pdf, RF/SVM/GLM/XGBoost 4곡선 범례 AUC
#   (3) 그림 AUC/CI 를 single_gene_ROC_AUC.csv / combined_model_AUC.csv 와 대조
#   ※ 모델 학습은 10_step34_paper_design.R 과 동일(set.seed(123)) → CSV 와 일치해야 함
# ---------------------------------------------------------------------------

suppressMessages({ library(pROC); library(randomForest); library(e1071); library(xgboost) })
library(here)
source(here::here("config.R"))
OUT <- file.path(RES_DIR, "step4_paper")
set.seed(123)

read_lab <- function(p){ m<-as.matrix(read.table(p,header=TRUE,sep="\t",check.names=FALSE,row.names=1)); storage.mode(m)<-"double"; m }
TRAIN <- read_lab(file.path(OUT_DIR,"GSE96804.labeled.txt"))
VALID <- read_lab(file.path(OUT_DIR,"data.valid.paper.txt"))
ytr <- ifelse(grepl("_Control$", colnames(TRAIN)),0,1)
yva <- ifelse(grepl("_Control$", colnames(VALID)),0,1)

## ---- (1) 단일유전자 ROC ----
draw_single <- function(mat, y, gene, tag) {
  set.seed(123)
  r <- roc(y, as.numeric(mat[gene,]), quiet=TRUE)
  ci <- ci.auc(r, method="bootstrap")           # (lower, auc, upper)
  pdf(file.path(OUT, paste0(tag,"_ROC.",gene,".pdf")), width=3.6, height=3.6)
  par(pty="s", mar=c(4,4,3,1))
  plot(r, col="red", lwd=2.5, legacy.axes=TRUE,
       main=paste0(gene, " (", tag, ")"),
       xlab="1 - Specificity", ylab="Sensitivity")
  abline(0,1, lty=3, col="grey60")
  text(0.55, 0.18, sprintf("AUC: %.3f", as.numeric(ci[2])), col="red", cex=0.95, font=2)
  text(0.55, 0.08, sprintf("95%% CI: %.3f-%.3f", as.numeric(ci[1]), as.numeric(ci[3])), col="red", cex=0.8)
  dev.off()
  data.frame(fig=paste0(tag,"_ROC.",gene,".pdf"), gene=gene, set=tag,
             AUC=round(as.numeric(ci[2]),3), CIlow=round(as.numeric(ci[1]),3), CIhigh=round(as.numeric(ci[3]),3))
}
single <- rbind(
  draw_single(TRAIN, ytr, "FN1",   "train"), draw_single(VALID, yva, "FN1",   "valid"),
  draw_single(TRAIN, ytr, "ALDH2", "train"), draw_single(VALID, yva, "ALDH2", "valid"))

## ---- (2) 결합모델 (FN1+ALDH2) ROC — 학습은 script10 과 동일 ----
dtr <- data.frame(FN1=TRAIN["FN1",], ALDH2=TRAIN["ALDH2",], y=ytr)
dva <- data.frame(FN1=VALID["FN1",], ALDH2=VALID["ALDH2",], y=yva)
set.seed(123); g <- glm(y~FN1+ALDH2, data=dtr, family=binomial)
set.seed(123); rf <- randomForest(x=dtr[,c("FN1","ALDH2")], y=factor(dtr$y), ntree=500)
set.seed(123); sv <- svm(x=dtr[,c("FN1","ALDH2")], y=factor(dtr$y), probability=TRUE, kernel="radial")
set.seed(123); xgb <- xgb.train(params=list(objective="binary:logistic",max_depth=3,eta=0.3),
                                 data=xgb.DMatrix(as.matrix(dtr[,c("FN1","ALDH2")]),label=dtr$y), nrounds=50, verbose=0)
prob <- function(which){
  X<-if(which=="tr") dtr else dva
  list(
    GLM = predict(g, X, type="response"),
    RF  = predict(rf, X[,c("FN1","ALDH2")], type="prob")[,2],
    SVM = attr(predict(sv, X[,c("FN1","ALDH2")], probability=TRUE),"probabilities")[,"1"],
    XGBoost = predict(xgb, as.matrix(X[,c("FN1","ALDH2")])))
}
cols <- c(GLM="#1b9e77", RF="#d95f02", SVM="#7570b3", XGBoost="#e7298a")
draw_combined <- function(which, y, tag) {
  pr <- prob(which); aucs <- c()
  pdf(file.path(OUT, paste0("combined_ROC_",tag,".pdf")), width=4.2, height=4.2)
  par(pty="s", mar=c(4,4,3,1)); first<-TRUE
  for(m in names(pr)){
    r <- roc(y, pr[[m]], quiet=TRUE); aucs[m]<-as.numeric(auc(r))
    plot(r, col=cols[m], lwd=2.2, legacy.axes=TRUE, add=!first,
         main=paste0("FN1+ALDH2 combined (", tag, ")"),
         xlab="1 - Specificity", ylab="Sensitivity"); first<-FALSE
  }
  abline(0,1, lty=3, col="grey60")
  legend("bottomright", bty="n", lwd=2.2, col=cols[names(pr)],
         legend=sprintf("%-8s AUC=%.3f", names(pr), aucs[names(pr)]), cex=0.8)
  dev.off()
  data.frame(fig=paste0("combined_ROC_",tag,".pdf"), model=names(pr), set=tag, AUC=round(aucs[names(pr)],3))
}
comb <- rbind(draw_combined("tr", ytr, "train"), draw_combined("va", yva, "valid"))

## ---- (3) CSV 대조 ----
cat("\n== 단일유전자 그림 AUC vs single_gene_ROC_AUC.csv ==\n")
csv1 <- read.csv(file.path(OUT,"single_gene_ROC_AUC.csv"))
for(i in 1:nrow(single)){
  g<-single$gene[i]; s<-single$set[i]; col<-if(s=="train")"AUC_train" else "AUC_valid"
  csvv<-round(csv1[csv1$gene==g, col],3)
  cat(sprintf("  %-6s %-5s 그림=%.3f csv=%.3f 일치=%s (CI %.3f-%.3f)\n",
      g,s,single$AUC[i],csvv, isTRUE(all.equal(single$AUC[i],csvv,tolerance=0.02)), single$CIlow[i], single$CIhigh[i]))
}
cat("\n== 결합모델 그림 AUC vs combined_model_AUC.csv ==\n")
csv2 <- read.csv(file.path(OUT,"combined_model_AUC.csv"))
for(i in 1:nrow(comb)){
  m<-comb$model[i]; s<-comb$set[i]; col<-if(s=="train")"AUC_train" else "AUC_valid"
  csvv<-round(csv2[csv2$model==m,col],3)
  cat(sprintf("  %-8s %-5s 그림=%.3f csv=%.3f 일치=%s\n", m,s,comb$AUC[i],csvv, isTRUE(all.equal(comb$AUC[i],csvv,tolerance=0.02))))
}
write.csv(single, file.path(OUT,"ROC_figures_index.csv"), row.names=FALSE)  # 그림 인덱스(01의 10유전자 CSV 는 보존)
cat("\n[완료] ROC 그림 저장 ->", OUT, "\n")
