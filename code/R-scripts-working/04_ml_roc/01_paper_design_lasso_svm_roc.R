# 10_step34_paper_design.R --------------------------------------------------
# 논문 본문(Figure 4) 확정 설계로 STEP 3-4 재실행.
#   훈련 = GSE96804 단독,  검증 = ComBat(GSE104948 + GSE104954)  (GSE30529 제외)
#   후보 = 저자 확정 10개 (Supp Table 8/9/10/12 "candidate genes"):
#          ALDH2,FN1,VNN2,CREB5,XAF1,CA2,CDKN1B,IFI44L,SYTL2,TSPYL5
#          (= MR risk∩DKD고발현 ∪ MR protective∩DKD저발현, MR 견고성 검증 통과)
#   4) LASSO(10) → 6개(CDKN1B,ALDH2,FN1,XAF1,TSPYL5,VNN2) 재현 확인
#   5) 단일유전자 ROC(train/valid) 논문 대조
#   6) FN1+ALDH2 결합모델 GLM/RF/SVM/XGBoost → Figure 4H,I 대조
#   출력: results/step4_paper/
# ---------------------------------------------------------------------------

suppressMessages({
  library(limma); library(sva); library(glmnet); library(pROC)
  library(randomForest); library(e1071); library(xgboost)
})
library(here)
source(here::here("config.R"))
OUT <- file.path(RES_DIR, "step4_paper"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

CAND10 <- c("ALDH2","FN1","VNN2","CREB5","XAF1","CA2","CDKN1B","IFI44L","SYTL2","TSPYL5")
EXPECT6 <- c("CDKN1B","ALDH2","FN1","XAF1","TSPYL5","VNN2")
writeLines(CAND10, file.path(OUT, "interGenes.List.txt"))

read_lab <- function(p) { m <- as.matrix(read.table(p, header=TRUE, sep="\t", check.names=FALSE, row.names=1)); storage.mode(m)<-"double"; m }

## ---- 훈련: GSE96804 단독 ----
TRAIN <- read_lab(file.path(OUT_DIR, "GSE96804.labeled.txt"))
message("[train] GSE96804: ", nrow(TRAIN), " x ", ncol(TRAIN))

## ---- 검증: ComBat(GSE104948 + GSE104954) ----
v1 <- read_lab(file.path(OUT_DIR, "GSE104948.labeled.txt")); colnames(v1) <- paste0("GSE104948_", colnames(v1))
v2 <- read_lab(file.path(OUT_DIR, "GSE104954.labeled.txt")); colnames(v2) <- paste0("GSE104954_", colnames(v2))
common <- intersect(rownames(v1), rownames(v2))
vall <- cbind(v1[common,], v2[common,]); vbatch <- c(rep(1,ncol(v1)), rep(2,ncol(v2)))
keep <- apply(vall, 1, function(z) all(tapply(z, vbatch, function(w) sd(w)>0)))
vall <- vall[keep,]
VALID <- ComBat(dat = vall, batch = vbatch, par.prior = TRUE)
write.table(rbind(geneNames=colnames(VALID), VALID), file.path(OUT_DIR,"data.valid.paper.txt"),
            sep="\t", quote=FALSE, col.names=FALSE)
message("[valid] ComBat(104948+104954): ", nrow(VALID), " x ", ncol(VALID),
        " (C/D: ", paste(table(ifelse(grepl("_Control$",colnames(VALID)),"Control","DKD")),collapse="/"), ")")

ytr <- ifelse(grepl("_Control$", colnames(TRAIN)), 0, 1)
yva <- ifelse(grepl("_Control$", colnames(VALID)), 0, 1)

## 후보 10개 존재 확인
cat("\n== 후보 10 존재 (train/valid) ==\n")
for(g in CAND10) cat(sprintf("  %-8s train=%s valid=%s\n", g, g%in%rownames(TRAIN), g%in%rownames(VALID)))
feats <- intersect(CAND10, rownames(TRAIN))

## ---- 4) LASSO(10) → ? ----
x <- t(TRAIN[feats,]); yf <- factor(ifelse(ytr==0,"Control","DKD"))
set.seed(123)
cv <- cv.glmnet(as.matrix(x), yf, family="binomial", alpha=1, type.measure="deviance", nfolds=10)
fit <- glmnet(as.matrix(x), yf, family="binomial", alpha=1)
co <- coef(fit, s=cv$lambda.min); lassoGene <- setdiff(rownames(co)[which(co!=0)], "(Intercept)")
writeLines(lassoGene, file.path(OUT,"LASSO.gene.txt"))
cat("\n== LASSO ==\n")
cat("  선택(",length(lassoGene),"):", paste(sort(lassoGene),collapse=", "), "\n")
cat("  기대 6개와 일치:", setequal(lassoGene, EXPECT6),
    "| 기대-선택 차:", paste(setdiff(EXPECT6,lassoGene),collapse=","),
    "| 선택-기대 차:", paste(setdiff(lassoGene,EXPECT6),collapse=","), "\n")

## ---- 5) 단일유전자 ROC (train/valid) ----
aucf <- function(mat, y, g) if(!(g%in%rownames(mat))) NA_real_ else as.numeric(auc(roc(y, as.numeric(mat[g,]), quiet=TRUE)))
roc_tab <- data.frame(gene=CAND10,
  AUC_train=sapply(CAND10, function(g) aucf(TRAIN,ytr,g)),
  AUC_valid=sapply(CAND10, function(g) aucf(VALID,yva,g)), row.names=NULL)
write.csv(roc_tab, file.path(OUT,"single_gene_ROC_AUC.csv"), row.names=FALSE)
cat("\n== 단일유전자 ROC-AUC ==\n"); print(roc_tab, row.names=FALSE)

## ---- 6) FN1+ALDH2 결합 모델 (GLM/RF/SVM/XGBoost) ----
dtr <- data.frame(FN1=TRAIN["FN1",], ALDH2=TRAIN["ALDH2",], y=ytr)
dva <- data.frame(FN1=VALID["FN1",], ALDH2=VALID["ALDH2",], y=yva)
comb <- function(train_p, valid_p) c(train=as.numeric(auc(roc(dtr$y, train_p, quiet=TRUE))),
                                     valid=as.numeric(auc(roc(dva$y, valid_p, quiet=TRUE))))
res <- list()
# GLM (logistic)
g <- glm(y~FN1+ALDH2, data=dtr, family=binomial)
res$GLM <- comb(predict(g, dtr, type="response"), predict(g, dva, type="response"))
# RF
set.seed(123); rf <- randomForest(x=dtr[,c("FN1","ALDH2")], y=factor(dtr$y), ntree=500)
res$RF <- comb(predict(rf, dtr[,c("FN1","ALDH2")], type="prob")[,2],
               predict(rf, dva[,c("FN1","ALDH2")], type="prob")[,2])
# SVM
set.seed(123); sv <- svm(x=dtr[,c("FN1","ALDH2")], y=factor(dtr$y), probability=TRUE, kernel="radial")
pp_tr <- attr(predict(sv, dtr[,c("FN1","ALDH2")], probability=TRUE),"probabilities")[,"1"]
pp_va <- attr(predict(sv, dva[,c("FN1","ALDH2")], probability=TRUE),"probabilities")[,"1"]
res$SVM <- comb(pp_tr, pp_va)
# XGBoost (안정 저수준 API: xgb.DMatrix + xgb.train)
set.seed(123)
dtrain <- xgb.DMatrix(data=as.matrix(dtr[,c("FN1","ALDH2")]), label=dtr$y)
xgb <- xgb.train(params=list(objective="binary:logistic", max_depth=3, eta=0.3),
                 data=dtrain, nrounds=50, verbose=0)
res$XGBoost <- comb(predict(xgb, as.matrix(dtr[,c("FN1","ALDH2")])),
                    predict(xgb, as.matrix(dva[,c("FN1","ALDH2")])))
model_tab <- data.frame(model=names(res),
  AUC_train=round(sapply(res,function(v)v["train"]),3),
  AUC_valid=round(sapply(res,function(v)v["valid"]),3), row.names=NULL)
write.csv(model_tab, file.path(OUT,"combined_model_AUC.csv"), row.names=FALSE)
cat("\n== FN1+ALDH2 결합모델 AUC ==\n"); print(model_tab, row.names=FALSE)
message("\n[완료] 산출 -> ", OUT)
