# 01_paper_design_lasso_svm_roc.R -------------------------------------------
# STEP 4 : LASSO + SVM-RFE 교집합 → 단일 ROC → 결합 진단모델 (논문 Figure 4).
#
# ▶ 원본 대조 (machine learning modeling 1·2.R):
#   · LASSO   : glmnet(family="binomial", alpha=1) + cv.glmnet(type.measure="deviance", nfolds=10),
#               coef at lambda.min, 절편 제거 → LASSO.gene.txt.  (원본과 동일)
#   · SVM-RFE : rfe(caretFuncs, method="cv", sizes=2:8, preProcess=center/scale) → SVM-RFE.gene.txt.
#               (원본 오타 methods= → 올바른 method="svmRadial" 로 정렬)
#   · Venn    : intersect(LASSO, SVM-RFE) → interGenes.txt.  (원본과 동일)
#   · ROC     : 유전자별 roc()+ci.auc(bootstrap), train/valid.  (원본 machine learning modeling 2.R 와 동일)
#
# ▶ 논문 Figure 4 설계 (원본 ML 경로와 다른 부분 → 논문 채택):
#   · 훈련 = GSE96804 단독,  검증 = STEP3 data.valid.paper.txt (=ComBat(104948+104954), 71).
#   · 후보 = 저자 확정 10개(Supp8/9/10/12).  · 논문 본문: LASSO→ROC(AUC>0.8 both)→결합모델(FN1+ALDH2).
#   · 결합모델 GLM/RF/SVM/XGBoost = Figure 4H,I.
#   (원본은 data.train/test.txt·SVM-RFE∩LASSO 를 씀. 논문 본문엔 SVM-RFE 없음 → 둘 다 산출, 논문을 정본으로.)
#
#   출력: results/step4_paper/. 원본 데이터·스크립트 무수정.
# ---------------------------------------------------------------------------

suppressMessages({
  library(glmnet); library(pROC); library(caret); library(kernlab); library(e1071)
  library(randomForest); library(xgboost)
})
library(here)
source(here::here("config.R"))
OUT <- file.path(RES_DIR, "step4_paper"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(123)

CAND10  <- c("ALDH2","FN1","VNN2","CREB5","XAF1","CA2","CDKN1B","IFI44L","SYTL2","TSPYL5")
EXPECT6 <- c("CDKN1B","ALDH2","FN1","XAF1","TSPYL5","VNN2")   # 논문 LASSO 6개
writeLines(CAND10, file.path(OUT, "interGenes.List.txt"))

read_lab <- function(p) { m <- as.matrix(read.table(p, header=TRUE, sep="\t", check.names=FALSE, row.names=1)); storage.mode(m)<-"double"; m }

## ---- 입력: 훈련 GSE96804 / 검증 STEP3 data.valid.paper.txt ----
TRAIN <- read_lab(file.path(OUT_DIR, "GSE96804.labeled.txt"))
vpath <- file.path(OUT_DIR, "data.valid.paper.txt")
stopifnot("STEP3 data.valid.paper.txt 없음 — STEP3 먼저 실행" = file.exists(vpath))
VALID <- read_lab(vpath)
ytr <- ifelse(grepl("_Control$", colnames(TRAIN)), 0, 1)
yva <- ifelse(grepl("_Control$", colnames(VALID)), 0, 1)
message("[train] GSE96804 ", nrow(TRAIN), "x", ncol(TRAIN), " (C", sum(ytr==0), "/D", sum(ytr==1), ")")
message("[valid] data.valid.paper ", nrow(VALID), "x", ncol(VALID), " (C", sum(yva==0), "/D", sum(yva==1), ")")

feats <- intersect(CAND10, rownames(TRAIN))
cat("\n== 후보 10 존재 (train/valid) ==\n")
for(g in CAND10) cat(sprintf("  %-8s train=%s valid=%s\n", g, g%in%rownames(TRAIN), g%in%rownames(VALID)))

## ---- (1) LASSO (원본 machine learning modeling 1.R 와 동일) ----
x  <- as.matrix(t(TRAIN[feats,]))
yf <- factor(ifelse(ytr==0,"Control","DKD"), levels=c("Control","DKD"))
set.seed(123)
fit   <- glmnet(x, yf, family="binomial", alpha=1)
cvfit <- cv.glmnet(x, yf, family="binomial", alpha=1, type.measure="deviance", nfolds=10)
co    <- coef(fit, s=cvfit$lambda.min)
lassoGene <- setdiff(rownames(co)[which(as.numeric(co)!=0)], "(Intercept)")
writeLines(lassoGene, file.path(OUT,"LASSO.gene.txt"))
cat("\n== LASSO ==\n  선택(",length(lassoGene),"):", paste(sort(lassoGene),collapse=", "),
    "\n  논문 6개 일치:", setequal(lassoGene, EXPECT6),
    "| 논문-우리:", paste(setdiff(EXPECT6,lassoGene),collapse=","),
    "| 우리-논문:", paste(setdiff(lassoGene,EXPECT6),collapse=","), "\n")

## ---- (2) SVM-RFE (원본 machine learning modeling 1.R 와 동일; methods→method 오타 정렬) ----
grp <- factor(ifelse(ytr==0,"Control","DKD"), levels=c("Control","DKD"))
set.seed(123)
svmGene <- tryCatch({
  Profile <- rfe(x = t(TRAIN[feats,]),
                 y = as.numeric(as.factor(grp)),              # 원본과 동일(회귀형 RFE)
                 sizes = c(2,3,4,5,6,7,8),
                 rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
                 method = "svmRadial",
                 preProcess = c("center","scale"))
  Profile$optVariables
}, error = function(e){ message("[SVM-RFE 실패] ", conditionMessage(e)); feats })
writeLines(svmGene, file.path(OUT,"SVM-RFE.gene.txt"))
cat("== SVM-RFE ==\n  선택(",length(svmGene),"):", paste(sort(svmGene),collapse=", "), "\n")

## ---- (3) Venn: LASSO ∩ SVM-RFE (원본과 동일) ----
interGenes <- intersect(lassoGene, svmGene)
writeLines(interGenes, file.path(OUT,"interGenes.txt"))
cat("== 교집합(LASSO ∩ SVM-RFE) ==\n  (",length(interGenes),"):", paste(sort(interGenes),collapse=", "),
    "\n  FN1:", "FN1"%in%interGenes, " ALDH2:", "ALDH2"%in%interGenes, "\n")

## ---- (4) 단일유전자 ROC + ci.auc(bootstrap) (원본 machine learning modeling 2.R 와 동일) ----
roc_ci <- function(mat, y, g){
  if(!(g %in% rownames(mat))) return(c(auc=NA, lo=NA, hi=NA))
  r <- roc(y, as.numeric(mat[g,]), quiet=TRUE)
  ci <- as.numeric(ci.auc(r, method="bootstrap"))   # ci[1]=lo, ci[2]=auc, ci[3]=hi
  c(auc=as.numeric(auc(r)), lo=ci[1], hi=ci[3])
}
set.seed(123)
rows <- lapply(CAND10, function(g){
  tr <- roc_ci(TRAIN,ytr,g); va <- roc_ci(VALID,yva,g)
  data.frame(gene=g, AUC_train=round(tr["auc"],3), tr_lo=round(tr["lo"],3), tr_hi=round(tr["hi"],3),
             AUC_valid=round(va["auc"],3), va_lo=round(va["lo"],3), va_hi=round(va["hi"],3), row.names=NULL)
})
roc_tab <- do.call(rbind, rows)
write.csv(roc_tab, file.path(OUT,"single_gene_ROC_withCI.csv"), row.names=FALSE)
write.csv(roc_tab[,c("gene","AUC_train","AUC_valid")], file.path(OUT,"single_gene_ROC_AUC.csv"), row.names=FALSE)
cat("\n== 단일유전자 ROC-AUC (train/valid, 95%CI) ==\n"); print(roc_tab, row.names=FALSE)

## 논문 본문 선택기준: AUC>0.8 in both train & valid
pass08 <- roc_tab$gene[roc_tab$AUC_train>0.8 & roc_tab$AUC_valid>0.8]
cat("\n  [논문기준] AUC>0.8 (train&valid 모두):", paste(pass08,collapse=", "), "\n")
writeLines(pass08, file.path(OUT,"genes_AUC_over_0.8_both.txt"))

## ---- (5) FN1+ALDH2 결합모델 (GLM/RF/SVM/XGBoost) = Figure 4H,I ----
dtr <- data.frame(FN1=TRAIN["FN1",], ALDH2=TRAIN["ALDH2",], y=ytr)
dva <- data.frame(FN1=VALID["FN1",], ALDH2=VALID["ALDH2",], y=yva)
comb <- function(tp, vp) c(train=as.numeric(auc(roc(dtr$y, tp, quiet=TRUE))),
                           valid=as.numeric(auc(roc(dva$y, vp, quiet=TRUE))))
res <- list()
g <- glm(y~FN1+ALDH2, data=dtr, family=binomial)
res$GLM <- comb(predict(g,dtr,type="response"), predict(g,dva,type="response"))
set.seed(123); rf <- randomForest(x=dtr[,c("FN1","ALDH2")], y=factor(dtr$y), ntree=500)
res$RF <- comb(predict(rf,dtr[,c("FN1","ALDH2")],type="prob")[,2], predict(rf,dva[,c("FN1","ALDH2")],type="prob")[,2])
set.seed(123); sv <- svm(x=dtr[,c("FN1","ALDH2")], y=factor(dtr$y), probability=TRUE, kernel="radial")
pt <- attr(predict(sv,dtr[,c("FN1","ALDH2")],probability=TRUE),"probabilities")[,"1"]
pv <- attr(predict(sv,dva[,c("FN1","ALDH2")],probability=TRUE),"probabilities")[,"1"]
res$SVM <- comb(pt, pv)
set.seed(123)
dtrain <- xgb.DMatrix(as.matrix(dtr[,c("FN1","ALDH2")]), label=dtr$y)
xgb <- xgb.train(list(objective="binary:logistic", max_depth=3, eta=0.3), dtrain, nrounds=50, verbose=0)
res$XGBoost <- comb(predict(xgb,as.matrix(dtr[,c("FN1","ALDH2")])), predict(xgb,as.matrix(dva[,c("FN1","ALDH2")])))
model_tab <- data.frame(model=names(res),
  AUC_train=round(sapply(res,function(v)v["train"]),3),
  AUC_valid=round(sapply(res,function(v)v["valid"]),3), row.names=NULL)
write.csv(model_tab, file.path(OUT,"combined_model_AUC.csv"), row.names=FALSE)
cat("\n== FN1+ALDH2 결합모델 AUC ==\n"); print(model_tab, row.names=FALSE)

## ---- (6) 논문값 대조 CSV ----
paperAUC <- data.frame(
  gene = c("FN1","ALDH2"),
  paper_train = c(0.911, 0.912), paper_valid = c(0.911, 0.815))
cmp <- merge(paperAUC, roc_tab[,c("gene","AUC_train","AUC_valid")], by="gene")
names(cmp)[names(cmp)=="AUC_train"] <- "ours_train"; names(cmp)[names(cmp)=="AUC_valid"] <- "ours_valid"
cmp$d_train <- round(cmp$ours_train - cmp$paper_train, 3)
cmp$d_valid <- round(cmp$ours_valid - cmp$paper_valid, 3)
write.csv(cmp, file.path(OUT,"compare_paper_vs_ours.AUC.csv"), row.names=FALSE)
cat("\n== 논문 대조 (FN1/ALDH2 AUC) ==\n"); print(cmp, row.names=FALSE)
message("\n[STEP4] 완료 -> ", OUT)
