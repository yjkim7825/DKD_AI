
#install.packages("pROC")


library(pROC)                  #引用包
expFile="data.train.txt"      #表达数据文件
geneFile="LASSO.gene.txt"      #交集基因列表文件
setwd("G:\\187geneMR\\20.lasso_SVM\\16.ROC\\GSE96804_104948")    #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
#rt=t(rt)


# 提取样本分组（Control=0, DKD=1）
y <- ifelse(grepl("_Control$", colnames(rt)), 0, 1)

#或者y = gsub(".*_(Control|DKD)$", "\\1", colnames(rt))
#y=ifelse(y=="Control", 0, 1)
#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
for(x in as.vector(geneRT[,1])){
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("training_ROC.",x,".pdf"), width=3.5, height=3.5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(x, " of Training")) 
	text(0.34, 0.38, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}

#########验证集###

library(pROC)                  #引用包
expFile="data.test.txt"      #表达数据文件
geneFile="LASSO.gene.txt"      #交集基因列表文件
#setwd("G:\\187geneMR\\20.lasso_SVM\\16.ROC")    #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
#rt=t(rt)


# 提取样本分组（Control=0, DKD=1）
y <- ifelse(grepl("_Control$", colnames(rt)), 0, 1)

#或者y = gsub(".*_(Control|DKD)$", "\\1", colnames(rt))
#y=ifelse(y=="Control", 0, 1)
#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
for(x in as.vector(geneRT[,1])){
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("test_ROC.",x,".pdf"), width=3.5, height=3.5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(x, " of Validation")) 
	text(0.34, 0.38, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}


########tiff#########
#install.packages("pROC")
library(pROC)                  #引用包
expFile="data.train.txt"      #表达数据文件
geneFile="LASSO.gene.txt"      #交集基因列表文件
setwd("G:\\187geneMR\\20.lasso_SVM\\16.ROC\\GSE96804_104948")    #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
#rt=t(rt)

# 提取样本分组（Control=0, DKD=1）
y <- ifelse(grepl("_Control$", colnames(rt)), 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
for(x in as.vector(geneRT[,1])){
  #绘制ROC曲线
  roc1=roc(y, as.numeric(rt[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  tiff(file=paste0("training_ROC.",x,".tiff"), width=3.5, height=3.5, units="in", res=300)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(x, " of Training")) 
  text(0.34, 0.38, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}

#########验证集###
library(pROC)                  #引用包
expFile="data.test.txt"      #表达数据文件
geneFile="LASSO.gene.txt"      #交集基因列表文件

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

# 提取样本分组（Control=0, DKD=1）
y <- ifelse(grepl("_Control$", colnames(rt)), 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
for(x in as.vector(geneRT[,1])){
  #绘制ROC曲线
  roc1=roc(y, as.numeric(rt[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  tiff(file=paste0("test_ROC.",x,".tiff"), width=3.5, height=3.5, units="in", res=300)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=paste0(x, " of Validation")) 
  text(0.34, 0.38, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
