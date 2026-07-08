################# lasso ###############
#install.packages("glmnet")
set.seed(123)
library(glmnet)                   

# 设置路径
inputFile <- "data.train.txt"       
interGeneFile <- "interGenes.List.txt"  # 新增：基因列表文件
setwd("G:\\187geneMR\\20.lasso_SVM\\12.lasso\\GSE96804_104948")       

# 读取输入文件并筛选目标基因
rt <- read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
interGenes <- scan(interGeneFile, what="character", quiet=TRUE)  # 读取基因列表
rt <- rt[rownames(rt) %in% interGenes, ]  # 筛选目标基因

rt=t(rt)

#构建模型
x=as.matrix(rt)
#y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))  #从行名提取样本类别（如"Tumor"或"Normal"）
y = gsub(".*_(Control|DKD)$", "\\1", row.names(rt))
y = factor(y)  # 确保y是因子类型
fit=glmnet(x, y, family = "binomial", alpha=1)
#family = "binomial"：指定逻辑回归（因变量为二分类）。
#alpha = 1：启用 L1 正则化（即 LASSO，若 alpha=0 为岭回归，alpha∈(0,1) 为弹性网络）
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
#10 折交叉验证：评估不同 λ 下模型的预测偏差（type.measure='deviance'）

pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#输出筛选的特征基因
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)


################SVM-RFE###########
#引用包
library(e1071)
library(kernlab)
library(caret)

set.seed(123)
inputFile="data.train.txt"        #输入文件

# 设置路径
inputFile <- "data.train.txt"       
interGeneFile <- "interGenes.List.txt"  # 新增：基因列表文件
    
# 读取输入文件并筛选目标基因
rt <- read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
interGenes <- scan(interGeneFile, what="character", quiet=TRUE)  # 读取基因列表
rt <- rt[rownames(rt) %in% interGenes, ]  # 筛选目标基因
rt=t(rt)

#构建模型
x=as.matrix(rt)
#y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))  #从行名提取样本类别（如"Tumor"或"Normal"）
group = gsub(".*_(Control|DKD)$", "\\1", row.names(rt))

# 确保group是因子类型且水平正确
group = factor(group, levels=c("Control", "DKD"))
#SVM-RFE分析
Profile=rfe(x=rt,
            y=as.numeric(as.factor(group)),
            #sizes = c(2,4,6,8, seq(10,40,by=3)),
            sizes = c(2, 3, 4, 5, 6, 7, 8),  # 最多测试 8 个基因
            #sizes = seq(2, 10)  # 测试 2~10 个基因
            #通常建议特征数 ≤ 样本量的 1/10（即 ≤ 3-4 个基因），但 SVM 可以稍放宽至 5-8 个
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial",
            preProcess = c("center", "scale"))  # 新增此行

#绘制图形
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#标注交叉验证误差最小的点
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

#输出选择的基因
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

###############venn##################

# 加载包
library(VennDiagram)
# 设置参数
outFile <- "interGenes.txt"        # 输出交集基因文件
# 初始化基因列表
geneList <- list()

# 1. 读取LASSO回归的基因
rt <- read.table("LASSO.gene.txt", header=F, sep="\t", check.names=F)
geneList[["LASSO"]] <- unique(as.vector(rt[,1]))  # 直接提取并去重

# 2. 读取SVM-RFE的基因
rt <- read.table("SVM-RFE.gene.txt", header=F, sep="\t", check.names=F)
geneList[["SVM-RFE"]] <- unique(as.vector(rt[,1]))

# 3. 绘制Venn图（完全保持原风格）
venn.plot <- venn.diagram(
  x = geneList,
  filename = NULL,  # 不自动保存文件
  fill = c("cornflowerblue", "darkorchid1"),  # 指定填充色
  alpha = 0.5,      # 透明度
  scaled = FALSE,   # 不按比例缩放
  cat.pos = c(-20, 20),  # 标签位置（-20和20度）
  cat.dist = 0.05,  # 标签距离圆圈的距离
  cat.col = c("cornflowerblue", "darkorchid1"),  # 标签颜色
  cat.cex = 1.2,    # 标签字体大小
  margin = 0.1,     # 图形边距
  cex = 1.2         # 交集数字大小
)

# 保存为PDF
pdf(file = "venn.pdf", width = 5, height = 5)
grid.draw(venn.plot)  # 绘制图形
dev.off()

# 4. 保存交集基因
intersectGenes <- Reduce(intersect, geneList)
write.table(
  intersectGenes, 
  file = outFile, 
  sep = "\t", 
  quote = FALSE, 
  col.names = FALSE, 
  row.names = FALSE
)
