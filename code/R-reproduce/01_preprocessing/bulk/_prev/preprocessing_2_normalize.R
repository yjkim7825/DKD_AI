
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


# 引用包
library(limma)

# 设置文件路径
inputFile <- "GSE142025_matrix.txt"   # 表达数据文件
controlFile <- "s1.txt"               # 对照组样品信息
earlyFile <- "s2.txt"                 # 早期组样品信息
lateFile <- "s3.txt"                  # 晚期组样品信息
geoID <- "GSE142025"
setwd("G:\\187geneMR\\05.normalize")

# 读取并预处理表达矩阵
rt <- read.table(inputFile, header=TRUE, sep="\t", check.names=FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt <- avereps(data)

# 自动log2转换判断
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)
                 
                 # 读取三组样品信息
                 sample_control <- read.table(controlFile, header=FALSE, sep="\t", check.names=FALSE)
                 sample_early <- read.table(earlyFile, header=FALSE, sep="\t", check.names=FALSE)
                 sample_late <- read.table(lateFile, header=FALSE, sep="\t", check.names=FALSE)
                 
                 # 清理样品名称
                 sampleName_control <- gsub("^ | $", "", as.vector(sample_control[,1]))
                 sampleName_early <- gsub("^ | $", "", as.vector(sample_early[,1]))
                 sampleName_late <- gsub("^ | $", "", as.vector(sample_late[,1]))
                 
                 # 提取三组数据
                 controlData <- data[, sampleName_control]
                 earlyData <- data[, sampleName_early]
                 lateData <- data[, sampleName_late]
                 
                 # 合并三组数据
                 combinedData <- cbind(controlData, earlyData, lateData)
                 controlNum <- ncol(controlData)
                 earlyNum <- ncol(earlyData)
                 lateNum <- ncol(lateData)
                 
 # 输出所有基因矫正后的表达量（带三组标签）
Type <- c(rep("Control", controlNum), rep("Early", earlyNum), rep("Late", lateNum))
outData <- rbind(id=paste0(colnames(combinedData), "_", Type), combinedData)
write.table(outData,file=paste0(geoID, "_threeGroups.normalize.txt"), sep="\t",quote=F,col.names=F)

 