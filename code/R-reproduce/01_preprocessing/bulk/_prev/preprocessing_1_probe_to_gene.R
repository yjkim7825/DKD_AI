
####----Basic_oprations----####
rm(list = ls())
gc()
getwd()
stringsASFactors = FALSE

##加载 R 包，根据理解，按照逻辑顺序逐个加载
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(limma)
library(stringr)

setwd("G:\\187geneMR\\03.download\\GSE37263")
####----001_GSE7476据的加载和数据转换判断----####
##读取表达矩阵
# Set a larger connection buffer size (in bytes)
Sys.setenv("VROOM_CONNECTION_SIZE" = 262144)  # Double the default
exprSet <- read_tsv("GSE37263_series_matrix.txt.gz", skip = 67) #根据实际情况修改
exprSet <- data.frame(exprSet)
class(exprSet)

#获得行为探针ID名称，列为样本名称的表达矩阵。
rownames(exprSet) <- exprSet$ID_REF
exprSet <- exprSet[,-1] 

####----003_probe_ID_translation_GPL570----####
library(readr)
library(dplyr)
probe <- read_tsv("GPL5175-3188.txt", skip = 12)  #根据实际情况修改
Probe_ID <- select(probe,c("ID","Gene Symbol"))  #探针文件改为

ids <- Probe_ID %>%
  mutate(`Gene Symbol` = str_split(`Gene Symbol`, " // ", simplify = TRUE)[, 2])
#保留多基因注释，但只取第一个基因作为代表,注意修改//,或者位置1
####----004_probe_ID_to_Expression_GPL570----####
library(dplyr)
library(tidyverse)

exprSet1 <- exprSet %>%
  rownames_to_column("ID")

ids$ID <- as.character(ids$ID)
exprset1 <- inner_join(ids,exprSet1,by="ID") 

length(exprset1$ID)
length(exprset1$'Gene Symbol')

exprset1 = avereps(exprset1[,-c(1,2)],      ##高频操作     
                   ID = exprset1$'Gene Symbol')
exprset1 <- data.frame(exprset1)  #获得行为基因名称，列为样本名称的表达举证。

# 确保 exprSet 是一个数据框，且行名是基因名（geneNames）
geneMatrix <- exprset1
geneMatrix <- cbind(geneNames = rownames(geneMatrix), geneMatrix)  # 添加 geneNames 列

# 导出为 TXT 文件（制表符分隔）
write.table(
  geneMatrix,
  file = "geneMatrix3.txt",
  sep = "\t",          # 使用制表符（TAB）分隔
  quote = FALSE,       # 不加引号
  row.names = FALSE,   # 不额外写入行名（因为已经作为第一列）
  col.names = TRUE     # 写入列名
)

