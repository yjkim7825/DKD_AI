library(limma)
library(dplyr)
library(pheatmap)
library(ggplot2)

logFCfilter = 0.585 # logFC的过滤条件
adj.P.Val.Filter = 0.05 # 矫正后p值的过滤条件
inputFile = "GSE142025_twoGroups.normalize.txt" # 表达数据文件
setwd("G:\\187geneMR\\08.Diffgene\\GSE142025_twogroup") # 设置工作目录

#读取输入文件，并对输入文件整理
rt = read.table(inputFile, header=T, sep="\t", check.names=F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data = avereps(data)

#获取样品的分组信息（修改为只有Control和DKD两组）
Type = sapply(strsplit(colnames(data), "_"), function(x) x[2])

#只保留Control和DKD样本（假设DKD样本标记为"DKD"）
keep_samples = Type %in% c("Control", "DKD")
data = data[, keep_samples]
Type = Type[keep_samples]
data = data[, order(Type)] # 根据样品的分组信息对样品进行排序

#检查分组
print(table(Type))

#进行差异分析
design <- model.matrix(~0 + factor(Type))
colnames(design) <- levels(factor(Type))
fit <- lmFit(data, design)

#设置对比矩阵（Control vs DKD）
cont.matrix <- makeContrasts(DKD - Control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#输出所有基因的差异情况
allDiff = topTable(fit2, adjust = 'fdr', number = 200000)
allDiffOut = rbind(id = colnames(allDiff), allDiff)
write.table(allDiffOut, file = "all_DKD_vs_Control.txt",
            sep = "\t", quote = F, col.names = F)

#输出显著的差异基因
diffSig = allDiff[with(allDiff, (abs(logFC) > logFCfilter & adj.P.Val < adj.P.Val.Filter)), ]
if (nrow(diffSig) > 0) {
  diffSigOut = rbind(id = colnames(diffSig), diffSig)
  write.table(diffSigOut, file = "diff_DKD_vs_Control.txt",
              sep = "\t", quote = F, col.names = F)
  
  #输出差异基因表达量
  diffGeneExp = data[row.names(diffSig), ]
  diffGeneExpOut = rbind(id = paste0(colnames(diffGeneExp), "_", Type), diffGeneExp)
  write.table(diffGeneExpOut, file = "diffGeneExp_DKD_vs_Control.txt",
              sep = "\t", quote = F, col.names = F)
  
 # 绘制差异基因热图
  geneNum = 50 # 定义展示基因的数目
  diffUp = diffSig[diffSig$logFC > 0, ]
  diffDown = diffSig[diffSig$logFC < 0, ]
  geneUp = row.names(diffUp)
  geneDown = row.names(diffDown)
  if (nrow(diffUp) > geneNum) { geneUp = row.names(diffUp)[1:geneNum] }
  if (nrow(diffDown) > geneNum) { geneDown = row.names(diffDown)[1:geneNum] }
  hmExp = data[c(geneUp, geneDown), ]
  
 # 准备注释文件
  annotation_col = data.frame(Group = Type)
  rownames(annotation_col) = colnames(data)
  
  #输出图形
  pdf(file = "heatmap_DKD_vs_Control.pdf", width = 10, height = 7)
  pheatmap(hmExp,
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue2", "white", "red2"))(50),
           cluster_cols = FALSE,
           show_colnames = FALSE,
           scale = "row",
           fontsize = 8,
           fontsize_row = 5.5,
           fontsize_col = 8)
  dev.off()
  
  #绘制火山图
  rt = allDiff
  Sig = ifelse((rt$adj.P.Val < adj.P.Val.Filter) & (abs(rt$logFC) > logFCfilter),
               ifelse(rt$logFC > logFCfilter, "Up", "Down"), "Not")
  rt = mutate(rt, Sig = Sig)
  p = ggplot(rt, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = Sig)) +
    scale_color_manual(values = c("blue2", "grey", "red2")) +
    labs(title = "DKD vs Control") +
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
  
  #输出图形
  pdf(file = "vol_DKD_vs_Control.pdf", width = 5.5, height = 4.5)
  print(p)
  dev.off()
} else {
  message("No significant genes found in DKD vs Control")
}

#由于只有两组比较，不需要Venn图部分,如果需要保存所有差异基因（上调和下调）
if (exists("diffSig") && nrow(diffSig) > 0) {
  
  #上调基因
  upGenes = rownames(diffSig[diffSig$logFC > 0, ])
  write.table(upGenes, file = "up_genes_DKD_vs_Control.txt",
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  #下调基因
  downGenes = rownames(diffSig[diffSig$logFC < 0, ])
  write.table(downGenes, file = "down_genes_DKD_vs_Control.txt",
              sep = "\t", quote = F, row.names = F, col.names = F)
}
