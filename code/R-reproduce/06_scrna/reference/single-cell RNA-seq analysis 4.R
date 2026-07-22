#这部分内容主要是教大家如何展示单细胞测序数据集中组别间的差异，相当于Seurat那里的“组间差异分析及可视化”
#这里的可视化适用于组间细胞类群比较类似的样本
#这部分的内容主要分析目的为：
#1、细胞通讯在不同组别中是否发生变化
#2、细胞通讯在不同细胞类型中是否发生变化
#3、细胞通讯的“发起者”与接收者是否随组别发生变化
setwd("G:\\187geneMR\\24.GSE209781\\7.cellchat")
#library(devtools)
#devtools::install_github("jokergoo/ComplexHeatmap")
#devtools::install_github("LTLA/BiocNeighbors")

#或者#if (!requireNamespace(“BiocManager”, quietly=TRUE))
# install.packages(“BiocManager”)
#BiocManager::install("ComplexHeatmap")
#ComplexHeatmap安装之后,cellchat才能安装成功
#devtools::install_github("sqjin/CellChat")

#if (!require("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
#BiocManager::install("Seurat")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(BiocNeighbors)
library(devtools)

#本地安装包,小红书教程
#install.packages("F:\\R-4.5.1\\library\\CellChat-main.zip", repos = NULL, type = "source")
#install.packages("F:\\R-4.5.1\\library\\CellChat-main", repos = NULL, type = "source") 

library(CellChat)
library(dplyr)
library(igraph)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratObject)
library(Seurat)
dir.create('./comparison')
setwd('./comparison')
# 创建cellchat.control对象

pbmc = Late

####cellchat.control预处理##############
Control = pbmc[, pbmc@meta.data$group %in% c("Control")]

rownames(Control@assays$RNA@layers$counts) <- rownames(Control)
# 获取归一化数据（log-normalized data）

data.input <- GetAssayData(Control, assay = "RNA", layer = "counts")  # Use layer = "counts" for raw counts
# 获取metadata
meta <- Control@meta.data  
celltypes <- unique(meta$celltype)  # 查看各聚类细胞类型

# 创建cellchat对象
cellchat.control <- createCellChat(object = as.matrix(data.input), 
                           meta = meta, 
                           group.by = "celltype")

# 保存cellchat对象
save(cellchat.control, file = "cellchat.control.RDA")

# 设置默认细胞标识
cellchat.control <- setIdent(cellchat.control, ident.use = "celltype") 
levels(cellchat.control@idents) # 显示细胞标签的因子水平

groupSize <- as.numeric(table(cellchat.control@idents)) 
# number of cells in each cell group
groupSize


#载入数据库并开始计算

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#展示细胞组成的比例:
dplyr::glimpse(CellChatDB$interaction)#展示互作的记录
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#取出相应分类(Secreted Signaling)用作分析数据库,相当于取子集

#CellChatDB.use <- CellChatDB   # simply use the default CellChatDB
cellchat.control@DB <- CellChatDB.use  #将数据库内容载入cellchat对象中

#表达量预处理
cellchat.control <- subsetData(cellchat.control,features = NULL)#取出表达数据

#devtools::install_github('immunogenomics/presto')

# 1. 识别高表达基因（相当于Seurat的FindVariableFeatures）
cellchat.control <- identifyOverExpressedGenes(cellchat.control)#寻找高表达的基因#


# 2. 识别高表达的配体-受体对
cellchat.control <- identifyOverExpressedInteractions(cellchat.control)#寻找高表达的通路

# 3. 计算细胞间通信概率
cellchat.control <- computeCommunProb(cellchat.control, raw.use = T)
#默认计算方式为type = "truncatedMean",
#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， 

###也可以自己设参数
cellchat.control <- computeCommunProb(
  cellchat.control,
  raw.use = TRUE,     # 使用原始count数据（若已归一化建议设为FALSE）
  type = "truncatedMean",  # 默认方法（表达量在trim阈值以下的视为0）
  trim = 0.1,         # 调整截断阈值（默认10%的细胞需表达该基因）
  population.size = TRUE  # 考虑细胞群大小的影响
)

#去掉通讯数量很少的细胞
cellchat.control <- filterCommunication(cellchat.control, min.cells = 10)

#将细胞通讯预测结果以数据框的形式取出
df.net <- subsetCommunication(cellchat.control)

class(df.net)  #data.frame
View(df.net)

write.csv(df.net,'01.df.control.net.csv')
#df.net <- subsetCommunication(cellchat,slot.name = "netP")
##这种方式只取通路，数据结构更简单
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#source可以以细胞类型的名称定义，也可以按照细胞名称中的顺序以数值向量直接取
#指定输入与输出的细胞集群
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#指定通路提取
cellchat.control <- computeCommunProbPathway(cellchat.control)
#每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
cellchat.control <- aggregateNet(cellchat.control)
#计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向

groupSize <- as.numeric(table(cellchat.control@idents))

groupSize

save(cellchat.control, file = "cellchat.controlprocessed.RDA")

#######cellchat.Late预处理#########
LateDKD = pbmc[,pbmc@meta.data$group %in% c("Late_DKD")]
rownames(LateDKD@assays$RNA@layers$counts) <- rownames(LateDKD)

data.input <- GetAssayData(LateDKD, assay = "RNA", layer = "counts")  # Use layer = "counts" for raw counts

meta <- LateDKD@meta.data  #与meta结构类似

celltypes <- unique(meta$celltype)  # 查看各聚类细胞类型

#
# 获取metadata


cellchat.Late <- createCellChat(object = as.matrix(data.input), 
                                   meta = meta, 
                                   group.by = "celltype")
cellchat.Late

save(cellchat.Late,file="cellchat.Late.RDA")
# set "labels" as default cell identity
cellchat.Late <- setIdent(cellchat.Late, ident.use = "celltype") 
#Idents(pbmc) <- 'group'
levels(cellchat.Late@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat.Late@idents)) 
# number of cells in each cell group
groupSize



#载入数据库并开始计算

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

#展示细胞组成的比例:
dplyr::glimpse(CellChatDB$interaction)#展示互作的记录
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#取出相应分类(Secreted Signaling)用作分析数据库,相当于取子集

#CellChatDB.use <- CellChatDB   # simply use the default CellChatDB
cellchat.Late@DB <- CellChatDB.use  #将数据库内容载入cellchat对象中

#表达量预处理
cellchat.Late <- subsetData(cellchat.Late,features = NULL)#取出表达数据

# 1. 识别高表达基因（相当于Seurat的FindVariableFeatures）
cellchat.Late <- identifyOverExpressedGenes(cellchat.Late)#寻找高表达的基因#

# 2. 识别高表达的配体-受体对
cellchat.Late <- identifyOverExpressedInteractions(cellchat.Late)#寻找高表达的通路

# 3. 计算细胞间通信概率
cellchat.Late <- computeCommunProb(cellchat.Late, raw.use = T)
#默认计算方式为type = "truncatedMean",
#默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0，

#去掉通讯数量很少的细胞
cellchat.Late <- filterCommunication(cellchat.Late, min.cells = 10)


df.net <- subsetCommunication(cellchat.Late)
#将细胞通讯预测结果以数据框的形式取出
class(df.net)  #data.frame
View(df.net)


write.csv(df.net,'01.df.net_Late.csv')
#df.net <- subsetCommunication(cellchat,slot.name = "netP")
##这种方式只取通路，数据结构更简单
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#source可以以细胞类型的名称定义，也可以按照细胞名称中的顺序以数值向量直接取
#指定输入与输出的细胞集群
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))#指定通路提取
cellchat.Late <- computeCommunProbPathway(cellchat.Late)
#每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
cellchat.Late <- aggregateNet(cellchat.Late)
#计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向

groupSize <- as.numeric(table(cellchat.Late@idents))

groupSize

save(cellchat.Late, file = "cellchat.Lateprocessed.RDA")



#######加载合并两个数据集##########
load("G:\\187geneMR\\24.GSE209781\\7.cellchat\\comparison\\cellchat.Lateprocessed.RDA")
load("G:\\187geneMR\\24.GSE209781\\7.cellchat\\comparison\\cellchat.controlprocessed.RDA")

# 正确做法：Control放第一位
object.list <- list(Control = cellchat.control, LateDKD = cellchat.Late)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#相当于Seurat对象的merge
cellchat
save(cellchat, file = "merged_cellchat.control_Late.rDA")

####################可视化#########################

#最简单的展示，查看细胞互作的数量在不同条件下是否有差异,有改动group = c(1,2)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(2,1),size.text = 20)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(2,1), measure = "weight",size.text = 20)
pdf("1.细胞互作的数量和强度1.pdf",width=8,height=8)
print(gg1 + gg2)
dev.off()

#查看细胞通路在两组间的富集程度,comparison = c(2, 1)可以改颜色


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(2, 1),do.stat = TRUE,font.size = 12) +  # 保持默认字体大小
  theme(
    legend.text = element_text(size = 14)  # 调大图例文字（LateDKD/Control）
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5)  # 调大颜色方块大小（默认是2-3）
    )
  )
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(2, 1),do.stat = TRUE,font.size = 12) +  # 保持默认字体大小
  theme(
    legend.text = element_text(size = 14)  # 调大图例文字（LateDKD/Control）
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 5)  # 调大颜色方块大小（默认是2-3）
    )
  )
pdf("2.细胞通路在两组间的富集程度2.pdf",width=7,height=8)
print(gg1)  #+ gg2)
dev.off()
#font.size = 11通路名称大小



#以circle plot的形式展示第二个组别中相较于第一个组别细胞通讯发生的变化，红色为上调蓝色为下调
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_diffInteraction(cellchat, weight.scale = T)
#netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


#以上是直接展示二组“相减”的结果，当然你也可以直接将两组分开展示
#weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
#par(mfrow = c(1,2), xpd=TRUE)
#for (i in 1:length(object.list)) {
#  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
#}


#同理，用heatmap也可以进行展示展示第二个组别中相较于
#第一个组别细胞通讯发生的变化，红色为上调蓝色为下调
#即组二减去组一,红色为高
gg1 <- netVisual_heatmap(cellchat,font.size = 16,font.size.title = 17)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 16,font.size.title = 17)
pdf("3.heatmap组别间细胞通讯发生的变化.pdf",width=10,height=8)
print(gg1 + gg2)
dev.off()

#####只展示一个图
gg1 <- netVisual_heatmap(cellchat,font.size = 16,font.size.title = 17)
#gg2 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 16,font.size.title = 17)
pdf("3.heatmap组别间细胞通讯发生的变化.pdf",width=5,height=5)
print(gg1)
dev.off()
####展示特定通路， 老三样##################


### 只需修改这里的通路名 ###
pathways.show <- "GDF"  # <<< 只需修改这个变量


# 检查通路是否存在
if(!pathways.show %in% cellchat.control@netP$pathways | 
   !pathways.show %in% cellchat.Late@netP$pathways) {
  stop("错误：通路 '", pathways.show, "' 在至少一个组别中不存在\n",
       "Control组通路: ", paste(cellchat.control@netP$pathways, collapse = ", "), "\n",
       "Late组通路: ", paste(cellchat.Late@netP$pathways, collapse = ", "))
}

# 1. 圆圈图
pdf(paste0("4.圆圈图_", pathways.show, ".pdf"), width=10, height=6)
weight.max <- getMaxWeight(
  object.list = list(Late = cellchat.Late, Control = cellchat.control),
  slot.name = "netP",
  attribute = pathways.show
)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# 2. 热图（自动使用上方定义的pathways.show）
pdf(paste0("4.2热图_", pathways.show, ".pdf"), width=10, height=6)
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "in", names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

# 3. 弦图
pdf(paste0("4.3弦图_", pathways.show, ".pdf"), width=10, height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, 
                      layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# 4. 基因表达
pdf(paste0("4.4基因表达_", pathways.show, ".pdf"), width=10, height=6)
plotGeneExpression(cellchat, signaling = pathways.show, 
                   split.by = "datasets", colors.ggplot = TRUE)
dev.off()

########批量通路分析#######
# 获取两组共有通路
common_pathways <- intersect(cellchat.control@netP$pathways, 
                             cellchat.Late@netP$pathways)

# 排除不需要的通路（可选）
excluded_pathways <- c("PARs", "CXCL")  # 示例：排除这两条通路
common_pathways <- setdiff(common_pathways, excluded_pathways)

cat("即将分析以下通路：\n", paste(common_pathways, collapse = ", "), "\n")

# 设置输出目录（自动创建）
output_dir <- "CellChat_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# 循环绘制每种可视化
for (pathway in common_pathways) {
  cat("\n正在分析通路：", pathway, "\n")
  
  # 1. 圆圈图
  pdf(file.path(output_dir, paste0("1.Circle_", pathway, ".pdf")), width=10, height=6)
  weight.max <- getMaxWeight(
    object.list = list(Control=cellchat.control, LateDKD=cellchat.Late),
    slot.name = "netP",
    attribute = pathway
  )
  par(mfrow=c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], 
                        signaling = pathway,
                        layout = "circle",
                        edge.weight.max = weight.max[1],
                        edge.width.max = 10,
                        signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()
  
  # 2. 弦图
  pdf(file.path(output_dir, paste0("2.Chord_", pathway, ".pdf")), width=10, height=6)
  par(mfrow=c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], 
                        signaling = pathway,
                        layout = "chord",
                        signaling.name = paste(pathway, names(object.list)[i]))
  }
  dev.off()
  
  # 3. 基因表达图
  pdf(file.path(output_dir, paste0("3.GeneExpression_", pathway, ".pdf")), 
      width=10, height=6)
  print(
    plotGeneExpression(cellchat, 
                       signaling = pathway,
                       split.by = "datasets",
                       colors.ggplot = TRUE)
  )
  dev.off()
  
  cat("已完成通路：", pathway, "的可视化\n")
}

cat("\n所有分析完成！结果已保存至：", normalizePath(output_dir))


#######气泡图########
# 加载必要的包
library(CellChat)
library(patchwork) # 用于图形组合

# 基础气泡图：显示所有细胞类型(1-12)之间的通讯差异
table(cellchat@meta[["celltype"]])

T      PCT       EC Mono-Mac     dPCT      MES  LOH-DCT        B     Mast 
5127     4940     2535     1902     1201     1106     1051      491      219 
Neut  Plasma1  Plasma2 
121       91       32 

gg1 <- netVisual_bubble(cellchat, 
                        sources.use = c(5), 
                        targets.use = c(1:12),  
                        comparison = c(1, 2), 
                        max.dataset = 2,         # 高亮第2组(LateDKD)增强的信号
                        title.name = "Increased signaling in LateDKD (vs Control)", # 修正标题
                        remove.isolate = TRUE,
                        thresh = 0.001) +           # 通信概率阈值（默认0.05）) +
  theme(
    # 调整坐标轴字体
    axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = 1, vjust = 0.5),  # X轴标签（靶细胞）
    axis.text.y = element_text(color = "black", size = 10),  # Y轴标签（配体-受体对）
    # 调整图例字体和标块大小
    legend.text = element_text(size = 10),    # 图例文字大小
    legend.title = element_text(size = 11))
pdf("dPCT-increased signaling in LateDKD.pdf", width=12, height=7.5)
print(gg1)  # 使用patchwork包合并两图
dev.off()
# Control组增强的信号（蓝色高亮）
gg2 <- netVisual_bubble(cellchat, 
                        sources.use = 1, 
                        targets.use = c(1:12),  
                        comparison = c(1, 2), 
                        max.dataset = 1,         # 高亮第1组(Control)增强的信号
                        title.name = "Increased signaling in Control (vs LateDKD)", # 修正标题
                        remove.isolate = TRUE)+
  theme(
    # 调整坐标轴字体
    axis.text.x = element_text(color = "black", size = 10, angle = 90, hjust = 1, vjust = 0.5),  # X轴标签（靶细胞）
    axis.text.y = element_text(color = "black", size = 10),  # Y轴标签（配体-受体对）
    # 调整图例字体和标块大小
    legend.text = element_text(size = 10),    # 图例文字大小
    legend.title = element_text(size = 11))
pdf("Increased signaling in Control.pdf", width=6, height=6)
print(gg2)  # 使用patchwork包合并两图
dev.off()


######下面的没用#########
##### 以上的计算依赖的都是细胞通讯的可能性(强度)，
#接下来我们将通过配受体对基因表达的角度进一步研究
#以通路为单位
pos.dataset = "Control"
features.name = pos.dataset
#差异计算：
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 1)

#提取细胞通讯预测数据
net <- netMappingDEG(cellchat, features.name = features.name)

#提取Control组上调的配体以及non.Control上调的受体，
net.up <- subsetCommunication(cellchat, net = net, datasets = "Control",
                              ligand.logFC = 0.2, receptor.logFC = NULL)

#反之亦然
net.down <- subsetCommunication(cellchat, net = net, datasets = "Control",
                                ligand.logFC = -0.2, receptor.logFC = -0.1)

#提取其中的差异基因
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


#可视化：
#气泡图
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = 1, targets.use = c(2:7), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", 
                                            names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use = 1, targets.use = c(2:7), 
                        comparison = c(1, 2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", 
                                            names(object.list)[2]))
pdf("6.Control中上调下调的通路.pdf",width=10,height=6)
print(gg1 + gg2)
dev.off()

#弦图
par(mfrow = c(1,2), xpd=TRUE)

netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), 
                     slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", 
                                         names(object.list)[1]))
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), 
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", 
                                         names(object.list)[2]))



#绘制基因表达量

cellchat@meta$datasets = factor(cellchat@meta$datasets, 
                                levels = c("LateDKD", "Control")#设定组别的level
                                )

cellchat@meta$datasets = factor(cellchat@meta$datasets, 
                                levels = c("Control","LateDKD")#设定组别的level
)
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)

#版本信息
sessionInfo()
