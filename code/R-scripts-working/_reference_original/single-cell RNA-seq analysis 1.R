setwd("G:\\187geneMR\\24.GSE209781\\1.Seurat\\6 sample")    
###################################04.04.数据前期处理和矫正###################################
#GITR (TNFRSF18) IL2RA仅存在/ILRα/Tac/IL2R/p55/ (CD25); IL7R (CD127); SELL (CD62L);
#"ITGAE"="CD103" ;TIM-3(CD366/HAVCR2); 
#读取文件，并对重复基因取均值
gc()
####################从这一部开始分析#####################
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(harmony)
library(cowplot)
library(patchwork)
library(devtools)
library(tidyverse)
#install.packages('devtools')
#devtools::install_github('immunogenomics/harmony')
#devtools::install_github('junjunlab/scRNAtoolVis')
#devtools::install_github('sajuukLyu/ggunchull')
#devtools::install_github('lydiaMyr/ImmuCellAI')
#devtools::install_github('mojaveazure/seurat-disk')
library(scRNAtoolVis) #为了clusterCornerAxes函数
library(ggunchull)
Late <- readRDS("G:\\187geneMR\\24.GSE209781\\1.Seurat\\6 sample\\late_merged_seurat.rda")

# 2. 合并数据层（关键步骤）
Late <- JoinLayers(Late)
# 确保assay的行名与counts层一致
rownames(Late@assays$RNA@layers$counts) <- rownames(Late)



Late[["percent.mt"]] <- PercentageFeatureSet(object = Late, pattern = "^MT-")
pdf(file="featureViolin.pdf",width=10,height=6)           #保存基因特征小提琴图
VlnPlot(object = Late, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

Late <- subset(Late,
               subset = nFeature_RNA > 300 &
                 nFeature_RNA < 5000 &
                 percent.mt < 10)
##美化,pt.size = 0
pdf(file="featureViolin2.pdf",width=10,height=6)
Late[["percent.mt"]] <- PercentageFeatureSet(object = Late, pattern = "^MT-")
VlnPlot(Late,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0,group.by = "orig.ident")&NoLegend()&labs(x = '') #&geom_hline(yintercept = 2500)
dev.off()



pdf(file="featureCor.pdf",width=10,height=6)             #保存基因特征相关性图
plot1 <- FeatureScatter(object = Late, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = Late, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

Late <- NormalizeData(object = Late, normalization.method = "LogNormalize", scale.factor = 10000)
Late <- FindVariableFeatures(object = Late, selection.method = "vst", nfeatures = 2000)

top10 <- head(x = VariableFeatures(object = Late), 10)  #在图片上标记出前十个
pdf(file="featureVar.pdf",width=10,height=6)                 #保存基因特征方差图
plot1 <- VariableFeaturePlot(object = Late)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()



table(Late$orig.ident)  #查看各类细胞数目
#Control1  Control2  Control3 Late_DKD1 Late_DKD2 Late_DKD3 
#3603      5383       949      1081      5188      2612 


###################################05.PCA主成分分析###################################
Late=ScaleData(Late)                   #PCA降维之前的标准预处理步骤
#Late=RunPCA(object= Late,npcs = 30,pc.genes=VariableFeatures(object = Late))     #PCA分析
Late <- RunPCA(Late,npcs = 30,verbose = FALSE)
#ElbowPlot(Late,ndims=50)
# 更健壮的分组方式（处理所有可能的Control和DKD样本）
Late@meta.data$patient <- ifelse(grepl("^Control", Late$orig.ident), 
                                 "Control",
                                 "Late_DKD")
Late@meta.data$group <- ifelse(grepl("^Control", Late$orig.ident), 
                                 "Control",
                                 "Late_DKD")
# 验证分组
table(Late@meta.data[["patient"]])


save(Late, file = "pbmc_after_pca.RDA")
#在object的metadata中定义细胞ID信息，变量名为stim,对应细胞个数
#Late@meta.data$stim <- c(rep("ccRCC", 20058), rep("NS", 8347))
options(repr.plot.height = 5, repr.plot.width = 12)

pdf(file="05.before harmony.pdf",width=6.5,height=6)
DimPlot(object = Late, reduction = "pca", pt.size = .1, group.by = "patient")
dev.off()

pdf(file="05.before harmony2.pdf",width=6.5,height=6)
VlnPlot(object = Late, features = "PC_1", group.by = "patient", pt.size = .1)
dev.off()

#
options(repr.plot.height = 2.5, repr.plot.width = 6)
Late <- Late %>% 
  RunHarmony("patient", plot_convergence = TRUE)

#获取Harmony 矫正之后的信息，使用Embeddings()函数
harmony_embeddings <- Embeddings(Late, 'harmony')
harmony_embeddings[1:5, 1:5]

##查看数据Harmony整合之后的前两个维度上数据是不是很好的整合，最好是很好的整合结果。
pdf(file="05.after harmony.pdf",width=12,height=5)
#options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = Late, reduction = "harmony", pt.size = .1, group.by = "patient")
p2 <- VlnPlot(object = Late, features = "harmony_1", group.by = "patient", pt.size = .1)
plot_grid(p1,p2)
dev.off()




#绘制每个PCA成分的相关基因,分为20个pc?#
pdf(file="051.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = Late, dims = 1:4, reduction = "pca",nfeatures = 20) 
#只显示出四个pc图#每个图展示20个基因
dev.off()

#主成分分析图形
pdf(file="05.PCA.pdf",width=6.5,height=6)
DimPlot(object = Late, reduction = "pca")
dev.off()

#主成分分析热图#只展示前四个热图#每张图显示30个基因#一行两个图#
pdf(file="05.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = Late, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#每个PC的p值分布和均匀分布
Late <- JackStraw(object = Late, num.replicate = 100)
Late <- ScoreJackStraw(object = Late, dims = 1:20)
pdf(file="05.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = Late, dims = 1:20)  #解螺旋是15,生信自学网20
dev.off()
save(Late, file = "pbmc_ScoreJackStraw_12636.rds")


Late <- readRDS("G:/187geneMR/24.GSE209781/1.Seurat/6 sample/pbmc_ScoreJackStraw_12636.rda")

###################################06.TSNE聚类分析和marker基因###################################
setwd("G:\\187geneMR\\24.GSE209781\\3.annotation\\6 sample")
####UMAP/tSNE聚类分析#综合解螺旋的代码#比较好用#
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)

pcSelect=30
Late <- FindNeighbors(Late,reduction = "harmony", dims = 1:pcSelect)       #计算邻接距离 #解螺旋取得值为10,生信自学网为20          
Late <- FindClusters(Late, reduction = "harmony",resolution = 0.2)        #对细胞分组,优化标准模块化 
####选择不同的resolution值可以获得不同的cluster数目，值越大cluster数目越多，默认值是0.5.         
Late <- RunTSNE(Late, reduction = "harmony",dims = 1:pcSelect)      #TSNE聚类 #解螺旋取得值为10,生信自学网为20               
#DimPlot(Late, reduction ="tsne")
TSNEPlot(Late, label = TRUE)
ggsave(filename = "06.TSNE2.pdf",width=6.5,height=6)
#write.table(Late$seurat_clusters,file="06.tsneCluster.txt",quote=F,sep="\t",col.names=F)
Late <- RunUMAP(Late,reduction = "harmony", dims = 1:pcSelect)    #解螺旋取得值为10,生信自学网为20     
#save(Late, file = "pbmc_for_14cluster_markers.RDA")


DimPlot(Late,reduction = "umap",label = TRUE)
 #UAMP可视化#标注序号# pt.size = 2为点的大小#pt.size = 2
ggsave(filename = "06.UMAP.pdf",width=6.5,height=6)
dev.off()

DimPlot(Late,pt.size = 1,reduction = "umap",group.by="orig.ident",label = TRUE,label.size = 3)
 #UAMP可视化#标注序号# pt.size = 2为点的大小#
ggsave(filename = "06.UMAP_sample.pdf",width=6.5,height=6)
dev.off()

# 寻找所有以聚类的差异基因`
cluster.markers <- FindAllMarkers(
  object = Late,
  min.pct = 0.25,
  logfc.threshold = 0.5,
  only.pos = TRUE,  # 只保留上调基因
  verbose = TRUE
)
significant.markers <- subset(cluster.markers, 
                              p_val_adj < 0.05 & abs(avg_log2FC) > 2)

# 3. 按cluster和log2FC排序输出
significant.markers <- significant.markers[order(
  significant.markers$cluster, 
  -abs(significant.markers$avg_log2FC)
), ]

sig50_markers <- subset(cluster.markers, p_val_adj < 0.05 & avg_log2FC > 1)

# 2. 提取每个cluster前50基因
top50 <- sig50_markers %>% 
  group_by(cluster) %>% 
  top_n(50, avg_log2FC) %>% 
  arrange(cluster, -avg_log2FC)

# 3. 输出到CSV
write.csv(top50, "top50_markers_per_cluster.csv", row.names = FALSE)
# 4. 输出完整结果和筛选结果
write.table(cluster.markers, 
            file = "all_cluster_markers.xls",
            sep = "\t", 
            row.names = TRUE, 
            quote = FALSE)

write.table(significant.markers,
            file = "significant_cluster_markers.xls",
            sep = "\t",
            row.names = TRUE,
            quote = FALSE)

# 9 kidney cell types (podocytes, fibroblasts, endothelial cells, proximal
#tubules, loop of Henle [LOH], distal convoluted tubules
# [DCT], DCT-LOH, collecting ducts [CD]–principal cells [PC],
# CD-intercalated cells [IC]) and 7 immune cell types (B cells, T
# cells, macrophages, monocytes, mast cells, neutrophils, and T
#natural killer cells) 
###################################07.注释细胞类型###############################
save(Late, file = "Late_for_annotation_18888.rds")

setwd("G:\\187geneMR\\24.GSE209781\\3.annotation\\6 sample")
pdf(file="02CD4_CD8.pdf",width=17,height=6)
#CD4_CD8
cluster4Marker=c("PDGFRB","FN1","COL4A2","NPHS1","PODXL","NPHS2","SLC4A11","AQP9","SLC34A1", "CUBN", "LRP2", "SLC22A6", "GPX3", "ALDOB", "SLC22A8","ALDH1A2","ALDH2") #uster4中logFC最正和最负值的五个基因
DotPlot(object = Late, features = cluster4Marker)
dev.off()


pdf(file="近端小管上皮细胞0,2,12 (PT).pdf",width=12,height=6)
cluster4Marker=c("COL1A1","COL3A1","DCN","RBP4", "AGXT","SLC22A8", "CYP4A11", "SLC34A1", "HAO2","ALDOB", "SLC22A6", "LRP2", "CUBN", "GPX3") #uster4中logFC最正和最负值的五个基因
DotPlot(object = Late, features = cluster4Marker)
dev.off()


pdf(file="髓袢升支粗段 (LOH).pdf",width=12,height=6)
cluster4Marker=c("UMOD", "SLC12A1", "CLDN16", "KCNJ1", "SLC12A3","CLCNKA", "CLCNKB")
DotPlot(object = Late, features = cluster4Marker)
dev.off()

pdf(file="集合管主细胞-远曲小管 (DCT)-.pdf",width=15,height=6)
cluster4Marker=c("AQP2","AQP3", "CALB1", "HSD11B2", "SCNN1G", "SCNN1B", "FXYD4", "RALBP1", "GATA3",
                 "SLC12A3", "WNK1", "KLHL3", "CNNM2", "TRPM6", "WNK4")
DotPlot(object = Late, features = cluster4Marker)
dev.off()


pdf(file="1损伤PCT.pdf",width=17,height=6)
cluster4Marker=c("FN1", "VCAM1", "CD74", "TYROBP", "CTSB", 
                                 "CD44", "HMOX1", "SERPINE1", "ANXA1", 
                                 "S100A11", "S100A10", "VIM", "LGALS1",
                                 "CHI3L1", "LCN2", "TIMP1", "B2M","HAVCR1")
DotPlot(object = Late, features = cluster4Marker)
dev.off()

pdf(file="2损伤PCT.pdf",width=17,height=6)
cluster4Marker=c("SLC4A11", "AQP9", "ALDH1A2", "LRP2", "CCN2", "VIM", 
                 "THBS1", "COL4A1", "COL4A2", "FNBP1", "PLOD2", "ITGB3", 
                 "ITGA3", "HAVCR1", "PDGFB", "ZEB2", "NFKB1", "C3", 
                 "CXCL1", "IL34", "PCNA", "PROM1", "EPCAM")
DotPlot(object = Late, features = cluster4Marker)
dev.off()




pdf(file="免疫细胞.pdf", width=25,height=6)
cluster4Marker=c("IL7R", "CD3E", "TRAC", "GZMK", "PTPRC",
"CD163", "C1QC", "LYZ", "APOE", "CSF1R", "ITGAX", "CD14", "TLR2", "TLR4", "IL1B",
"GNLY", "NKG7", "GZMB", "KLRD1",
"CD79A", "CD79B", "MS4A1", "IGKC", "BANK1",
"TPSAB1","TPSAB2", "CPA3", "KIT", "HDC",
"IGHG", "IGHA", "JCHAIN", "CD38",
"FCGR3B", "S100A8", "S100A9", "S100A12", "FPR1", "FPR2", "TLR1",  "CSF3R")
DotPlot(object = Late, features = cluster4Marker)
dev.off()
#"T细胞", "巨噬细胞", "NK-T细胞", "B细胞", "肥大细胞", "浆细胞", "中性粒细胞"

pdf(file="内皮_系膜_足_FN1_ALDH2.pdf", width=19,height=6)
cluster4Marker=c("VWF", "PECAM1", "FLT1", "PLVAP",
                 "NPHS1", "NPHS2", "PLA2R1", "MAFB", "SYNPO", "PODXL", "PLCE1",
                 "PDGFRB", "ITGA8","COL6A2" , "ITGA1", "MYL9","ACTA2", "COL1A1", "RGS5", "NOTCH3", "TAGLN", "LMOD1", "THBS1",
                 "FN1","ALDH2")

DotPlot(object = Late, features = cluster4Marker)
dev.off()

pdf(file="代谢活跃_PCT.pdf",width=12,height=6)
cluster4Marker=c("CYP4A11", "PCK1", "FABP1", "SLC5A12", "HAO2", 
                 "ACSM2A", "SLC13A3", "UGT1A8", "DPEP1", "SLC22A8",
                 "SLC4A11","ALDOB","LRP2", "CUBN","SLC34A1", "SLC22A6", "GPX3")
DotPlot(object = Late,features = cluster4Marker)+coord_flip()
dev.off()

#######分组展示
DimPlot(Late,
        pt.size = 1,
        reduction = "tsne",
        group.by = "seurat_clusters",
        split.by = "patient",  # 按group列拆分
        label = TRUE,
        label.size = 3) +
  plot_annotation(title = "Normal vs DKD (Clusters 0-16)")

ggsave("06.UMAP_Normal_vs_DKD_split.pdf", width = 12, height = 6)
dev.off()

library(Seurat)
library(patchwork)  # 用于拼图


#重新对注释后的细胞可视化
#Late <- subset(Late, idents = c("21"), invert = TRUE)#去掉低质量细胞群
new.cluster.ids <- c("0"="T", 
                     "1"="PCT", 
                     "2"="EC", 
                     "3"="Mono-Mac", 
                     "4"="dPCT", 
                     "5"="MES", 
                     "6"="LOH-DCT", 
                     "7"="B", 
                     "8"="Mast", 
                     "9"="Neut", 
                     "10"="Plasma1", 
                     "11"="Plasma2"
                    )
#save(Late, file = "pbmc_for_cell_markers.RDA") #按照上面注释,生成此文件


Late <- RenameIdents(Late, new.cluster.ids)                        
Late$celltype <- Late@active.ident #将注释结果添加到metadata
#也可以从这一步开始进行monocle分析,参见简书代码

cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
               "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272',"#fc8542") #,,
DimPlot(Late, reduction = "umap",group.by = "celltype",label = T,pt.size = 0.2,cols=cbPalette)#
#DimPlot(Early,pt.size = 0.2,reduction = "umap",label = TRUE)
#UAMP可视化#标注序号# pt.size = 2为点的大小#
ggsave(filename = "注释后.sample_cellUMAP.pdf",width=7.5,height=6)
dev.off()
TSNEPlot(Late, pt.size = 0.2, label = TRUE,cols=cbPalette)
ggsave(filename = "注释后.ample_cellTSNE.pdf",width=7.5,height=6)
dev.off()
save(Late, file = "Late_annotated.RDA")
#按照sample对cluster进行注释,改成late

#biomamama
## umap/tsne
library(RColorBrewer)
color_ct=c(brewer.pal(12, "Set3"),"#b3b3b3",
           brewer.pal(5, "Set1"),
           brewer.pal(3, "Dark2"),
           "#fc4e2a","#fb9a99","#f781bf","#e7298a")
clusterCornerAxes(Late,reduction = 'umap',clusterCol = 'celltype',pSize = 0.05,cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(cbPalette,0.65)) + NoLegend() +
  scale_fill_manual(values = alpha(cbPalette,0.65))
#pSize = 0.5为细胞点的大小;cellLabelSize = 5文字大小
ggsave("bioma_umap.pdf",width = 5,height = 5,dpi = 600)

clusterCornerAxes(object = Late,reduction = 'tsne',clusterCol = 'celltype',cellLabel = T,cellLabelSize = 5,
                  noSplit = T) +  scale_color_manual(values = alpha(cbPalette,0.65)) +  NoLegend() +
  scale_fill_manual(values = alpha(cbPalette,0.65))
ggsave("bioma_tsne.pdf",width = 5,height = 5,dpi = 600)

clusterCornerAxes(object = Late,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("bioma_umap_indi.pdf",width = 5.5,height = 5,dpi = 600)


clusterCornerAxes(object = Late,reduction = 'tsne',clusterCol = 'orig.ident',pSize = 0.01,
                  noSplit = T) 
ggsave("bioma_tsne_indi.pdf",width = 5.5,height = 5,dpi = 600)

clusterCornerAxes(object = Late,reduction = 'umap',clusterCol = 'orig.ident',pSize = 0.01,groupFacet = 'orig.ident',noSplit = F) + NoLegend()
ggsave("bioma_umap_indi_s.pdf",width = 10,height = 4,dpi = 600)

cluster4Marker=c(
  # PCT         
  "ALDOB","LRP2", "CUBN","SLC34A1", "SLC22A6","HAO2",  "GPX3", 
  ##内皮
  "VWF", "PECAM1","EMCN", "FLT1", "PLVAP", "SYNPO", "PODXL", "FN1", "ITGA1",
  ##系膜
  "PDGFRB", "ITGA8","COL6A2" ,  "MYL9","ACTA2", "COL1A1", "RGS5", "NOTCH3", "TAGLN", "LMOD1", "THBS1",
  ###LOH
  "UMOD", "SLC12A1", "CLDN16", "KCNJ1", "SLC12A3","CLCNKA", "CLCNKB",
  
  #集合管主细胞-DCT
  "SCNN1G", "SCNN1B", "AQP2",  
  ##T
  "IL7R", "CD3E", "TRAC", "GZMK", "PTPRC","CD3D",
  ##单核巨噬
  "CD163", "C1QC", "LYZ", "APOE", "CSF1R", "ITGAX", "CD14", "TLR2", "TLR4", "IL1B",
  #NK-T细胞
  "GNLY", "NKG7", "GZMB", "KLRD1","NCAM1",
  #"B细胞"
  "CD79A", "CD79B", "MS4A1", "IGKC", "BANK1",
  #"肥大细胞"
   "TPSAB1","TPSAB2", "CPA3", "KIT", "HDC",
  #"浆细胞"
  "IGHG", "IGHA", "JCHAIN", "CD38",
  # 中性粒细胞"
  "FCGR3B", "S100A8", "S100A9", "S100A12", "FPR1", "FPR2", "TLR1",  "CSF3R")


#没有集合管细胞
cluster4Marker=c(
    ##T最优组合：CD3E + TRAC
  "CD3E", "TRAC", "NKG7",
  
  # PCT         
  "ALDOB","LRP2", "CUBN",

  ##内皮
 "PECAM1","EMCN", 
 ##单核CD14,CD16/FCGR3A 巨噬CD163,C1QC
  "CD14", "CD163", 
  ##系膜
  "PDGFRB", "ACTA2",  "RGS5",
  ###LOH-DCT
  "UMOD", "SLC12A1","WNK1","PVALB",
 
  #"B细胞"
 "CD79B", "MS4A1",
  #"肥大细胞"最优组合：TPSAB1 + CPA3
  "TPSAB1", "CPA3",
 # 中性粒细胞"最优组合（推荐）：S100A8 + S100A9
 "S100A8", "S100A9",
 
 #"浆细胞"
    "CD38","SDC1","JCHAIN",
 #GZMB
  "GZMB","PTPRC","ALDH2","FN1"
 
)
 
 
# 使用更美观的渐变色方案（深蓝-浅蓝-白-橙-红）
pdf("包括ALDH2_FN1正确渐变色气泡图.pdf", width=12, height=6)

DotPlot(
  object = Late,
  features = cluster4Marker,
  cols = c("#ffffbf", "#d73027"),  # 浅黄→深红渐变
  col.min = 0,                     # 匹配图例最小值
  col.max = 10,                    # 匹配图例最大值
  dot.scale = 10                    # 适当调大气泡
) + 
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=12),
    panel.grid.major = element_line(color="grey90"),)  # 添加浅灰色网格线

dev.off()


#ͬ同上绘制marker的小提琴图

VlnPlot(Late,features = c("ALDH2", "FN1"),pt.size = 0,group.by = "celltype")&NoLegend()&labs(x = '') #&geom_hline(yintercept = 2500)
ggsave(filename = "ALDH2_FN1_violin.pdf",width=10,height=6)  

VlnPlot(Late, features = c("ALDH2", "FN1"), slot = "counts", log = TRUE)
ggsave(filename = "ALDH2_FN1_violin_Violin log(count).pdf",width=10,height=6) #取log值


# 寻找聚类4和聚类7的差异基因
clusterL_NE.markers <- FindMarkers(Late, ident.1 = "Luminal", ident.2 = "Luminal/NE", min.pct = 0.25,logfc.threshold = 0.5)
head(clusterL_NE.markers, n = 5)
write.table(clusterL_NE.markers,file="clusterL_NE.markers.xls",sep="\t",row.names=T,quote=F)
# 寻找聚类1和所以聚类的差异基因
cluster3.markers <- FindMarkers(Late, ident.1 = "EX CD8+T", min.pct = 0.25,logfc.threshold = 0.5)
head(cluster3.markers, n = 5)
write.table(cluster3.markers,file="EX CD8+T_markers.xls",sep="\t",row.names=T,quote=F)


#NEPC_CRPC比较
marker <- FindMarkers(Late, ident.1 = c(12,21), ident.2= c(7,10),min.pct = 0.25,only.pos = T,logfc.threshold = 0.5)
head(marker, n = 5)
marker<- subset(marker,p_val_adj<0.05)
marker=rbind(id=colnames(marker),marker)
write.table(marker,file="NEPC_CRPC_Luminal_0.25.markers.xls",sep="\t",quote=F,col.names=F)
#write.table(marker,file="NEPC_CRPC_Luminal_0.25.markers.xls",sep="\t",row.names=T,quote=F)
###min.pct表示基因在多少细胞中表达的阈值，only.pos = TRUE表示只求高表达的基因，

#NEPC_NS比较
marker <- FindMarkers(Late, ident.1 = c(12,21), ident.2= 3,min.pct = 0.25,only.pos = T,logfc.threshold = 0.5)
head(marker, n = 5)
marker<- subset(marker,p_val_adj<0.05)
marker=rbind(id=colnames(marker),marker)
write.table(marker,file="NEPC_NS_Luminal_0.25.markers.xls",sep="\t",quote=F,col.names=F)
#write.table(marker,file="NEPC_CRPC_Luminal_0.25.markers.xls",sep="\t",row.names=T,quote=F)
###min.pct表示基因在多少细胞中表达的阈值，only.pos = TRUE表示只求高表达的基因，


###组内差异分析##############
setwd("G:\\RCC\\GSE181061_ATAC_TCR_RNA\\scRNAseq_Tcell\\3.annotation")
library(dplyr)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
Late.markers <- FindAllMarkers(Late, min.pct = 0.25,logfc.threshold = 1) 
#因基因太多logfc.threshold = 1,基因少改为0.5#only.pos = TRUE,
Late.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)   ##n=2与n=10运行结果没多大区别?##
write.table(Late.markers,file="all_1.xls",sep="\t",row.names=F,quote=F)
Late = Late[, Idents(Late) %in% c( "EX CD8+T" , "C1-RM CD8+T","C2-RM CD8+T","Cytotoxic CD8+T")]#

###此代码出的结果过有部分PadjValue值不满足小于等于0.05,需要自己手工筛选#
Late.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#输出tsneTop10差异基因热图
DoHeatmap(Late, features = top10$gene) + NoLegend()
ggsave(filename = "02Top5gene_heatmap.pdf",width=15,height=15)
#save(Late, file = "pbmc_for_markers.rda")



#结合解螺旋的代码进行14个cluster注释#可成功运行
setwd("G:\\RCC\\GSE181061_ATAC_TCR_RNA\\scRNAseq_Tcell\\3.annotation") 
test <- GetAssayData(Late, slot="data") 
clusters<-Late@meta.data$seurat_clusters   
library(celldex)
library(SingleR)
HumanCell <- celldex::HumanPrimaryCellAtlasData()  #加载人类的
singler <- SingleR(test, ref = HumanCell,
                   labels = HumanCell$label.main,clusters = clusters)  
#test为需要注释的矩阵,ref = HumanCell参考矩阵为人类细胞文件#label.main主标签#  clusters为之前得到的聚类结果14个              
write.table(singler,file="03.17个cluster聚类注释2.txt",quote=F,sep="\t",col.names=F)
#包含微调前（first.labels）、微调后（labels）以及修剪后（pruned.labels),选择labels即可###

#对4504个细胞进行细胞注释
singler2 <- SingleR(test, ref = HumanCell, labels = HumanCell$label.main) #对每个细胞进行细胞注释
write.table(singler2,file="07.各个细胞类型注释-解.txt",quote=F,sep="\t",col.names=F)  
#需比对生信自学网产生的数据类型,看是否一致#

# 比较Seurat聚类和SignleR聚类的结果
table(singler2$labels,clusters)
write.table(table(singler2$labels,clusters),file="07.SingleR-Seurat聚类和SignleR聚类的结果对比.txt",quote=F,sep="\t",col.names=T)


#SingleR-绘制细胞类型结果热图Heatmap
pdf(file="07.SingleR-绘制细胞类型结果Heatmap.pdf",width=10,height=8)
plotScoreHeatmap(singler)
dev.off()


#同上,选择性绘制不同cluster中某些基因的散点图#和上面没啥区别?
pdf(file="06.PDCD1标记基因.pdf",width=7.5,height=6)
cluster4Marker=c("TIGIT", "CTLA4", "CD200R1", "LAG3", "TOX2", "HAVCR2", "PDCD1") 
cluster4Marker=c("LAG3", "HAVCR2", "PDCD1") 
FeaturePlot(Late, features =  cluster4Marker)#cols = c("green", "red")
dev.off()

#组间差异分析########
logFCfilter=0.5
adjPvalFilter=0.05

Late.markers=FindMarkers(Late, ident.1 = "Responder", ident.2 = "Non.Responder", group.by = 'patient')
Late.markers=cbind(Gene=row.names(Late.markers), Late.markers)
write.table(Late.markers,file="0.5allGene.txt",sep="\t",row.names=F,quote=F)
###vg_log2FC：记录两组之间平均表达的倍数变化
#。正值表示该基因在第一组中表达更高。


sig.markers=Late.markers[(abs(as.numeric(as.vector(Late.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(Late.markers$p_val_adj))<adjPvalFilter),]
sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
write.table(sig.markers,file="0.5diffGene.txt",sep="\t",row.names=F,quote=F)

