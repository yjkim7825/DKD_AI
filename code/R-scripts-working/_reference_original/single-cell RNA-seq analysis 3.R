
setwd("G:\\187geneMR\\24.GSE209781\\6.monocle")                 #设置工作目录

Late <- readRDS("G:\\187geneMR\\24.GSE209781\\6.monocle\\Late_with_geneset_scores.rds")
#导入注释好的seurat对象（已注释）


responder <- subset(Late, idents = c("PCT", "dPCT"))
table(idents(responder))  # 检查细胞数量

library(monocle)
library(dplyr)
library(Matrix)
library(Seurat)
library(ggplot2)  

library(libcoin)
library(partykit)


##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
#expr_matrix <- as(as.matrix(responder@assays[["RNA"]]@counts), 'sparseMatrix')
expr_matrix <- as.sparse(responder@assays[["RNA"]]@layers[["counts"]])
##提取表型信息到p_data(phenotype_data)里面 
p_data <- responder@meta.data 
p_data$celltype <- responder@active.ident  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(responder),row.names = row.names(responder))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
# 在运行任何Monocle函数前执行

cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit=.5,
                      expressionFamily = negbinomial.size())
save(cds, file = "before_estimateSizeFactors.rda")

#伪时间分析流程分支1和分支2
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)#,cores=4,relative_expr=T
cds=detectGenes(cds,min_expr = 0.1) #计算每个基因在多少细胞中表达
print(head(fData(cds)))  

save(cds, file = "select_gene_for_monocle.rda")

load("select_gene_for_monocle.rda")
#下一步选择要构建轨迹所用的基因,对轨迹影响很大
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10)) 
#过滤掉在小于10个细胞中表达的基因
#也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(responder) 

diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~celltype") 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)
##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
##差异基因的结果文件保存
write.table(deg,file="gene_for_trajectory.xls",col.names=T,row.names=F,sep="\t",quote=F)

## 轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
save(cds, file = "cds-for_plot.rda")
#load("G:/RCC/sc_immune_therapy/6.monocle/cds-for_plot.rda")


pdf("train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()

cds=reduceDimension(cds,reduction_method = "DDRTree",max_components = 2) 

cds=orderCells(cds)

save(cds, file = "cds-after_orderCells.rda")

# 加载ggplot2包（如果尚未加载）
library(ggplot2)

# 1. 绘制伪时间轨迹图 -----------------------------------------------------------------
pdf("train.monocle.pseudotime.pdf", width = 5, height = 5)  # 创建PDF文件，设置宽高为7英寸
plot_cell_trajectory(
  cds, 
  color_by = "Pseudotime",  # 按伪时间值着色
  size = 1,                 # 点的大小
  show_backbone = TRUE      # 显示主干轨迹线
) + 
  theme(
    text = element_text(size = 12),        # 全局文本大小
    axis.title = element_text(size = 14),   # 坐标轴标题大小
    axis.text = element_text(size = 12),    # 坐标轴刻度标签大小
    legend.title = element_text(size = 14), # 图例标题大小
    legend.text = element_text(size = 12)   # 图例项目文本大小
  )
dev.off()  # 关闭图形设备


# 2. 绘制细胞类型轨迹图 --------------------------------------------------------------
pdf("train.monocle.celltype.pdf", width = 4, height = 4)
plot_cell_trajectory(
  cds,
  color_by = "celltype",  # 按细胞类型着色
  size = 1,
  show_backbone = TRUE
) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()


# 3. 绘制状态轨迹图 ------------------------------------------------------------------
pdf("train.monocle.state.pdf", width = 5, height = 5)
plot_cell_trajectory(
  cds, 
  color_by = "State",  # 按Monocle状态着色
  size = 1,
  show_backbone = TRUE
) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()


# 4. 绘制Seurat聚类轨迹图 ------------------------------------------------------------
pdf("seurat.clusters.pdf", width = 5, height = 5)
plot_cell_trajectory(
  cds, 
  color_by = "seurat_clusters"  # 按Seurat聚类结果着色
) +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()


# 5. 分面绘制状态轨迹图（横向排列）--------------------------------------------------
pdf("train.monocle.state.faceted.pdf", width = 10, height = 5)
plot_cell_trajectory(
  cds, 
  color_by = "State"
) + 
  facet_wrap("~State", nrow = 1) +  # 按State分面，单行排列
  theme(
    text = element_text(size = 14),         # 增大基础字号
    axis.title = element_text(size = 16),    # 增大坐标轴标题
    axis.text = element_text(size = 14),     # 增大刻度标签
    legend.title = element_text(size = 16),  # 增大图例标题
    legend.text = element_text(size = 14),   # 增大图例文本
    strip.text = element_text(size = 14)     # 分面标签大小
  )
dev.off()

# 加载包
library(monocle)
library(ggplot2)
library(tidyr)

# 目标基因"FN1", "ALDH2"拟时序分析####
target_genes <- c("FN1", "ALDH2")  # 确保基因名与数据匹配

# 检查基因是否存在
missing_genes <- setdiff(target_genes, rownames(cds))
if (length(missing_genes) > 0) {
  stop(paste("以下基因不在数据中:", paste(missing_genes, collapse = ", ")))
}

# 1. 拟时序热图 (PDF输出)
pdf("FN1_ALDH2_pseudotime_heatmap.pdf", width = 8, height = 5)
plot_genes_in_pseudotime(cds[target_genes, ], 
                         color_by = "celltype",
                         ncol = 2) +
  ggtitle("FN1 and ALDH2 Expression along Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

library(tidyr)
# 2. 拟时序趋势曲线 (PDF输出)
pdf("FN1_ALDH2_trend_curves.pdf", width = 3, height = 3)
exprs <- exprs(cds)[target_genes, ]
pdata <- pData(cds)
df <- data.frame(
  Pseudotime = pdata$Pseudotime,
  State = pdata$State,
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df, aes(x = Pseudotime, y = Expression, color = Gene)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("FN1" = "red", "ALDH2" = "#0072B2")) +
  labs(x = "Pseudotime", y = "Expression", 
       title = "Expression Trends") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


#####3. 拟时序趋势曲线 (PDF输出)
pdf("FN1_ALDH2_trend_curves_by_celltype.pdf", width = 6, height = 4)
exprs <- exprs(cds)[target_genes, ]
pdata <- pData(cds)
df <- data.frame(
  Pseudotime = pdata$Pseudotime,
  CellType = pdata$celltype,  # 确保列名与您的metadata一致
  State = pdata$State,
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df, aes(x = Pseudotime, y = Expression, color = Gene)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("FN1" = "red", "ALDH2" = "#0072B2")) +
  labs(x = "Pseudotime", y = "Expression", 
       title = "FN1 and ALDH2 Expression Trends by Cell Type") +
  facet_wrap(~ CellType, ncol = 3) +  # 按细胞类型分面
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())
dev.off()




# 5. 小提琴图 (PDF输出)
pdf("FN1_ALDH2_violin_by_celltype.pdf", width = 8, height = 5)
df_violin <- data.frame(
  celltype = as.factor(pdata$celltype),
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df_violin, aes(x = celltype, y = Expression, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("FN1" = "#d9352a", "ALDH2" = "#0072B2")) +
  labs(x = "celltype", y = "Expression", 
       title = "FN1 and ALDH2 Expression by celltype") +
  theme_minimal(base_size = 12) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# 5. 小提琴图 (PDF输出)
# 5. 小提琴图 (PDF输出)
pdf("FN1_ALDH2_violin_by_state.pdf", width = 8, height = 5)
df_violin <- data.frame(
  State = as.factor(pdata$State),
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
) %>% 
  pivot_longer(cols = c("FN1", "ALDH2"), 
               names_to = "Gene", 
               values_to = "Expression")

ggplot(df_violin, aes(x = State, y = Expression, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("FN1" = "red", "ALDH2" = "#0072B2")) +
  labs(x = "State", y = "Expression", 
       title = "FN1 and ALDH2 Expression by State") +
  theme_minimal(base_size = 12) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


## 目标基因ALDOB_LRP2_SLC34A1拟时序小提琴图####

target_genes <- c("ALDOB", "LRP2", "SLC34A1")
stopifnot(all(target_genes %in% rownames(exprs(cds))))

# 创建小提琴图数据
df_violin <- data.frame(
  celltype = as.factor(pData(cds)$celltype),
  ALDOB = exprs(cds)["ALDOB", ],
  LRP2 = exprs(cds)["LRP2", ],
  SLC34A1 = exprs(cds)["SLC34A1", ]
) %>% 
  pivot_longer(cols = all_of(target_genes),
               names_to = "Gene", 
               values_to = "Expression")

# 绘制小提琴图
pdf("ALDOB_LRP2_SLC34A1_violin_by_celltype.pdf", width = 8, height = 6)
ggplot(df_violin, aes(x = celltype, y = Expression, fill = Gene)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
  scale_fill_manual(values = c("ALDOB" = "#4979b6", "LRP2" = "#d98d2a", "SLC34A1" = "#d9352a")) +
  labs(x = "Cell Type", y = "Expression") +
  theme_minimal(base_size = 14) +
  facet_wrap(~ Gene, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
dev.off()


######绘制基因表达与基因集评分的关系(没啥用)########
# 确保已加载必要的包
library(monocle)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 提取关键数据
pdata <- pData(cds)
exprs <- exprs(cds)

# 检查FN1和ALDH2是否在表达矩阵中
if (!all(c("FN1", "ALDH2") %in% rownames(exprs))) {
  stop("FN1或ALDH2不在表达矩阵中，请检查基因名大小写或是否过滤。")
}

# 合并表达数据和基因集评分
plot_data <- data.frame(
  pdata[, c("Pseudotime", "State", "celltype",
            "EPITHELIAL_MESENCHYMAL_TRANSITION_1", "CELL_AGING_1", 
            "AUTOPHAGY_1", "OXIDATIVE_STRESS_1", 
            "APOPTOSIS_1", "INFLAMMATORY_RESPONSE_1")],
  FN1 = exprs["FN1", ],
  ALDH2 = exprs["ALDH2", ]
)

# 查看数据结构
head(plot_data)


#散点图 + 趋势线（分面显示）
pdf("FN1_ALDH2_vs_Geneset_Scores.pdf", width = 15, height = 20)

# FN1与6个基因集评分的关系
fn1_plots <- lapply(colnames(plot_data)[4:9], function(geneset) {
  ggplot(plot_data, aes_string(x = "FN1", y = geneset)) +
    geom_point(aes(color = celltype), alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(x = "FN1 Expression", y = geneset, 
         title = paste0("FN1 vs ", geneset)) +
    theme_bw() +
    theme(legend.position = "bottom")
})

# ALDH2与6个基因集评分的关系
aldh2_plots <- lapply(colnames(plot_data)[4:9], function(geneset) {
  ggplot(plot_data, aes_string(x = "ALDH2", y = geneset)) +
    geom_point(aes(color = celltype), alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(x = "ALDH2 Expression", y = geneset, 
         title = paste0("ALDH2 vs ", geneset)) +
    theme_bw() +
    theme(legend.position = "bottom")
})

# 合并所有图
(wrap_plots(fn1_plots, ncol = 3) / wrap_plots(aldh2_plots, ncol = 3)) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

dev.off()

# 计算Pearson相关系数
cor_results <- sapply(colnames(plot_data)[4:9], function(geneset) {
  c(
    FN1_cor = cor(plot_data$FN1, plot_data[[geneset]], method = "pearson"),
    ALDH2_cor = cor(plot_data$ALDH2, plot_data[[geneset]], method = "pearson")
  )
})

# 输出相关系数表格
write.csv(cor_results, "Gene_Geneset_Correlation_Results.csv")

# 热图可视化
pdf("Correlation_Heatmap.pdf", width = 6, height = 4)
pheatmap::pheatmap(cor_results,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   display_numbers = TRUE,
                   number_format = "%.2f",
                   main = "Pearson Correlation: FN1/ALDH2 vs Gene Sets")
dev.off()

pdf("Geneset_Scores_by_Celltype_and_Gene_Expression.pdf", width = 14, height = 10)

# 将FN1/ALDH2表达分为高/低两组（按中位数）
plot_data <- plot_data %>%
  mutate(
    FN1_group = ifelse(FN1 > median(FN1), "High", "Low"),
    ALDH2_group = ifelse(ALDH2 > median(ALDH2), "High", "Low")
  )

# 对每个基因集绘制分组箱线图
for (geneset in colnames(plot_data)[4:9]) {
  p1 <- ggplot(plot_data, aes_string(x = "celltype", y = geneset, fill = "FN1_group")) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x = "Cell Type", y = geneset, 
         title = paste0(geneset, " Score by Cell Type and FN1 Expression")) +
    theme_bw() +
    scale_fill_manual(values = c("High" = "#E64B35", "Low" = "#7AA6DC"))
  
  p2 <- ggplot(plot_data, aes_string(x = "celltype", y = geneset, fill = "ALDH2_group")) +
    geom_boxplot(outlier.size = 0.5) +
    labs(x = "Cell Type", y = geneset, 
         title = paste0(geneset, " Score by Cell Type and ALDH2 Expression")) +
    theme_bw() +
    scale_fill_manual(values = c("High" = "#4DBBD5", "Low" = "#FFC77D"))
  
  print(p1 + p2 + plot_layout(ncol = 2))
}
dev.off()

######可视化基因集分数沿伪时间的变化,以celltype#######

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 提取关键数据
library(stringr)
plot_data <- pData(cds) %>%
  select(Pseudotime, celltype, 
         EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs_1:INFLAMMATORY_RESPONSE.v2024.1.Hs_1) %>%
  rename_with(~ str_remove(., "\\.v2024\\.1\\.Hs"), 
              EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs_1:INFLAMMATORY_RESPONSE.v2024.1.Hs_1) %>%
  pivot_longer(cols = -c(Pseudotime, celltype),
               names_to = "Geneset",
               values_to = "Score")
# 移除Geneset名称中的"_1"

# 查看数据结构
head(plot_data)
plot_data <- plot_data %>%
  mutate(Geneset = stringr::str_remove(Geneset, "_1$"))


pdf("Geneset_Scores_along_Pseudotime_by_Celltype.pdf", width = 10, height = 8)

ggplot(plot_data, aes(x = Pseudotime, y = Score, color = celltype)) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  facet_wrap(~ Geneset, scales = "free_y", ncol = 3) +
  labs(x = "Pseudotime", y = "Geneset Score", 
       title = "Geneset Scores along Pseudotime by Cell Type") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("PCT" = "#1F77B4", 
                                "dPCT" = "red"))

dev.off()

pdf("Geneset_Scores_vs_Pseudotime_Scatter.pdf", width = 12, height = 8)

ggplot(plot_data, aes(x = Pseudotime, y = Score, color = celltype)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  facet_grid(Geneset ~ celltype, scales = "free_y") +
  labs(x = "Pseudotime", y = "Score", 
       title = "Geneset Scores Distribution by Cell Type") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8)) +
  scale_color_manual(values = c("PCT" = "#1F77B4", 
                                "dPCT1" = "#FF7F0E", 
                                "dPCT2" = "red"))

dev.off()

# 将拟时序分为5段
plot_data <- plot_data %>%
  mutate(Pseudotime_bin = cut(Pseudotime, breaks = 5, labels = FALSE))

pdf("Geneset_Scores_by_Pseudotime_Bins.pdf", width = 14, height = 10)

ggplot(plot_data, aes(x = as.factor(Pseudotime_bin), y = Score, fill = celltype)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Geneset, scales = "free_y", ncol = 3) +
  labs(x = "Pseudotime Bins", y = "Score", 
       title = "Geneset Scores by Pseudotime Segments") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("PCT" = "#1F77B4", 
                               "dPCT1" = "#FF7F0E", 
                               "dPCT2" = "red"))
dev.off()


#######以下是旧代码#########

#寻找拟时间相关的基因
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
stares_de <- differentialGeneTest(cds[ordergene,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
stares_de <- stares_de[,c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
stares_de <- stares_de[order(stares_de$qval), ]
write.csv(stares_de, "states_diff_gene.csv", row.names = F)  

Time_genes <- stares_de %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)



#寻找以依赖于分支的方式调控的基因
BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 2) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res_branch_point_one.csv", row.names = F)
#branch_point = 1产生的是node1处的差异基因
pdf("genes_branched_heatmap",width = 10,height = 7)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了
dev.off()

BEAM_res <- BEAM(cds[ordergene,], branch_point = 2, cores = 2) 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
#BEAM_res <- BEAM(cds, branch_point = 1, cores = 2) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res, "BEAM_res_branch_point_two.csv", row.names = F)

library(ggpubr)
df <- pData(cds) 
## pData(cds)取出的是cds对象中cds@phenoData@data的内容
#View(df)

pdf("geom_density.pdf",width = 10,height = 7)
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()
dev.off()

#手动设置颜色
ClusterName_color_panel <- c(
  "Epithelial Luminal" = "#DC143C", "Epithelial Other" = "#0000FF", "Endothelial" = "#20B2AA",
  "Fibroblast" = "#FFA500", "Monocyte" = "#9370DB", "B cell" = "#98FB98",
  "Epithelial Basal" = "#F08080", "Epithelial NE" = "#0000FF")  #, "Platelet" = "#20B2AA"
  pdf("geom_density2.pdf",width = 10,height = 7)
ggplot(df, aes(Pseudotime, colour = celltype, fill=celltype)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ scale_fill_manual(name = "", values = ClusterName_color_panel)+scale_color_manual(name = "", values = ClusterName_color_panel)
dev.off()
