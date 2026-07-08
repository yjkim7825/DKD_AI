setwd("G:\\187geneMR\\24.GSE209781\\4.Cell ratio\\6 sample")    
load("G:\\187geneMR\\24.GSE209781\\4.Cell ratio\\6 sample\\Late_annotated.RDA")

unique(Late$orig.ident)
table(Late$celltype)
table(Late$patient)

###########堆积柱状图#############
library(ggplot2)
library(ggpubr)
cellnum <- table(Late$celltype,Late$orig.ident)
cell.prop<-as.data.frame(prop.table(cellnum))
colnames(cell.prop)<-c("Celltype","Group","Proportion")
View(cell.prop)
cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
               "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272',"#fc8542")

pdf(file="6样本堆积柱状图.pdf",width=6,height=6)   
ggplot(cell.prop,aes(Group,Proportion,fill=Celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=cbPalette)+#自定义fill的颜色
  ggtitle("cell proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))
dev.off()


# 直接使用patient列分组计算（关键修改）

library(ggplot2)
library(ggpubr)

cellnum <- table(Late$celltype, Late$patient)  # 将orig.ident改为patient
cell.prop <- as.data.frame(prop.table(cellnum))
colnames(cell.prop) <- c("Celltype", "Group", "Proportion")

# 保持您的颜色和绘图代码完全不变
cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
               "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272',"#fc8542")

pdf(file="两样本堆积柱状图.pdf", width=4.5, height=6)  # 调整宽度适应两列   
ggplot(cell.prop, aes(Group, Proportion, fill=Celltype))+
  geom_bar(stat="identity", position="fill")+
  scale_fill_manual(values=cbPalette)+
  ggtitle("cell proportion")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'))+
  guides(fill=guide_legend(title=NULL))
dev.off()

##########展示特定细胞堆积柱状图####
library(ggplot2)
library(ggpubr)

# 计算细胞比例
cellnum <- table(Late$celltype, Late$patient)  # 统计各细胞类型在各组的数量
cell.prop <- as.data.frame(prop.table(cellnum))  # 计算比例
colnames(cell.prop) <- c("Celltype", "Group", "Proportion")

# 仅保留 PCT, dPCT1, dPCT2, dPCT3
target_cells <- c("PCT", "dPCT")
cell.prop_filtered <- subset(cell.prop, Celltype %in% target_cells)

# 自定义颜色（确保颜色数量足够）
cbPalette <- c("#009E73","#E69F00" )  # 4种颜色对应4种细胞类型

# 绘图
pdf(file="PCT_proportion堆积柱状图.pdf", width=4.5, height=6)  
ggplot(cell.prop_filtered, aes(Group, Proportion, fill=Celltype)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=cbPalette) +
  ggtitle("Proximal Tubule Cell Proportions") +  # 修改标题
  theme_bw() +
  theme(
    axis.ticks.length=unit(0.1,'cm'),
    plot.title = element_text(hjust=0.5, size=12),  # 标题居中
    legend.position="right"  # 图例在右侧
  ) +
  guides(fill=guide_legend(title="Cell Type"))  # 图例标题
dev.off()

library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)

# ============= 用户只需修改这里 =============
gene <- "FN1"  # 只需修改这个基因名即可

# ===========================================

# 颜色定义（全局统一）
group_colors <- c("Control" = "#3B4992", "Late_DKD" = "#EE0000")
celltype_palette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
                     "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75',
                     '#fac858','#ee6666','#73c0de','#3ba272',"#fc8542")

###### 1. marker在Normal/DKD分组散点图 ########
p_scatter <- FeaturePlot(
  Late,
  features = gene,
  reduction = "tsne",
  split.by = "patient",
  pt.size = 1,
  order = TRUE
) +
  plot_annotation(title = paste0(gene, " Expression in Normal vs DKD"))

ggsave(paste0(gene, "_scatter_Normal_vs_DKD.pdf"), p_scatter, width = 10, height = 5)

###### 2. ALDH2在正常和对照组小提琴图（所有细胞）#####
expr_data <- FetchData(Late, vars = c(gene, "patient"))
colnames(expr_data) <- c("expression", "group")

# 统计检验
wilcox_test <- wilcox.test(expression ~ group, data = expr_data)
p_value <- format.pval(wilcox_test$p.value, digits = 3)
control_mean <- mean(expr_data$expression[expr_data$group == "Control"])
dkd_mean <- mean(expr_data$expression[expr_data$group == "Late_DKD"])
log2fc <- log2(dkd_mean / control_mean)

# 绘图
p_violin_all <- ggplot(expr_data, aes(x = group, y = expression, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = paste0(gene, " Expression in Single Cells"),
    subtitle = paste0("Control (n=", sum(expr_data$group == "Control"), 
                     ") vs Late_DKD (n=", sum(expr_data$group == "Late_DKD"), ")"),
    y = "Expression Level",
    x = "",
    caption = paste0("p = ", p_value, " | Log2FC = ", round(log2fc, 2))
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

ggsave(paste0(gene, "_single_cell_expression.pdf"), plot = p_violin_all, width = 6, height = 5)

###### 3. ALDH2在所有细胞类型小提琴图(不分组) ####

#####方法1首选

VlnPlot(Late,features = c("ALDH2", "FN1"),pt.size = 0,group.by = "celltype")&NoLegend()&labs(x = '') #&geom_hline(yintercept = 2500)
ggsave(filename = "ALDH2_FN1_violin(所有细胞不分组).pdf",width=10,height=6)  

# 绘制 ALDOB、CUBN 和 SLC34A1 在 PCT、dPCT中的表达小提琴图
VlnPlot(
  Late,
  features = c("ALDOB", "CUBN","SLC34A1"),
  pt.size = 0,
  group.by = "celltype",
  idents = c("PCT", "dPCT"),
  cols = c("#0072B2", "#FF7F0E","red")  # 自定义颜色（顺序对应 PCT, dPCT1, dPCT2）
) &
  NoLegend() &
  labs(x = '') &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 保存图片
ggsave(filename = "ALDOB_SLC34A1_CUBN_violin(PCT_dPCTs).pdf", width=4, height=3)

# 绘制 "ALDH2", "FN1"在 PCT、dPCT中的表达小提琴图
VlnPlot(
  Late,
  features = c("ALDH2", "FN1"),
  pt.size = 0,
  group.by = "celltype",
  idents = c("PCT", "dPCT"),
  cols = c("#0072B2", "#FF7F0E","red")  # 自定义颜色（顺序对应 PCT, dPCT1, dPCT2）
) &
  NoLegend() &
  labs(x = '') &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 保存图片
ggsave(filename = "ALDH2_FN1_violin(PCT_dPCTs).pdf", width=4, height=3)


####方法2,次选
p_celltypes <- VlnPlot(
  object = Late,
  features = gene,
  group.by = "celltype",
  pt.size = 0,
  cols = celltype_palette
) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  labs(
    title = paste0(gene, " Expression Across Cell Types"),
    y = "Expression Level",
    x = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
  )

ggsave(paste0(gene, "_all_celltypes_ungrouped.pdf"), plot = p_celltypes, width = 8, height = 5)

###### 4. ALDH2在所有细胞类型小提琴图(分组) ####
# 计算每个细胞类型的p值
p_values <- sapply(unique(Idents(Late)), function(celltype) {
  subset_cells <- subset(Late, idents = celltype)
  if(sum(subset_cells$patient == "Control") >= 3 && 
     sum(subset_cells$patient == "Late_DKD") >= 3) {
    wilcox.test(FetchData(subset_cells, vars = gene)[subset_cells$patient == "Control", 1],
                FetchData(subset_cells, vars = gene)[subset_cells$patient == "Late_DKD", 1])$p.value
  } else NA
})

# 绘图
p_celltypes_grouped <- VlnPlot(
  object = Late,
  features = gene,
  split.by = "patient",
  group.by = "celltype",
  pt.size = 0,
  cols = group_colors
) +
  labs(
    title = paste0(gene, " Expression by Cell Type"),
    y = "Expression Level",
    x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right"
  )

# 添加p值标注
y_max <- max(FetchData(Late, vars = gene)[, 1]) * 1.05
for (i in seq_along(unique(Idents(Late)))) {
  if (!is.na(p_values[i])) {
    p_celltypes_grouped <- p_celltypes_grouped + 
      annotate("text", x = i, y = y_max, 
               label = paste0("p=", format.pval(p_values[i], digits = 3)), 
               size = 3, vjust = 0)
  }
}

ggsave(paste0(gene, "_violin_by_celltype_grouped.pdf"), plot = p_celltypes_grouped, width = 10, height = 6)


###### 4. marker在两种细胞类型小提琴图(不分组) ####
library(Seurat)
library(ggplot2)
library(ggpubr)

# 设置参数
target_celltypes <- c("PCT", "dPCT")  # 目标细胞类型
group_colors <- c("#0072B2", "#FF7F0E")  # 颜色
# 提取目标细胞
target_cells <- subset(Late, idents = target_celltypes)

# 准备数据（确保列名正确）
plot_data <- FetchData(target_cells, vars = c(gene, "celltype"))
colnames(plot_data) <- c("expression", "celltype")

# 计算p值（Wilcoxon检验）
p_value <- wilcox.test(expression ~ celltype, data = plot_data)$p.value
p_label <- ifelse(p_value < 0.001, 
                  paste0("p = ", format(p_value, scientific = TRUE, digits = 2)),
                  paste0("p = ", round(p_value, 3)))  # 强制显示为 p = 数值

# 计算标注位置
y_max <- max(plot_data$expression) * 0.9  # 提高标注位置

# 绘制小提琴图（关键修改点）
p <- VlnPlot(
  object = target_cells,
  features = gene,
  group.by = "celltype",
  pt.size = 0,
  cols = group_colors
) +
  theme_classic() +
  labs(
    title = paste0(gene, " Expression in PCT vs. dPCT"),
    y = "Expression Level",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # 彻底移除图例
  ) +
  annotate("text", 
           x = 1.5, 
           y = y_max, 
           label = p_label, 
           size = 5)

# 保存PDF
ggsave(
  filename = paste0(gene, "_PCT_vs_dPCT_expression.pdf"),
  plot = p,
  width = 4,
  height = 6
)

print(p)
###### 5. ALDH2在特定细胞类型小提琴图(分组) ####
target_celltype <- "MES"  # 设置目标细胞类型（用于特定细胞分析）
target_cells <- subset(Late, idents = target_celltype)
group_colors <- c("#0072B2", "#FF7F0E")  # 颜色
# 统计检验
control_expr <- FetchData(target_cells, vars = gene)[target_cells$patient == "Control", 1]
dkd_expr <- FetchData(target_cells, vars = gene)[target_cells$patient == "Late_DKD", 1]
wilcox_test <- wilcox.test(control_expr, dkd_expr)

# 科学计数法格式化p值
format_p_value <- function(p) {
  if (p < 0.001) {
    return(format(p, scientific = TRUE, digits = 2))
  } else {
    return(round(p, 3))
  }
}

p_value <- format_p_value(wilcox_test$p.value)
y_max <- max(FetchData(target_cells, vars = gene)[, 1]) * 0.9

# 绘图
p_specific <- VlnPlot(
  object = target_cells,
  features = gene,
  group.by = "patient",
  pt.size = 0,
  cols = group_colors
) +
  labs(
    title = paste0(gene, " in ", target_celltype),
    y = "Expression Level",
    x = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  annotate("text", x = 1.5, y = y_max, 
           label = paste0("p = ", p_value), size = 5)

ggsave(paste0(gene, "_", target_celltype, "_expression_violin.pdf"), 
       plot = p_specific, width = 2.5, height = 3)
print(p_specific)  # 注意这里改为p_specific
###############基因集富集评分###########
cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
               "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75','#fac858','#ee6666','#73c0de','#3ba272',"#fc8542")
library(Seurat)
library(msigdbr)  # 用于加载GMT基因集

# 1. 加载基因集（示例，需替换为你的GMT文件）
library(GSEABase)

# 确保文件路径正确（使用绝对路径）
gmt_file <- "EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt"
gene_sets <- getGmt(gmt_file)

# 提取第一个基因集的基因列表
geneset_name <- names(gene_sets)[1]
geneset_list <- geneIds(gene_sets)[[1]]

# 强行计算（假设geneset_list可能包含不匹配的基因）
# 不检查基因匹配情况，直接计算
EMT_score <- AddModuleScore(
  object = Late,  # 假设Late是你的Seurat对象
  features = list(EMT_geneset = geneset_list),
  name = "EMT_Score",
  ctrl = 5,  # 控制基因数量
  seed = 123  # 设置随机种子保证可重复性
)

# 将分数添加到metadata中方便后续使用
EMT_score$EMT_Score1 <- EMT_score$EMT_Score1  # AddModuleScore会自动在名字后加数字


# 3.1 箱线图展示不同分组间的EMT分数差异
# 3.1 优化箱线图（保留标题+标注TOP3）
# 计算各细胞类型的中位数评分
median_scores <- EMT_score@meta.data %>%
  group_by(celltype) %>%
  summarise(median_score = median(EMT_Score1)) %>%
  arrange(desc(median_score)) 

top3_celltypes <- median_scores$celltype[1:3]

p1 <- ggplot(EMT_score@meta.data, aes(x = reorder(celltype, -EMT_Score1, median), y = EMT_Score1)) +
  geom_boxplot(aes(fill = celltype), width = 0.6, outlier.shape = NA, lwd = 0.5) +
  scale_fill_manual(values = cbPalette) +
  theme_classic(base_size = 14) +
  labs(title = paste(geneset_name),
       x = "Cell Type", 
       y = "AddModuleScore") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  ) +
  # 标注TOP3
  geom_text(data = subset(median_scores, celltype %in% top3_celltypes),
            aes(x = celltype, y = max(EMT_score$EMT_Score1) * 1.05, 
                label = paste0("Top ", match(celltype, top3_celltypes))),
            size = 5, color = "red", fontface = "bold")
# 保存
ggsave("1_EMT_score_with_top3.pdf", p1, width = 9, height = 6, dpi = 300)


# 3.2 UMAP展示EMT分数分布
p2 <- FeaturePlot(EMT_score, 
                  features = "EMT_Score1",
                  reduction = "tsne",
                  pt.size = 0.8) +  # 适当增大点尺寸
  scale_colour_gradientn(
    colours = rev(rainbow(4)),  # 保持原始彩虹配色
    name = "EMT Score",  # 修改图例标题
    guide = guide_colorbar(
      barwidth = 0.8,
      barheight = 5,
      frame.colour = "black",  # 添加图例外框
      ticks.colour = "black"   # 刻度线颜色
    )
  ) +
  ggtitle(geneset_name) +  # 简化标题
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 16, 
      face = "bold",
      margin = margin(b = 10)  # 增加标题底部间距
    ),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(
      size = 12,
      face = "bold",
      hjust = 0.5  # 图例标题居中
    ),
    panel.background = element_rect(fill = "white")  # 纯白背景
  )

# 保存
ggsave("2_EMT_score_UMAP_simple.pdf", 
       plot = p2, 
       width = 6.5, 
       height = 5.5,
       dpi = 300)

##########批量基因集评分#########
############### 完整基因集评分分析流程 ###########
############### 完整基因集评分分析流程（修正版）###########
library(Seurat)
library(GSEABase)
library(ggplot2)
library(dplyr)  # 确保dplyr加载
library(writexl)
library(RColorBrewer)

# 1. 初始化设置 ------------------------------------------------------------
cbPalette <- c("#999999","#009E73","#56B4E9", "#E69F00", "#F0E442", 
               "#CC79A7","#D55E00","#0072B2",'#5470c6','#91cc75',
               '#fac858','#ee6666','#73c0de','#3ba272',"#fc8542")

# 2. 基因集文件定义 --------------------------------------------------------
geneset_files <- c(
  "EPITHELIAL_MESENCHYMAL_TRANSITION.gmt",
  "CELL_AGING.txt",
  "AUTOPHAGY.txt",
  "OXIDATIVE_STRESS.txt",
  "APOPTOSIS.gmt",
  "INFLAMMATORY_RESPONSE.gmt"
)

# 3. 基因集读取函数 -------------------------------------------------------
read_geneset <- function(file) {
  if (grepl("\\.gmt$", file)) {
    gmt <- getGmt(file)
    genes <- geneIds(gmt)[[1]]  # 取第一个基因集
  } else {
    genes <- readLines(file) %>% trimws()
    genes <- genes[genes != ""]
  }
  return(genes)
}

# 4. 基因匹配检查函数 -----------------------------------------------------
check_gene_match <- function(seurat_obj, geneset_list) {
  all_genes <- rownames(seurat_obj)
  
  lapply(geneset_list, function(genes) {
    total <- length(genes)
    found <- sum(genes %in% all_genes)
    percentage <- round(found/total * 100, 1)
    missing <- setdiff(genes, all_genes)
    
    list(
      total = total,
      found = found,
      percentage = percentage,
      missing = missing
    )
  })
}

# 5. 主分析流程 -----------------------------------------------------------

# 5.1 读取所有基因集
all_genesets <- lapply(geneset_files, function(file) {
  set_name <- tools::file_path_sans_ext(basename(file))
  genes <- read_geneset(file)
  list(name = set_name, genes = genes)
}) 
names(all_genesets) <- sapply(all_genesets, function(x) x$name)

# 5.2 检查基因匹配
match_results <- check_gene_match(Late, lapply(all_genesets, function(x) x$genes))

# 打印匹配结果
cat("\n=== 基因集匹配情况 ===\n")
match_summary <- data.frame(
  Geneset = names(match_results),
  Total = sapply(match_results, function(x) x$total),
  Found = sapply(match_results, function(x) x$found),
  Percentage = sapply(match_results, function(x) x$percentage)
)
print(match_summary)

# 5.3 计算模块评分
calculate_scores <- function(seurat_obj, geneset_list) {
  for(set in geneset_list) {
    clean_name <- gsub("\\s+", "_", set$name)
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(set$genes),
      name = paste0(clean_name, "_"),
      ctrl = min(100, length(set$genes)/2),
      seed = 123
    )
  }
  return(seurat_obj)
}

Late <- calculate_scores(Late, all_genesets)

# 5.4 可视化函数（修正select问题）
create_plots <- function(seurat_obj, score_name, title) {
  # 准备数据（使用明确的dplyr::select）
  plot_data <- seurat_obj@meta.data %>%
    dplyr::select(celltype, score = paste0(score_name, "_1"))
  
  # 按中位数排序
  med_order <- plot_data %>%
    group_by(celltype) %>%
    summarise(med = median(score)) %>%
    arrange(desc(med))
  
  # 小提琴图
  p_violin <- ggplot(plot_data, 
                     aes(x = factor(celltype, levels = med_order$celltype), 
                         y = score, 
                         fill = celltype)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = cbPalette) +
    labs(title = paste(title, "Score"), 
         subtitle = paste0("Genes: ", match_results[[title]]$found, "/", match_results[[title]]$total),
         x = "Cell Type (sorted by median)", 
         y = "Module Score") +
    theme_classic() +
    theme(
      # 主标题：加大到16pt，加粗
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      
      # 副标题：加大到12pt
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      
      # 坐标轴标题：加大到14pt
      axis.title = element_text(size = 14),
      
      # 坐标轴刻度：x轴加大到12pt并保持45度倾斜
      axis.text = element_text(size = 12),  # 影响x和y轴
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      
      # 图例：完全隐藏（保持不变）
      legend.position = "none"
    )
  
  # t-SNE图
  p_tsne <- FeaturePlot(seurat_obj,
                        features = paste0(score_name, "_1"),
                        reduction = "tsne",
                        pt.size = 0.8) +
    scale_color_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
    labs(title = paste(title, "Distribution"))
  
  list(violin = p_violin, tsne = p_tsne)
}

# 5.5 生成所有结果
results <- lapply(all_genesets, function(set) {
  set_name <- set$name
  clean_name <- gsub("\\s+", "_", set_name)
  message("Processing: ", set_name)
  
  plots <- create_plots(Late, clean_name, set_name)
  
  # 保存图片
  ggsave(paste0(clean_name, "_celltype2.pdf"), 
         plots$violin, width = 5, height = 3, dpi = 600)
  ggsave(paste0(clean_name, "_tsne2.pdf"), 
         plots$tsne, width = 7, height = 6, dpi = 600)
  
  # 统计结果
  stats <- Late@meta.data %>%
    group_by(celltype) %>%
    summarise(
      mean = mean(get(paste0(clean_name, "_1"))),
      median = median(get(paste0(clean_name, "_1"))),
      sd = sd(get(paste0(clean_name, "_1")))
    ) %>%
    arrange(desc(median))
  
  list(plots = plots, stats = stats)
})

# 6. 结果保存 -------------------------------------------------------------

# 6.1 保存统计结果
write_xlsx(
  lapply(results, function(x) x$stats) %>% set_names(names(all_genesets)),
  "all_geneset_scores.xlsx"
)

# 6.2 保存匹配详情
write_xlsx(
  list(
    Summary = match_summary,
    Missing_Genes = stack(lapply(match_results, function(x) paste(x$missing, collapse = ", ")))
  ),
  "geneset_matching_details.xlsx"
)

# 6.3 保存Seurat对象
saveRDS(Late, "Late_with_geneset_scores.rds")

# 6.4 生成匹配情况图
ggplot(match_summary, aes(x = reorder(Geneset, -Percentage), y = Percentage)) +
  geom_col(fill = "#56B4E9") +
  geom_text(aes(label = paste0(Found, "/", Total)), vjust = -0.5) +
  labs(title = "Gene Set Matching Rate", x = NULL, y = "Percentage (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("geneset_matching_rate.pdf", width = 8, height = 6)

cat("\n=== 分析完成 ===\n")
cat("结果已保存到：\n",
    "- all_geneset_scores.xlsx\n",
    "- geneset_matching_details.xlsx\n",
    "- Seurat_with_geneset_scores.rds\n",
    "- 各基因集的PDF图表\n")
