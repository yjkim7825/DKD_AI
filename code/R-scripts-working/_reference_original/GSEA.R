# 加载必要的包
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 设置工作目录
setwd("G:\\187geneMR\\08.GSEA")

# 定义基因集文件路径
gmtFile <- "c2.cp.kegg.Hs.symbols.gmt"

# 读取三组差异基因结果文件
early_vs_control <- read.table("diff_Early_vs_Control.txt", header=TRUE, row.names=1, sep="\t")
late_vs_control <- read.table("diff_Late_vs_Control.txt", header=TRUE, row.names=1, sep="\t")
late_vs_early <- read.table("diff_Late_vs_Early.txt", header=TRUE, row.names=1, sep="\t")

# 定义执行GSEA分析的函数
run_gsea_analysis <- function(diff_data, comparison_name, gmt, min_terms=3) {
  # 准备logFC向量
  logFC <- diff_data$logFC
  names(logFC) <- rownames(diff_data)
  logFC <- sort(logFC, decreasing = TRUE)
  
  # 执行GSEA分析
  message("\nRunning GSEA for: ", comparison_name)
  kk <- GSEA(logFC, TERM2GENE = gmt, pvalueCutoff = 1)
  
  # 保存完整结果
  kkTab <- as.data.frame(kk)
  write.table(kkTab, 
              file = paste0("GSEA.result.", comparison_name, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 筛选显著结果 (p.adjust < 0.05)
  kkTab_sig <- kkTab[kkTab$p.adjust < 0.05, ]
  message("Significant pathways found: ", nrow(kkTab_sig))
  
  # 可视化设置
  plot_pathways <- function(kk, pathways, plot_title, file_suffix) {
    if (length(pathways) > 0) {
      gseaplot <- gseaplot2(kk, pathways, base_size = 10, 
                            title = plot_title,
                            color = "firebrick")
      ggsave(paste0("GSEA.", comparison_name, ".", file_suffix, ".pdf"),
             plot = gseaplot, width = 8, height = 6)
      message("Generated plot: ", file_suffix)
    } else {
      message("No pathways to plot for ", file_suffix)
    }
  }
  
  # 高表达组富集通路
  up_terms <- rownames(kkTab_sig[kkTab_sig$NES > 0, ])
  plot_pathways(kk, 
                pathways = head(up_terms, min_terms),
                plot_title = paste0(comparison_name, " - Up-regulated Pathways"),
                file_suffix = "Up")
  
  # 低表达组富集通路
  down_terms <- rownames(kkTab_sig[kkTab_sig$NES < 0, ])
  plot_pathways(kk, 
                pathways = head(down_terms, min_terms),
                plot_title = paste0(comparison_name, " - Down-regulated Pathways"),
                file_suffix = "Down")
  
  return(kkTab_sig)
}

# 读取基因集文件
gmt <- read.gmt(gmtFile)

# 对三个比较组分别运行GSEA分析（设置min_terms=3降低阈值）
early_vs_control_gsea <- run_gsea_analysis(early_vs_control, "Early_vs_Control", gmt, min_terms=2)
late_vs_control_gsea <- run_gsea_analysis(late_vs_control, "Late_vs_Control", gmt, min_terms=3)
late_vs_early_gsea <- run_gsea_analysis(late_vs_early, "Late_vs_Early", gmt, min_terms=3)

# 合并所有显著结果
all_gsea_results <- rbind(
  cbind(Comparison = "Early_vs_Control", early_vs_control_gsea),
  cbind(Comparison = "Late_vs_Control", late_vs_control_gsea),
  cbind(Comparison = "Late_vs_Early", late_vs_early_gsea)
)

# 保存合并后的结果
write.table(all_gsea_results, 
            file = "All_GSEA_Results_Combined.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\nAnalysis completed! Check output files in: ", getwd())