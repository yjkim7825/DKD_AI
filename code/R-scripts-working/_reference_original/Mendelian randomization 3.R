library(forestploter)
library(grid)  # 提供gpar()函数
library(forestplot) 
# 参数设置
selectMethod = c("Inverse variance weighted")
setwd("G:\\187geneMR\\16.forest\\GSE96804_104948")
files = dir()
files = grep("csv$", files, value = TRUE)

# 读取并筛选数据
data = data.frame()
for(i in files){
  rt = read.csv(i, header = TRUE, sep = ",", check.names = FALSE)
  rt = rt[rt$method %in% selectMethod, ]  # 提前筛选方法
  data = rbind(data, rt)
}

# 整理数据
data$' ' <- paste(rep(" ", 10), collapse = " ")
data$'OR(95% CI)' = ifelse(is.na(data$or), "", sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))
data$pval = ifelse(data$pval < 0.001, "<0.001", sprintf("%.3f", data$pval))
data$exposure = ifelse(is.na(data$exposure), "", data$exposure)
data$nsnp = ifelse(is.na(data$nsnp), "", data$nsnp)

# 安全地处理重复的exposure
if(any(duplicated(data$exposure))) {
  data$exposure[duplicated(data$exposure)] <- ""
}

# 图形主题
tm <- forest_theme(
  base_size = 18,
  ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2,
  refline_lty = "dashed", refline_lwd = 1, refline_col = "grey20",
  xaxis_cex = 0.8,
  footnote_cex = 0.6, footnote_col = "blue"
)

# 绘制图
# 创建显示用的数据框时转换列名
display_data <- data[, c("exposure","nsnp","method","pval"," ","OR(95% CI)")]

# 将列名首字母大写
colnames(display_data) <- c("Exposure", "Nsnp", "Method", "Pval", " ", "OR(95% CI)")

# 绘制图形（使用修改后的列名）
plot <- forest(
  display_data,  # 使用修改后的数据框
  est = data$or,
  lower = data$or_lci95,
  upper = data$or_uci95,
  ci_column = 5,  
  ref_line = 1,
  xlim = c(0, 2),
  theme = tm
)
# 设置颜色（循环使用颜色）
boxcolor <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4",
              "#91D1C2","#DC0000","#7E6148","#E64B35","#4DBBD5","#00A087",
              "#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")

# 正确应用置信区间颜色
for(i in 1:nrow(data)){
  plot <- edit_plot(plot, 
                    col = 5,  # CI列的位置
                    row = i, 
                    which = "ci", 
                    gp = gpar(fill = boxcolor[(i-1)%%length(boxcolor)+1], 
                              col = boxcolor[(i-1)%%length(boxcolor)+1]))
}

# 统一设置method列字体（不加粗）
plot <- edit_plot(plot, col = 3, row = 1:nrow(data), which = "text", 
                  gp = gpar(fontface = "plain"))

# 标记显著p值（只加粗pval列）
pos_bold_pval = which(as.numeric(gsub('<', "", data$pval)) < 0.05)
if(length(pos_bold_pval) > 0){
  for(i in pos_bold_pval){
    plot <- edit_plot(plot, col = 4, row = i, which = "text",  # pval列现在是第4列
                      gp = gpar(fontface = "bold"))
  }
}

# 调整图形格式
plot <- add_border(plot, part = "header", row = 1, where = "top", gp = gpar(lwd = 2))
plot <- edit_plot(plot, col = 1:ncol(data), row = 1:nrow(data), which = "text", 
                  gp = gpar(fontsize = 12))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text", 
                  hjust = unit(0.5, "npc"), part = "header", x = unit(0.5, "npc"))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text", 
                  hjust = unit(0.5, "npc"), x = unit(0.5, "npc"))

# 输出图形
pdf("forest.pdf", width = 10, height = 8)
print(plot)
dev.off()