#devtools::install_github("mrcieu/gwasglue",force = TRUE)
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")

#install.packages("R.utils")


#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposureFile="exposure.F.csv"             #暴露数据文件
outcomeFile="finngen_R12_DM_NEPHROPATHY_EXMORE.gz"     #结局数据文件(需修改)
outcomeName="Diabetic Nephrology"           #图形中展示疾病的名称
setwd("G:\\187geneMR\\12.MR\\finngen_R12_DM_NEPHROPATHY_EXMORE")     #设置工作目录
#Zhou W(199+388756),PMID30104761
#读取暴露数据
exposure_dat=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

# 查看前几行数据（不加载全部）
#outcome_lines <- readLines(gzfile(outcomeFile), n = 5)
#cat(outcome_lines, sep = "\n")

# 或用 data.table 快速查看（适用于大文件）
library(data.table)
gwas_data <- fread(outcomeFile, nrows = 5)
colnames(gwas_data)

#读取结局数据文件
outcomeData=read_outcome_data(snps=exposure_dat$SNP,
                 filename=outcomeFile, sep = "\t",
                 snp_col = "rsids",
                 beta_col = "beta",
                 se_col = "sebeta",
                 effect_allele_col = "alt",
                 other_allele_col = "ref",
                 pval_col = "pval",
                 eaf_col = "af_alt")
write.csv(outcomeData, file="finngen_R12_DM_NEPHROPATHY_EXMORE.csv", row.names=F)
head(outcomeData)

#outcomeData <- read.csv("GCST90435706_outcome.csv")
#head(outcomeData_reloaded)  # 检查是否和之前一样

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcomeData)
dat=dat[dat$pval.outcome>5e-06,]
#dat <- dat[dat$pval.outcome > 1e-05, ]  # 更宽松,可能增加多效性风险

#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="table.SNP.csv", row.names=F)
nrow(dat) 
#孟德尔随机化分析
mr_method_list()
#mrResult=mr(dat)
mrResult <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="table.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file="table.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file="table.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="pic.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="pic.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="pic.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="pic.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()



