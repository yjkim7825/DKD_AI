
mrFile="table.MRresult.csv"        #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"     #多效性的结果文件
setwd("G:\\187geneMR\\13.IVWfilter\\finngen_R12_DM_NEPHROPATHY_EXMORE")     #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)
#提取IVW方法pvalue<0.05的基因
ivw=data.frame()
for(geneName in unique(rt$exposure)){
	geneData=rt[rt$exposure==geneName,]
	#提取3种方法OR方向一致的基因
	if(nrow(geneData)==3){
		if(geneData[geneData$method=="Inverse variance weighted","pval"]<0.05){
			if(sum(geneData$or>1)==nrow(geneData) | sum(geneData$or<1)==nrow(geneData)){
				ivw=rbind(ivw, geneData)
			}
		}
	}
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除多效性pvalue小于0.05的基因
pleRT=pleRT[pleRT$pval>0.05,]
immuneLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% immuneLists,]
write.csv(outTab, file="IVW.filter.csv", row.names=F)
