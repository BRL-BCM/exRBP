#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(data.table)

RBP_in_files = list.files(path=paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",args[1],"/cov5_sample30/",args[3],"perc_cov",args[2],sep=""),
#                          pattern = "_pValTable_subset.txt",full.names = T)
                          pattern = "_pValTable_cov30.txt.gz",full.names = T)

output=data.frame(p.value=NA,
                  d.value=NA)

rownames_use="x"

for (RBP in 1:length(RBP_in_files)) {
  RBP_in=as.data.frame(fread(RBP_in_files[RBP],sep="\t",header=T))
  name = gsub(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",args[1],"/cov5_sample30/",args[3],"perc_cov",args[2],"/",sep=""),"",RBP_in_files[RBP])
#  name = as.character(gsub("_pValTable_subset.txt","",name))
  name = as.character(gsub("_pValTable_cov30.txt.gz","",name))
 
  RBP_in$FDR=RBP_in$p.value/dim(RBP_in)[1]
  RBP_in$FDR_other=RBP_in$random_p.value/dim(RBP_in)[1]
  KS.p.val=ks.test(RBP_in$FDR,RBP_in$FDR_other)$p.value
  KS.d.val=as.numeric(ks.test(RBP_in$FDR,RBP_in$FDR_other)$statistic)
  
  output=rbind.data.frame(output,
                          c(KS.p.val,KS.d.val))
  rownames_use=c(rownames_use,name)
  
  pdf(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",args[1],"/cov5_sample30/",args[3],"perc_cov",args[2],"/pdfs/",name,"_pVal.pdf",sep=""))
  plot(ecdf(RBP_in$p.value), xlim = range(c(RBP_in$p.value, RBP_in$random_p.value)), lty = "dashed",col="green",main=name)
  plot(ecdf(RBP_in$random_p.value), add = TRUE,col="red")
  dev.off()
}

rownames(output)=rownames_use
output=output[-1,]

gz1=gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_KS_pVal_results_withFDR_",args[1],"_cov",args[2],"_sample_",args[3],".txt.gz",sep=""),"w")
write.table(output,gz1,
            sep="\t",quote=F)
close(gz1)

