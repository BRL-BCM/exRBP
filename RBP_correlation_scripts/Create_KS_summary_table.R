library(data.table)

RBP_files = list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines",
                          pattern = "_pValTable_subset.txt",full.names = T)
#RBP_files = list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/cov5_sample30",
#                          pattern = "_pValTable_cov30.txt",full.names = T)

output=data.frame(p.value=NA,
                  d.value=NA)

rownames_use="x"

for (RBP in 1:length(RBP_files)) {
  RBP_in=as.data.frame(fread(RBP_files[RBP],sep="\t",header=T))
  name = gsub("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/","",RBP_files[RBP])
#  name = gsub("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/cov5_sample30/","",RBP_files[RBP])
  name = as.character(gsub("_pValTable_subset.txt","",name))
#  name = as.character(gsub("_pValTable_cov30.txt","",name))
 
  RBP_in$FDR=RBP_in$p.value/dim(RBP_in)[1]
  RBP_in$FDR_other=RBP_in$random_p.value/dim(RBP_in)[1]
  KS.p.val=ks.test(RBP_in$FDR,RBP_in$FDR_other)$p.value
  KS.d.val=as.numeric(ks.test(RBP_in$FDR,RBP_in$FDR_other)$statistic)
  
  output=rbind.data.frame(output,
                          c(KS.p.val,KS.d.val))
  rownames_use=c(rownames_use,name)

  pdf(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/pdfs/",name,"_pVal.pdf",sep=""))
  plot(ecdf(RBP_in$p.value), xlim = range(c(RBP_in$p.value, RBP_in$random_p.value)), lty = "dashed",col="green",main=name)
  plot(ecdf(RBP_in$random_p.value), add = TRUE,col="red")
  dev.off()
  
  #pdf(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/cov5_sample30/pdfs/",name,"_pVal.pdf",sep=""))
  #plot(ecdf(RBP_in$p.value), xlim = range(c(RBP_in$p.value, RBP_in$random_p.value)), lty = "dashed",col="green",main=name)
  #plot(ecdf(RBP_in$random_p.value), add = TRUE,col="red")
  #dev.off()

}

rownames(output)=rownames_use
output=output[-1,]

write.table(output,"/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_KS_pVal_results_withFDR_CellLines.txt",
            sep="\t",quote=F)

#write.table(output,"/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_KS_pVal_results_withFDR_CellLines_cov5_sample_30.txt",
#            sep="\t",quote=F)

