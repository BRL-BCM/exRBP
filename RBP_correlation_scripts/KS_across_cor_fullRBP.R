library(data.table)

RBP_in_files = list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/PHT",
                          pattern = "_corComplete_subset.txt",full.names = T)


output=data.frame(p.value=NA,
                  d.value=NA)

rownames_use="x"

for (RBP in 1:length(RBP_in_files)) {
  RBP_in=as.data.frame(fread(RBP_in_files[RBP],sep="\t",header=T))
  name = gsub("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/PHT/","",RBP_in_files[RBP])
  name = as.character(gsub("_corComplete_subset.txt","",name))
  
  KS.p.val=ks.test(RBP_in$cor,RBP_in$otherCor)$p.value
  KS.d.val=as.numeric(ks.test(RBP_in$cor,RBP_in$otherCor)$statistic)
  
  output=rbind.data.frame(output,
                          c(KS.p.val,KS.d.val))
  rownames_use=c(rownames_use,name)
  
  #pdf(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/PHT/pdfs/",name,".pdf",sep=""))
  #plot(ecdf(RBP_in$p.value), xlim = range(c(RBP_in$p.value, RBP_in$random_p.value)), lty = "dashed",col="green",main=name)
  #plot(ecdf(RBP_in$random_p.value), add = TRUE,col="red")
  #dev.off()
}

rownames(output)=rownames_use
output=output[-1,]
output=cbind(output,output[,1]/dim(output)[1])

write.table(output,"/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_KS_overCor_withFDR_PHT.txt",
            sep="\t",quote=F)

