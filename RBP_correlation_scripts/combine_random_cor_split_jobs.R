args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)

#totalN=c("KHSRP","AQR","PRPF8","BCLAF1")
#totalN="DHX30"

RBP_pairs_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/breaks_E/CellLines/",RBP,"_breaks.txt",sep="")
rows_use=read.table(RBP_pairs_files,sep="\t",header=T)
print(dim(rows_use))
if (dim(rows_use)[1]>2){
  subset_loc=c()
  for(subset in 2:dim(rows_use)[1]) {
    subset_loc=c(subset_loc,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/CellLines/",RBP,"_corComplete_subset",rows_use[subset,1],".txt",sep=""))
  } 

  total_file=as.data.frame(fread(subset_loc[1],sep="\t",header=T))
  for (subset_file in subset_loc[2:length(subset_loc)]) {
    print(subset_file)
    new_file=as.data.frame(fread(subset_file,sep="\t",header=T))
    total_file=rbind.data.frame(total_file,new_file[-1,])
  } 

  print("done")
  write.table(total_file,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/CellLines/",RBP,"_corComplete_subset.txt",sep=""),quote=F,row.names=F,sep="\t")
} else {
  files=list.files("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/CellLines/",pattern=paste("^",RBP,sep=""),full.names = T)
  tmp=as.data.frame(fread(files[1],sep="\t",header=T))
  write.table(tmp,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/CellLines/",RBP,"_corComplete_subset.txt",sep=""),quote=F,row.names=F,sep="\t")
}
