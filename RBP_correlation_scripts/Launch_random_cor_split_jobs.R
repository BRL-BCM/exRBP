args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)
totalN=c("BCLAF1", "PRPF8", "DDX3X", "AQR", "GRWD1", "CSTF2T", "PPIG", "KHSRP", "HNRNPM", "SUGP2", "SF3B4", "ZNF622", "UPF1", "EXOSC5", "QKI", "PTBP1", "U2AF2", "PUM2", "PRPF4", "ILF3", "BUD13", "HNRNPL", "DDX24", "SND1", "HNRNPK", "MATR3")
totalN=c("BUD13", "DDX24", "EFTUD2", "EXOSC5", "FAM120A", "FXR2", "HNRNPK", "IGF2BP1", "LIN28B", "MATR3", "NCBP2", "PABPC4", "PRPF4", "PUM2", "RBFOX2", "RPS3", "XPO5", "ZNF800")
for (n in totalN) {
  files_exists=list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/csf",pattern = paste("*",n,"*",sep=""),full.names=T)
  RBP_pairs_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/breaks_E/CSF/",n,"_breaks.txt.gz",sep="")
  rows_use=read.table(RBP_pairs_files,sep="\t",header=T)
  for(subset in 2:dim(rows_use)[1]) {
    test=paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/csf/",n,"_corComplete_subset",rows_use[subset,1],".txt.gz",sep="")
    if (! test%in%files_exists) {
      str = sprintf("qsub -v RBP=%s,start=%i,end=%i -N %s random_cor_split_jobs.pbs", n, rows_use[subset,1], rows_use[subset,2], paste(n,"randomCor",rows_use[subset,1],sep="_") )
      system(str, wait=FALSE)
      system("sleep 0.5")
      #print(paste(n,"_corComplete_subset",rows_use[[subset]][1],".txt",sep=""))
    }
  }
  print(paste(n,"_corComplete_subset",rows_use[subset,1],".txt",sep=""))
}

