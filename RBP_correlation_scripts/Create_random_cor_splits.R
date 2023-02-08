library(data.table)

totalN=c("BCLAF1", "PRPF8", "DDX3X", "AQR", "GRWD1", "CSTF2T", "PPIG", "KHSRP", "HNRNPM", "SUGP2", "SF3B4", "ZNF622", "UPF1", "EXOSC5", "QKI", "PTBP1", "U2AF2", "PUM2", "PRPF4", "ILF3", "BUD13", "HNRNPL", "DDX24", "SND1", "HNRNPK", "MATR3", "RBFOX2")
totalN=c("BUD13", "DDX24", "EFTUD2", "EXOSC5", "FAM120A", "FXR2", "HNRNPK", "IGF2BP1", "LIN28B", "MATR3", "NCBP2", "PABPC4", "PRPF4", "PUM2", "RBFOX2", "RPS3", "XPO5", "ZNF800")
for (n in totalN) {
  RBP_pairs_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/pairSelected/CSF/",n,"_pairSelected_subset.txt.gz",sep="")
  RBP_cor=fread(RBP_pairs_files,sep="\t",header=T)
  #rows_use=split(1:dim(RBP_cor)[1], ceiling(seq_along(2:dim(RBP_cor)[1])/3000))
  rows_use=split(2:dim(RBP_cor)[1], ceiling(seq_along(2:dim(RBP_cor)[1])/30000))
  output_table=data.frame(start=NA,
 			  end=NA)
  for(subset in 1:length(rows_use)) {
    output_table=rbind.data.frame(output_table,
				c(rows_use[[subset]][1], rows_use[[subset]][length(rows_use[[subset]])]))
  }
  gz1 = gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/breaks_E/CSF/",n,"_breaks.txt.gz",sep=""),"w")
  write.table(output_table,gz1,sep="\t",quote=F)
  close(gz1)
}


