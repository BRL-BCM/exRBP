args = commandArgs(trailingOnly=TRUE)

totalN = list.files(path=paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",args[1],sep=""),pattern = "*_pValTable_subset.txt.gz")
totalN = gsub("_pValTable_subset.txt.gz","",totalN)

done=list.files(path=paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",args[1],"/cov5_sample30/",args[3],"perc_cov",args[2],sep=""),
                pattern = "*_pValTable_cov30.txt.gz")
done = gsub("_pValTable_cov30.txt.gz","",done)

totalN=totalN[!totalN%in%done]
#totalN=c("UPF1","PRPF8","AQR","BCLAF1","CSTF2T","YBX3")

for (n in totalN) {
  str = sprintf("qsub -v RBP=%s,SET=%s,BIO=%s,Cor=%s -N %s Create_KS_values_subset.pbs", n, paste(args[3],"perc_cov",args[2],sep=""), args[1], args[4], paste(n,"KS_sub",args[1],sep="_") )
  system(str, wait=FALSE)
  system("sleep 0.5")
}

