args = commandArgs(trailingOnly=TRUE)
sample = args[2]
# getwd()

totalN=list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/csf",pattern = "*_corComplete_subset.txt.gz",full.names = T)
totalN = gsub("_corComplete_subset.txt.gz","",totalN)
totalN = gsub("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/csf/","",totalN)
#totalN=c("UPF1","PRPF8","AQR","BCLAF1","CSTF2T","YBX3")

done=list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CSF",pattern="*_pValTable_subset.txt.gz")
done=gsub("_pValTable_subset.txt.gz","",done)

totalN=totalN[!totalN%in%done]
for (n in totalN) {
  str = sprintf("qsub -v RBP=%s -N %s Create_KS_values.pbs", n,  paste(n,"KS_csf",sep="_") )
  system(str, wait=FALSE)
  system("sleep 0.5")
}

