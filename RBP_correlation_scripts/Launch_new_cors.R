args = commandArgs(trailingOnly=TRUE)
sample = args[2]
# getwd()

totalN = list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/CSF",pattern = "*_RBP_subset.txt.gz")
totalN = gsub("_RBP_subset.txt.gz","",totalN)
done =  list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/pairSelected/CSF",pattern="*_pairSelected_subset.txt.gz")
done=gsub("_pairSelected_subset.txt.gz","",done)
totalN=totalN[!totalN%in%done]
#totalN=totalN[6:length(totalN)]
#totalN=totalN[1:5]
#totalN=c("AQR","PRPF8","BCLAF1")

for (n in totalN) {
  str = sprintf("qsub -v RBP=%s -N %s Create_correlatedPairs.pbs", n,  paste(n,"withinCor_CSF",sep="_") )
  system(str, wait=FALSE)
  system("sleep 0.5")
}

