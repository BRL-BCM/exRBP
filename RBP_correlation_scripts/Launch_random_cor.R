args = commandArgs(trailingOnly=TRUE)
sample = args[2]
# getwd()

totalN = list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/urine",pattern = "*_RBP_subset.txt",full.names = T)
totalN = gsub("_RBP_subset.txt","",totalN)
totalN = gsub("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/urine/","",totalN)

totalN=c("EFTUD2", "LIN28B", "FXR2", "RPS3", "HNRNPC", "NCBP2", "FAM120A", "IGF2BP1", "EWSR1", "XPO5", "SLTM", "AGGF1", "FUS", "PABPC4", "AKAP1", "EIF4G2", "PABPN1", "NKRF", "SUB1", "ZNF800", "HLTF", "EIF3H", "LARP4", "DROSHA", "RBM15", "CSTF2", "FMR1", "ZC3H11A", "YBX3", "GEMIN5", "LSM11", "IGF2BP3", "GTF2F1", "PCBP1", "SRSF1", "TRA2A", "DGCR8", "SF3A3", "NOLC1", "DDX55", "GRSF1", "RBM22", "GPKOW", "XRN2", "UCHL5", "TIA1", "XRCC6", "DDX6", "AKAP8L", "TAF15", "CPSF6", "DHX30", "FTO", "TBRG4", "U2AF1", "PPIL4", "APOBEC3C", "MTPAP", "TROVE2", "LARP7", "SDAD1", "SLBP", "FASTKD2", "CDC40", "BCCIP", "DDX52", "AATF", "SSB", "HNRNPU", "SAFB", "HNRNPUL1", "EIF3G", "SRSF7", "CPEB4", "NIPBL", "DDX51", "NSUN2", "SMNDC1", "PHF6", "NOL12", "WDR43", "DDX21", "HNRNPA1", "SUPV3L1", "POLR2G", "DDX42", "RPS11")
done=list.files("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/Urine")
done=gsub("_corComplete_subset.txt.gz","","done")

totalN=rev(totalN)
for (n in totalN) {
  if(!n%in%done){
    str = sprintf("qsub -v RBP=%s -N %s Create_random_correlatedPairs.pbs", n,  paste(n,"randomCor_urine",sep="_") )
    system(str, wait=FALSE)
    system("sleep 0.5")
  }
}

