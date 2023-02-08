args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)

#totalN=c("KHSRP","PRPF8")
#totalN=c("KHSRP","AQR","PRPF8","BCLAF1")
#totalN="DHX30"
totalN=c("AKAP1", "BCLAF1", "DDX3X", "GRWD1", "HNRNPK", "IGF2BP3", "LIN28B", "NIPBL", "PPIL4", "PRPF8", "QKI", "SAFB2", "SRSF1", "SUB1", "XPO5", "AQR", "CSTF2T", "FKBP4", "HLTF", "HNRNPM", "KHSRP", "MTPAP", "POLR2G", "PRPF4", "PUM1", "RPS11", "SERBP1", "STAU2", "WDR3", "YWHAG")

totalN=c("SRSF7", "SF3A3", "GPKOW", "NONO", "QKI", "LARP4", "XPO5", "XRN2", "NIPBL", "NSUN2","KHSRP", "RBM5", "SUPV3L1", "FXR2","HNRNPL", "ZRANB2", "NOLC1", "DDX51", "PUM2", "TROVE2", "FKBP4")
for (n in totalN) {
    str = sprintf("qsub -v RBP=%s -N %s combine_random_cor_split_jobs2.pbs", n, paste(n,"collapse",sep="_") )
    system(str, wait=FALSE)
    system("sleep 0.5")
}

