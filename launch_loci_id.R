
#totalN = c("EXR-DGALA1CSFBQ8113-BS", "EXR-DERLE1Dj6liR-AN")
totalN = as.character(read.table("../RBP.done.txt",sep="\t")[,1])

for (n in totalN) {
  str = sprintf("qsub -N %s -v study=%s study_loci_presence.pbs", paste(n,"countLoci",sep="_"), n )
  system(str, wait=FALSE)
  system("sleep 0.2")
}

