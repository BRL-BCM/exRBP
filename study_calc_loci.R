filenames_bed = list.files(path = ".", pattern = "*.bed$", full.names = TRUE)
print(getwd())
print(filenames_bed[1])

Loci_present=c()

for (i in filenames_bed) {
  bed_sample = read.table(i,sep="\t",header=T)
  bed_sample=bed_sample[,c(1:7, 10, 13:15, 18:20, 23:26, 29:32, 35:37, 40:43, 46, 49, 52, 53, 56, 59, 60:64, 67, 70, 73, 74, 75, 78, 79, 82, 83, 86:91, 94, 97, 100, 103, 106, 109, 112, 115, 118, 121, 124, 125, 126, 129, 130, 133, 136, 139, 142, 145, 148:150, 153:157, 160:165, 168:174, 177, 180:183, 186, 189, 192, 195:197, 200, 201, 204, 205, 208:211, 214:216, 219, 222, 225, 228, 231, 232, 235:238, 241, 244, 245, 248, 251, 252, 255, 258, 261, 264, 267, 270, 273:275, 278:280, 283, 286, 289, 290, 293:295, 298)]
  bed_sample=bed_sample[rowSums(bed_sample[,4:154])>0,]
  loci_present_sample=paste(bed_sample$chrom,bed_sample$start,bed_sample$end,sep="_")
  Loci_present=c(Loci_present,loci_present_sample)
  Loci_present=Loci_present[!duplicated(Loci_present)]
}


write.table(Loci_present,"Loci_in_study.txt",sep="\t",quote=F)
