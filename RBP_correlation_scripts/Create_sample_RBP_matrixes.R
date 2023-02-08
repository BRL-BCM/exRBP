library(data.table)

RBP_regions_files = list.files(path="/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/Cell_line_all/RBP_loci_each",pattern = "*_new.txt",full.names = T)
#df_DERLE = fread("serum_summed_filtered.txt",sep="\t",header=T,fill=T)
set=list.files(pattern = "*_filtered.txt")
df_DERLE = fread(set,sep="\t",header=T,fill=T)
#df_DERLE = fread("serum_healthy_filtered.txt",sep="\t",header=T,fill=T)

print(dim(df_DERLE))
print(colnames(df_DERLE))
print(head(df_DERLE[,1:5]))
totalNum=dim(df_DERLE)[2]

cols=c("chrom","start",'end')
df_DERLE[,loci:=Reduce(function(...) paste(..., sep = "_"), .SD[, mget(cols)])]

dtnew <- df_DERLE[, lapply(.SD, as.integer), .SDcols=4:totalNum]
dtnew[, SUM := Reduce(`+`, .SD,init=1), .SDcols=1:dim(dtnew)[2]]

dtnew=cbind(df_DERLE[,.(loci)],dtnew)
rm(df_DERLE)

for (RBP in RBP_regions_files) {
  RBP_regions=read.table(RBP,sep="\t",header=T)[,1]
  RBP_short=dtnew[dtnew$loci%in%RBP_regions,]
  print(dim(RBP_short)[1])
  if (dim(RBP_short)[1]>4){
   name = gsub("_new.*","",RBP)
   name = gsub("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/Cell_line_all/RBP_loci_each/","",name)
   print(name) 
   write.table(RBP_short,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/serum/",name,"_RBP_subset.txt",sep=""),
              sep="\t",quote=F, row.names = F)
  }
  rm(RBP_regions)
  rm(RBP_short)
  
}
