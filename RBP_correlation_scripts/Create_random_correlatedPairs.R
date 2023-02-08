args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)

RBP_pairs_files = paste(RBP,"_pairSelected_subset.txt",sep="")

df_DERLE = fread("urine_healthy_filtered.txt",sep="\t",header=T,fill=T)
#df_DERLE = fread("cellLine_all_filtered.txt",sep="\t",header=T,fill=T)
totalNum=dim(df_DERLE)[2]

cols=c("chrom","start",'end')
df_DERLE[,loci:=Reduce(function(...) paste(..., sep = "_"), .SD[, mget(cols)])]

dtnew <- df_DERLE[, lapply(.SD, as.integer), .SDcols=4:totalNum]

print("before")
dtnew=cbind(df_DERLE[,.(loci)],dtnew)
rm(df_DERLE)
print("after")

f=function(row_in,dt){
  pairA=as.character(row_in[1])
  pairB=as.character(row_in[5])
  
  rows_dt=dt[c(which(dt$loci==pairA),
               which(dt$loci==pairB)),2:dim(dt)[2]]
  correlation=cor(as.numeric(rows_dt[1,]),
                  as.numeric(rows_dt[2,]))
  return(correlation)
}

print("correlation start")
RBP_cor=fread(RBP_pairs_files,sep="\t",header=T)

new_Cor=apply(RBP_cor, 1, f, dt=dtnew)
RBP_cor$otherCor=new_Cor

#write.table(RBP_cor,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/Urine/",RBP,"_corComplete_subset.txt",sep=""),
#            sep="\t",quote=F, row.names = F)
gz1 = gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/Urine/",RBP,"_corComplete_subset.txt.gz",sep=""),"w")
write.table(RBP_cor,gz1,
            sep="\t",quote=F, row.names = F)
close(gz1)


