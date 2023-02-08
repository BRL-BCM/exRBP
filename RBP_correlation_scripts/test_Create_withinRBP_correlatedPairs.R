args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)
library(Hmisc)

#From http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}


RBP_regions_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/CSF/",RBP,"_RBP_subset.txt.gz",sep="")
print(RBP_regions_files)
RBP_rows=fread(RBP_regions_files,sep="\t",header=T)

#RBP_regions_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/Cell_line_all/RBP_loci_each/",RBP,"_loci_present.txt",sep="")
#RBP_regions=read.table(RBP_regions_files,sep="\t",header=T)[,1]

print("orig")
df_DERLE = fread("CSF_healthy_filtered.txt",sep="\t",header=T,fill=T) ##cellLine_all_filtered.txt for cell lines
totalNum=dim(df_DERLE)[2]

cols=c("chrom","start",'end')
df_DERLE[,loci:=Reduce(function(...) paste(..., sep = "_"), .SD[, mget(cols)])]
#print("new")
#print(head(df_DERLE))

#print("here2")
dtnew <- df_DERLE[, lapply(.SD, as.integer), .SDcols=4:totalNum]
dt_use=as.data.frame(cbind(df_DERLE[,.(loci)],dtnew))
print(head(dtnew[,1:3]))
print("next")
dtnew[, SUM := Reduce(`+`, .SD,init=1), .SDcols=1:dim(dtnew)[2]]
print(head(dtnew))

#print("3")
#print(head(df_DERLE[,.(loci)]))
#print("4")
#print(head(dtnew[,.(SUM)]))

print("5")
df_DERLE_sum=as.data.frame(cbind(df_DERLE[,.(loci)],dtnew[,.(SUM)]))
print(length(df_DERLE_sum$SUM))
print(length(df_DERLE_sum$SUM[is.na(df_DERLE_sum$SUM)]))

df_DERLE_sum=df_DERLE_sum[!is.na(df_DERLE_sum$SUM),]

print(head(df_DERLE_sum))

rm(dtnew)
rm(df_DERLE)

RBP_t=transpose(RBP_rows,make.names = "loci")
res2<-rcorr(as.matrix(RBP_t))

print("test")
print(head(res2))
rm(RBP_t)

output=flattenCorrMatrix(res2$r)
output$total_reads_lociB=rep(0,dim(output)[1])
output$replacement_loci=rep("None",dim(output)[1])

print("test2")
print(head(df_DERLE_sum$loci))
print(head(RBP_rows$loci))

dt_study=df_DERLE_sum[df_DERLE_sum$loci%in%RBP_rows$loci,]
dt_study$replacement_loci=rep("None",dim(dt_study)[1])

print(head(df_DERLE_sum))
use_Df=df_DERLE_sum[!df_DERLE_sum$loci%in%c(output$row,output$column),]

use_Df=df_DERLE_sum[!df_DERLE_sum$loci%in%dt_study$loci,]


for (lociB in unique(output$column)) {
  #print(lociB)
  coverage = dt_study[dt_study$loci==lociB,2]
  
  possible=use_Df[use_Df$SUM==coverage,1]
  
  #print(head(use_Df))
  # print(head(possible))
  #print(range(use_Df[,2]))
  # print(coverage)
  
  if (length(possible)<1 && coverage<max(use_Df[,2])) {
    possible_df=use_Df[use_Df$SUM>coverage,]
    possible_df=possible_df[order(possible_df$SUM,decreasing = FALSE),]
    possible=possible_df[1,1]
  }
  if (length(possible)<1 && coverage>max(use_Df[,2])) {
    coverage_use=max(use_Df[,2])
    possible=use_Df[use_Df$SUM==coverage_use,1]
  }
  
  selected=sample(possible,1)
  rm(possible)
  
  if (length(coverage)>1) {
    coverage=coverage[1]
  }
  
  output[output$column==lociB,4]=coverage
  output[output$column==lociB,5]=selected
  use_Df=use_Df[use_Df$loci!=selected,]
  
}

#write.table(use_Df,"/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/df_DERLE.txt",sep="\t",quote=F,row.names=F)
gz1=gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/pairSelected/CSF/",RBP,"_pairSelected_subset.txt.gz",sep=""),"w")
write.table(output,gz1,
            sep="\t",quote=F,row.names = F)
close(gz1)

#lociA=as.character(output[,1])
#lociA=lociA[!duplicated(lociA)]
#lociB=output[,5]
#lociB=lociB[!duplicated(lociB)]

#RBP_t=transpose(dt_use[dt_use$loci%in%c(lociA,lociB),],make.names = "loci")
#Each column represents the correlation to all loci in that RBP by a random loci
#output_cor=cor(RBP_t[,colnames(RBP_t)%in%lociA],
#               RBP_t[,colnames(RBP_t)%in%lociB])

#output$otherCor=rep("None",dim(output)[1])

#for (row_loci in 1:dim(output)[1]) {
#    output[row_loci,6]=output_cor[output[row_loci,1],output[row_loci,5]]
#}

#write.table(output,"/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/DGCR8_pairSelected_cor.txt",sep="\t",quote=F)
#write.table(output,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/",RBP,"_corComplete_subset.txt",sep=""),
#            sep="\t",quote=F, row.names = F)
