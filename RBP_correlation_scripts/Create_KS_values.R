args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)
library(Hmisc)
#library(sm)

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
# RBP_rows=fread("DGCR8_RBP_subset.txt",sep="\t",header=T)

#df_DERLE = fread("EXR-DERLE1Dj6liR-AN_filtered.txt",sep="\t",header=T,fill=T)
#df_DERLE = fread("cellLine_all_filtered.txt",sep="\t",header=T,fill=T)
df_DERLE = fread("CSF_healthy_filtered.txt",sep="\t",header=T,fill=T)

totalNum=dim(df_DERLE)[2]

cols=c("chrom","start",'end')
df_DERLE[,loci:=Reduce(function(...) paste(..., sep = "_"), .SD[, mget(cols)])]

print("here2")
dtnew <- df_DERLE[, lapply(.SD, as.integer), .SDcols=4:totalNum]
dt_use=as.data.frame(cbind(df_DERLE[,.(loci)],dtnew))
print("here3")
dtnew[, SUM := Reduce(`+`, .SD,init=1), .SDcols=1:dim(dtnew)[2]]

df_DERLE_sum=as.data.frame(cbind(df_DERLE[,.(loci)],dtnew[,.(SUM)]))
df_DERLE_sum=df_DERLE_sum[!is.na(df_DERLE_sum$SUM),]


dt_study=df_DERLE_sum[df_DERLE_sum$loci%in%RBP_rows$loci,]
dt_study$replacement_loci=rep("None",dim(dt_study)[1])
rm(RBP_rows)

rm(dtnew)
rm(df_DERLE)

print(head(df_DERLE_sum))
use_Df=df_DERLE_sum[!df_DERLE_sum$loci%in%dt_study$loci,]


for (lociB in unique(dt_study$loci)) {
  print(lociB)
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
  
  dt_study[dt_study$loci==lociB,3]=selected
  use_Df=use_Df[use_Df$loci!=selected,]
  
}

lociA=dt_study[,1]
lociA=lociA[!duplicated(lociA)]
lociB=dt_study[,3]

RBP_t=transpose(dt_use[dt_use$loci%in%c(lociA,lociB),],make.names = "loci")
rm(dt_use)
#Each column represents the correlation to all loci in that RBP by a random loci
output_cor=cor(RBP_t[,colnames(RBP_t)%in%lociB],
               RBP_t[,colnames(RBP_t)%in%lociB])

gz2=gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CSF/",RBP,"_corTable_subset.txt.gz",sep=""),"w")
write.table(output_cor,gz2,
            sep="\t",quote=F)
close(gz2)


#For testing with the DERLE table
# RBP_in=read.table("DGCR8_corComplete_subset.txt",sep="\t",header=T)
# RBP_in=RBP_in[RBP_in$column%in%lociA,]
# RBP_in=RBP_in[RBP_in$row%in%lociA,]
RBP_in_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/csf/",RBP,"_corComplete_subset.txt.gz",sep="")
RBP_in=as.data.frame(fread(RBP_in_files,sep="\t",header=T))
RBP_in[is.na(RBP_in)]=1

all_loci=unique(c(RBP_in$row,RBP_in$column))
print("here")
output=data.frame(locus="x",
                  p.value=1,
                  random_p.value=1,
		  d.value=0)

for (loci_num in 1:length(all_loci)) {
  loci=all_loci[loci_num]
  print(paste(loci,"p val loop"))
  
  rows_to_use=RBP_in[RBP_in$row==loci,]
  rows_to_use=rbind.data.frame(rows_to_use,RBP_in[RBP_in$column==loci,])
  
  KS.p.val=ks.test(rows_to_use$cor,rows_to_use$otherCor)$p.value
  KS.d.val=as.numeric(ks.test(rows_to_use$cor,rows_to_use$otherCor)$statistic)
  #output_cor[-loci_num,loci_num] because we want the same number of correlations between a random one that replaces the loci we're testing
  KS.p.val_otherLoci=ks.test(output_cor[-loci_num,loci_num],rows_to_use$otherCor)$p.value
  
  
  output=rbind.data.frame(output,
                          data.frame(locus=loci,
				     p.value=KS.p.val,
				     random_p.value=KS.p.val_otherLoci,
				     d.value=KS.d.val))
}

output=output[-1,]
gz1=gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CSF/",RBP,"_pValTable_subset.txt.gz",sep=""),"w")
write.table(output,gz1,
            sep="\t",quote=F,row.names = F)
close(gz1)
# KS.p.val_RBP=ks.test(as.numeric(output$p.value),as.numeric(output$random_p.value))$p.value
# plot(ecdf(output$p.value), xlim = range(c(output$p.value, output$random_p.value)), lty = "dashed",col="green")
# plot(ecdf(output$random_p.value), add = TRUE,col="red")



