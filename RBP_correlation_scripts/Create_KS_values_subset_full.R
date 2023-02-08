args = commandArgs(trailingOnly=TRUE)
RBP = args[1]

library(data.table)


RBP_regions_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/CellLines/",RBP,"_RBP_subset.txt",sep="")
print(RBP_regions_files)
RBP_rows=fread(RBP_regions_files,sep="\t",header=T)
# RBP_rows=fread("DGCR8_RBP_subset.txt",sep="\t",header=T)

# df_DERLE = fread("EXR-DERLE1Dj6liR-AN_filtered.txt",sep="\t",header=T,fill=T)
df_DERLE = fread("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/biosample_RBP_out_unzipped/cellLines_all_summed_filtered.txt",sep="\t",header=T,fill=T)
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

rm(dtnew)
rm(df_DERLE)

print(head(df_DERLE_sum))
use_Df=df_DERLE_sum[!df_DERLE_sum$loci%in%dt_study$loci,]


for (lociB in dt_study$loci[!duplicated(dt_study$loci)]) {
  print(lociB)
  coverage = dt_study[dt_study$loci==lociB,2]
 
  possible=use_Df[use_Df$SUM==coverage,1]
  
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
lociB=dt_study[,3]

#print(length(dt_use$loci%in%c(lociA,lociB)))
print(dt_use$loci[!dt_use$loci%in%c(lociA,lociB)])

RBP_t=transpose(dt_use[dt_use$loci%in%c(lociA,lociB),],make.names = "loci")
dim(RBP_t)
#Each column represents the correlation to all loci in that RBP by a random loci
output_cor=cor(RBP_t[,colnames(RBP_t)%in%lociA],
               RBP_t[,colnames(RBP_t)%in%lociB])
print(length(lociA%in%colnames(RBP_t)))
print(length(lociB%in%colnames(RBP_t)))
print(dim(output_cor))
print(dim(cor(RBP_t[,colnames(RBP_t)%in%lociA],
               RBP_t[,colnames(RBP_t)%in%lociB])))

write.table(output_cor,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/",RBP,"_corTable_subset_fixed.txt",sep=""),
            sep="\t",quote=F)


#output_cor=as.data.frame(fread(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/",RBP,"_corTable_subset.txt",sep=""),sep="\t"))
#rownames(output_cor)=output_cor[,1]
#output_cor=output_cor[,-1]

RBP_in_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/CellLines/",RBP,"_corComplete_subset.txt",sep="")
RBP_in=as.data.frame(fread(RBP_in_files,sep="\t",header=T))
RBP_in[is.na(RBP_in)]=1

#all_loci=unique(c(RBP_in$row,RBP_in$column))

all_loci=read.table(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/cover5_sample30/CellLines/",RBP,"_regions.txt",sep=""))[,1]

if (length(all_loci)>2) {

 RBP_in=RBP_in[RBP_in$row%in%all_loci,]
 RBP_in=RBP_in[RBP_in$column%in%all_loci,]

 col_take=which(rownames(output_cor) %in% all_loci)
 print(max(col_take))
 output_cor=output_cor[all_loci,col_take]
 
 output=data.frame(locus="x",
                   p.value=1,
                   random_p.value=1,
    		  d.value=0)

 for (loci_num in 1:length(all_loci)) {
   loci=all_loci[loci_num]
   print(paste(loci,"p val loop"))
  
   rows_to_use=RBP_in[RBP_in$row==loci,]
   rows_to_use=rbind.data.frame(rows_to_use,RBP_in[RBP_in$column==loci,])
  
   #vISUALIZE DENSITY PLOTS
   # test_density=cbind.data.frame(c(rep("within",dim(rows_to_use)[1]),rep("random",dim(rows_to_use)[1])),
   #                               c(rows_to_use$cor,rows_to_use$otherCor))
   # sm.density.compare(test_density[,2],
   #                    as.factor(test_density[,1]),
   #                    xlab="Correlations")
   # colfill<-c(2:(2+length(levels(as.factor(test_density[,1])))))
   # legend(x = 0.5,y=2, levels(as.factor(test_density[,1])), fill=colfill)
   
   KS.p.val=ks.test(rows_to_use$cor,rows_to_use$otherCor)$p.value
   KS.d.val=as.numeric(ks.test(rows_to_use$cor,rows_to_use$otherCor)$statistic)
   #output_cor[-loci_num,loci_num] because we want the same number of correlations between a random one that replaces the loci we're testing
   KS.p.val_otherLoci=ks.test(output_cor[-loci_num,loci_num],rows_to_use$otherCor)$p.value
   
   #vISUALIZE KS PLOTS
   # plot(ecdf(rows_to_use$cor), xlim = range(c(rows_to_use$cor, rows_to_use$otherCor)), lty = "dashed",col="green")
   # plot(ecdf(rows_to_use$otherCor), add = TRUE,col="red")
   
   output=rbind.data.frame(output,
                           data.frame(locus=loci,
 				     p.value=KS.p.val,
 				     random_p.value=KS.p.val_otherLoci,
 				     d.value=KS.d.val))
 }
 
 output=output[-1,]
 write.table(output,paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/CellLines/cov5_sample30/",RBP,"_pValTable_cov30.txt",sep=""),
             sep="\t",quote=F,row.names = F)
}

# KS.p.val_RBP=ks.test(as.numeric(output$p.value),as.numeric(output$random_p.value))$p.value
# plot(ecdf(output$p.value), xlim = range(c(output$p.value, output$random_p.value)), lty = "dashed",col="green")
# plot(ecdf(output$random_p.value), add = TRUE,col="red")



