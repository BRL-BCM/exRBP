args = commandArgs(trailingOnly=TRUE)
RBP = args[1]
SET=args[2]
BIO=args[3]
bio=args[4]

library(data.table)


output_cor=as.data.frame(fread(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",BIO,"/",RBP,"_corTable_subset.txt.gz",sep=""),sep="\t"))
rownames(output_cor)=output_cor[,1]
output_cor=output_cor[,-1]
RBP_in_files = paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/CompletedCors/",bio,"/",RBP,"_corComplete_subset.txt.gz",sep="")
RBP_in=as.data.frame(fread(RBP_in_files,sep="\t",header=T))
RBP_in[is.na(RBP_in)]=1

#print(head(RBP_in))
#all_loci=unique(c(RBP_in$row,RBP_in$column))

all_loci=read.table(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/cover5_sample30/",BIO,"/",SET,"/",RBP,"_regions.txt.gz",sep=""))[,1]
#print(all_loci)

if (length(all_loci)>4) {
 #print("in loop")
 RBP_in=RBP_in[RBP_in$row%in%all_loci,]
 RBP_in=RBP_in[RBP_in$column%in%all_loci,]

 col_take=which(rownames(output_cor) %in% all_loci)
 
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
 #print(output)
 gz1=gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/Cor_tables/",BIO,"/cov5_sample30/",SET,"/",RBP,"_pValTable_cov30.txt.gz",sep=""),"w")
 write.table(output,gz1,
             sep="\t",quote=F,row.names = F)
 close(gz1)
}

#print("out loop")

# KS.p.val_RBP=ks.test(as.numeric(output$p.value),as.numeric(output$random_p.value))$p.value
# plot(ecdf(output$p.value), xlim = range(c(output$p.value, output$random_p.value)), lty = "dashed",col="green")
# plot(ecdf(output$random_p.value), add = TRUE,col="red")



