args = commandArgs(trailingOnly=TRUE)

RBP_regions_files = list.files(path=paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/",args[1],sep=""),
                               pattern = "*_RBP_subset.txt.gz",full.names = T)

for (RBP in RBP_regions_files) {
  RBP_regions=read.table(RBP,sep="\t",header=T)
  RBP_regions=RBP_regions[!duplicated(RBP_regions[,1]),]
  rownames(RBP_regions)=RBP_regions[,1]
  RBP_regions=RBP_regions[,-1]
  
  table_in=RBP_regions
  coverage=as.numeric(args[2])
  table_in[table_in<coverage]=0
  table_in[table_in>coverage]=1
  table_in$total=rowSums(table_in)
  num_use=as.numeric(args[3])/100
  num=round(dim(table_in)[2]*num_use)
  #num=30
  RBP_regions=RBP_regions[table_in$total>num,]
  rm(table_in)
  if (dim(RBP_regions)[1]>0) {
   dat <- t(RBP_regions)
   cMat <- abs(cor(dat)) >= (1 - .Machine$double.eps^0.5)
   whichKeep <- which(rowSums(lower.tri(cMat) * cMat) == 0)
   test_out=RBP_regions[whichKeep,]
  

   name = gsub("_RBP_subset.txt.gz","",RBP)
   name = gsub(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/RBP_subsets/",args[1],"/",sep=""),"",name)
   print(name)
   gz1=gzfile(paste("/mnt/brlstor/Vol6_SP/exrna/readCoverage/EmilyTest/RBP_correlations/cover5_sample30/",args[1],"/",args[3],"perc_cov",args[2],"/",
                    name,"_regions.txt.gz",sep=""),"w")
   write.table(rownames(test_out),gz1,
               sep="\t",quote=F)
   close(gz1)
 }
  
}

