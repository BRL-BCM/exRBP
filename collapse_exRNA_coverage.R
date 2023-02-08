#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#print(args[1])
cov=as.numeric(args[1])
files=list.files(".",pattern = "*.bed")
df=read.table(files[1],sep="\t",header = T)
#df_bed=df[,1:3]
print(head(df[,1:5]))

row_set=c()
total=c()
for (n in 1:dim(df)[1]) {
  res=sum(5 <= df[n,-c(1:3)])
  if (res>0) {
    row_set=c(row_set,n)  
    total=c(total,res)
  }
}
print("out")

#df_rm=df[,-c(1:3)]
#rm(df)
#print(head(df_rm))
#df_rm[df_rm>cov]=cov
#df_rm[df_rm<cov]=0

#sets=apply(df_rm[,-c(1:3)], 1, function(r) any(r %in% c(cov)))
#sets=rowSums(df[,-c(1:3)]*(df[,-c(1:3)]>=cov))
#sets=rowSums(df_rm)
#print(head(sets))

#out=cbind.data.frame(df_bed,sets)
out=df[row_set,c(1:3)]
rm(df)
out$sets=total
#print(head(out))
#out=out[out[,4]>1,]
print(head(out))
write.table(out,"out.txt",sep="\t",quote=F,row.names=F)
