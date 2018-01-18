fam=read.table("IBDMerge_tmp.fam", as.is=T)
phe=read.table("IBD_UC.phe", as.is=T, header=T)

for(i in 1:nrow(fam))
{
  fam[i,1] = fam[i,2]
  idx = which(phe[,2] == fam[i,2])
  if(length(idx) > 0)
  {
    if(length(idx) > 1)
    {
      print(idx)
    }
    fam[i,5] = phe$Gender[idx]
    fam[i,6] = phe$CC[idx]
  }
}
write.table(fam, "IBDMerge_tmp.fam", row.names=F, col.names=F, quote=F)

unknown=which(phe$Batch == -9)
write.table(phe[unknown, c(1,2)], "UnknownInd.txt", row.names=F, col.names=F, quote=F)
