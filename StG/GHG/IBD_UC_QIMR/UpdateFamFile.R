phe=read.table("IBD_UC.phe", as.is=T, header=T)
fam=read.table("IBD_UC_tmp2.fam", as.is=T)

for(i in 1:nrow(fam))
{
  idx=which(phe$FID == fam[i,2])
  if(length(idx) > 0)
  {
    fam[i, 5] = phe$Gender[idx]
    fam[i, 6] = phe$CC[idx]
  }
}
write.table(fam, "IBD_UC_tmp2.fam", row.names=F, col.names=F, quote=F)
