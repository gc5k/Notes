fam1=read.table("IBDBatch1_rev.fam", as.is=T)
fam2=read.table("IBDBatch2_rev.fam", as.is=T)
phe=read.table("IBD_UC_nodup.phe", as.is=T, header=T)
phe$Batch=-9

for(i in 1:nrow(fam1))
{
  idx=which(phe[,2] == fam1[i,2])
  if(length(idx) > 0)
  {
    phe$Batch[idx] = 1
  }
}

for(i in 1:nrow(fam2))
{
  idx=which(phe[,2] == fam2[i,2])
  if(length(idx) > 0)
  {
    phe$Batch[idx] = 2
  }
}

write.table(phe, "IBD_UC.phe", row.names=F, col.names=T, quote=F)
