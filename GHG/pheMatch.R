args=commandArgs(TRUE)
sc=as.numeric(args[1])
real=args[2]
phe1=args[3]
phe2=args[4]
out1=args[5]
out2=args[6]

sc=0.9
real="ichip2cd.real"
phe1= "plink.cd.txt" 
phe2= "cd0.dat3"
out1= "cd_ichip_keep.phe"
out2= "cd_gwas_keep.phe"

Match=read.table(real, as.is=T, header=T)
idx=which(Match$Score > sc)
Match=Match[idx,]
Match$order=seq(1, nrow(Match))

L1=Match[,c(1,2, ncol(Match))]
L1$UID=paste0(L1[,1], ".", L1[,2])
PHE1=read.table(phe1, as.is=T)
PHE1$UID=paste0(PHE1[,1], ".", PHE1[,2])

m1=merge(L1, PHE1, by.x="UID", by.y="UID", sort=F)

L2=Match[,c(3,4, ncol(Match))]
L2$UID=paste0(L2[,1], ".", L2[,2])
PHE2=read.table(phe2, as.is=T)
PHE2$UID=paste0(PHE2[,1], ".", PHE2[,2])
m2=merge(L2, PHE2, by.x="UID", by.y="UID", sort=F)

commID=intersect(m1$order, m2$order)
filter1=which(m1$order %in% commID)
filter2=which(m2$order %in% commID)
m1=m1[filter1,]
m2=m2[filter2,]

idx=which(m1[,ncol(m1)] == 0 | m2[,ncol(m2)] == 0)
if(length(idx) > 0)
{
  m1=m1[-idx,]
  m2=m2[-idx,]
}

idx1=which(m1[,ncol(m1)] == m2[,ncol(m2)])

write.table(m1[idx1,c(2,3,ncol(m1))], out1, row.names=F, col.names=F, quote=F)
write.table(m2[idx1,c(2,3,ncol(m2))], out2, row.names=F, col.names=F, quote=F)

write.table(m1[-idx1,c(2,3,ncol(m1))], paste0(out1, ".dis"), row.names=F, col.names=F, quote=F)
write.table(m2[-idx1,c(2,3,ncol(m1))], paste0(out2, ".dis"), row.names=F, col.names=F, quote=F)
