set.seed(2014)
args=commandArgs(TRUE)
iBKfamF=args[1]
gBKfamF=args[2]
ifamF=args[3]
gfamF=args[4]

#common
#iBKfamF="OBKcd1iChip.fam"
#gBKfamF="OBKcd1gChip.fam"
#ifamF="Ocd1iChip.fam"
#gfamF="Ocd1gChip.fam"

iBKfam=read.table(iBKfamF, as.is=T)
gBKfam=read.table(gBKfamF, as.is=T)

TrT=sample(nrow(iBKfam), nrow(iBKfam))
iTmat=matrix(1, nrow=nrow(iBKfam), ncol=8)
iTmat[,1] = iBKfam[,1]
iTmat[,2] = iBKfam[,2]
iTmat[,3] = iBKfam[,6]

gTmat=matrix(1, nrow=nrow(gBKfam), ncol=8)
gTmat[,1] = gBKfam[,1]
gTmat[,2] = gBKfam[,2]
gTmat[,3] = gBKfam[,6]

len=ceiling(nrow(iBKfam)/5)
for(i in 1:4)
{
  iTmat[TrT[(1+(i-1)*len):(i*len)],3+i]="0"
  gTmat[TrT[(1+(i-1)*len):(i*len)],3+i]="0"
}
iTmat[TrT[(1+i*len):length(TrT)],8]="0"
gTmat[TrT[(1+i*len):length(TrT)],8]="0"

iBKfamOutF=sub("fam$", "trt", iBKfamF)
gBKfamOutF=sub("fam$", "trt", gBKfamF)
ifamOutF=sub("fam$", "trt", ifamF)
gfamOutF=sub("fam$", "trt", gfamF)

write.table(iTmat, iBKfamOutF, row.names=F, col.names=F, quote=F)
write.table(gTmat, gBKfamOutF, row.names=F, col.names=F, quote=F)
write.table(iTmat, ifamOutF, row.names=F, col.names=F, quote=F)
write.table(gTmat, gfamOutF, row.names=F, col.names=F, quote=F)

