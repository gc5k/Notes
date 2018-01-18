args=commandArgs(TRUE)
iBKtrtF=args[1]
gBKtrtF=args[2]
itrtF=args[3]
gtrtF=args[4]
covF=args[5]

#iBKtrtF="OBKcd1iChip.trt"
#gBKtrtF="OBKcd1gChip.trt"
#itrtF="Ocd1iChip.trt"
#gtrtF="Ocd1gChip.trt"
#covF="Ocd1gChip.profile"

iBKtrt=read.table(iBKtrtF, as.is=T)
gBKtrt=read.table(gBKtrtF, as.is=T)
cov=read.table(covF, as.is=T, header=T)

iMat=matrix(0, nrow=nrow(iBKtrt), ncol=ncol(iBKtrt))
gMat=matrix(0, nrow=nrow(gBKtrt), ncol=ncol(gBKtrt))
iMat[,1]=iBKtrt[,1]
iMat[,2]=iBKtrt[,2]
iMat[,3]=iBKtrt[,3]

gMat[,1]=gBKtrt[,1]
gMat[,2]=gBKtrt[,2]
gMat[,3]=gBKtrt[,3]

for(i in 4:8)
{
  idxTr=which(iBKtrt[,i]==1)

  yTr=as.numeric(iBKtrt[idxTr,3]) - 1
  xTr=as.matrix(cov[idxTr,3:22])
  modTr=glm(yTr~xTr, family=binomial)
  iMat[idxTr,i]=yTr - modTr$fitted.values

  yT=as.numeric(iBKtrt[-idxTr,3]) - 1
  xT=as.matrix(cov[-idxTr,3:22])
  modT=glm(yT~xT, family=binomial)
#  iMat[-idxTr,i]=yT - modT$fitted.values
  iMat[-idxTr,i]=-9
  gMat[,i]=iMat[,i]
}

AiBKtrtF=sub("trt$", "trtA", iBKtrtF)
AgBKtrtF=sub("trt$", "trtA", gBKtrtF)
AitrtF=sub("trt$", "trtA", itrtF)
AgtrtF=sub("trt$", "trtA", gtrtF)

write.table(iMat, AiBKtrtF, row.names=F, col.names=F, quote=F)
write.table(gMat, AgBKtrtF, row.names=F, col.names=F, quote=F)
write.table(iMat, AitrtF, row.names=F, col.names=F, quote=F)
write.table(gMat, AgtrtF, row.names=F, col.names=F, quote=F)
