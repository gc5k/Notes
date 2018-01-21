args=commandArgs(TRUE)
iBKtrtF=args[1]
covF=args[2]

iBKtrt=read.table(iBKtrtF, as.is=T)
cov=read.table(covF, as.is=T, header=T)

iMat=matrix(0, nrow=nrow(iBKtrt), ncol=ncol(iBKtrt))
iMat[,1]=iBKtrt[,1]
iMat[,2]=iBKtrt[,2]
iMat[,3]=iBKtrt[,3]

for(i in 4:8)
{
  idxTr=which(iBKtrt[,i]==1)
  yTr=as.numeric(iBKtrt[idxTr,3]) - 1
  xTr=as.matrix(cov[idxTr,3:22])
  modTr=glm(yTr~xTr, family=binomial)
  iMat[idxTr,i]=yTr - modTr$fitted.values
  iMat[-idxTr,i]=-9
}

AiBKtrtF=sub("trt$", "trtA", iBKtrtF)
write.table(iMat, AiBKtrtF, row.names=F, col.names=F, quote=F)
