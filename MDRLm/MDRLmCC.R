source("MDRFun.R")
loci=20
Mlocus=2

frq=runif(loci, min=0.05, max=0.5)
ld=array(0, dim=length(frq))
NTr=1000
Nt=500
g=GenerateGeno(frq, ld, NTr)
gt=GenerateGeno(frq, ld, Nt)
yTr=rbinom(NTr, 1, 0.5)
yT=rbinom(Nt, 1, 0.5)
comb=t(combn(loci, Mlocus))

StatTr=array(0, dim=nrow(comb))
TRA=array(0, dim=nrow(comb))

ngBTR=array(0, dim=nrow(comb))
ngPTR=array(0, dim=nrow(comb))
ngRsqTR=array(0, dim=nrow(comb))
ngRsqTR=array(0, dim=nrow(comb))


StatT=array(0, dim=nrow(comb))
TA=array(0, dim=nrow(comb))

ngBT=array(0, dim=nrow(comb))
ngPT=array(0, dim=nrow(comb))
ngRsqT=array(0, dim=nrow(comb))
ngRsqNgT=array(0, dim=nrow(comb))

for(i in 1:nrow(comb))
{
  print(comb[i,])
  gCode=MultiCode(g, comb[i,])
  gCodeT=MultiCode(gt, comb[i,])
  
  dat=GetScheme(gCode, yTr)
  scheme=dat[[1]]
  MG=dat[[2]]
  print(scheme)
  
  xTr=MDR(scheme, MG, gCode)
  TRA[i] = MDRScoreTrain(scheme, MG, gCode, yTr)
  #TR linear
  gmodTrFull=glm(yTr~xTr, family='binomial')
  ngBTR[i]=summary(gmodTrFull)$coefficients[2,1]
  ngPTR[i]=summary(gmodTrFull)$coefficients[2,4]
  gmodTrNull=glm(yTr~1, family='binomial')
  LLKnullTr=logLik(gmodTrNull)
  LLKfullTr=logLik(gmodTrFull)
  R_cs=1-10^((LLKnullTr-LLKfullTr)*2/NTr)
  ngRsqTR[i]=R_cs/(1-10^(LLKnullTr*2/NTr))
  
  Loss=LossFun(gCode, gCodeT)
  print(Loss)
  
  xT=MDR(scheme, MG, gCodeT)
  TA[i] = MDRScoreTest(scheme, MG, gCodeT, yT)
  gmodTFull=glm(yT~xT, family='binomial')
  ngBT[i]=summary(gmodTFull)$coefficients[2,1]
  ngPT[i]=summary(gmodTFull)$coefficients[2,4]
  gmodTNull=glm(yT~1, family='binomial')
  LLKnullT=logLik(gmodTNull)
  LLKfullT=logLik(gmodTFull)
  R_cs=1-10^((LLKnullT-LLKfullT)*2/Nt)
  ngRsqT[i]=R_cs/(1-10^(LLKnullT*2/Nt))
}
layout(matrix(1:6, 3, 2, byrow=F))
plot(ngBTR, TRA, xlab="Regression coefficient", ylab="Training accuracy")
plot(ngBTR, -log10(ngPTR), xlab="Regression coefficient", ylab="p-values")
plot(TRA, -log10(ngPTR), xlab="Training accuracy", ylab="p-values")

plot(ngBT, TA, xlab="Regression coefficient", ylab="Testing accuracy")
plot(ngBT, -log10(ngPT), xlab="Regression coefficient", ylab="p-values")
plot(TA, -log10(ngPT), xlab="Testing accuracy", ylab="p-values")

#pvalue
layout(matrix(1:2, 1, 2))
hist(ngPTR, breaks=25, xlab="p-values", main="Distribtuion of the p-value \nfor the training set")
hist(ngPT, breaks=25, xlab="p-values", main="Distribution of the p-value \nfor the test set")

pCut=0.05/nrow(comb)
layout(matrix(1:2, 1, 2))
plot(TRA, ngRsqTR, xlab="Training accuracy", ylab=expression(paste(R[N]^{2})), xlim=c(0.35, 0.65), ylim=c(0, 0.07), col=ifelse( ngPTR < pCut, "red", "blue"))
plot(TA, ngRsqT, xlab="Test accuracy", ylab=expression(paste(R[N]^{2})), xlim=c(0.35, 0.65), ylim=c(0, 0.07), col=ifelse( ngPT < pCut, "red", "blue"))

#plot(TRA, RsqTR, xlab="Training accuracy", ylab=expression(R[l]^{2}), xlim=c(0.35, 0.65), ylim=c(0, 0.05), col=ifelse( PTR < pCut, "red", "blue"))
#plot(TA, RsqT, xlab="Test accuracy", ylab=expression(R[l]^{2}), xlim=c(0.35, 0.65), ylim=c(0, 0.05), col=ifelse(PT < pCut, "red", "blue"))

plot(ngRsqTR, ngRsqT)
#plot(RsqTR, RsqT)


layout(matrix(1:2, 1, 2))
plot(TRA, ngBTR, xlab="Training accuracy", ylab=expression(beta), xlim=c(0.35, 0.65), ylim=c(-1.1*max(ngBTR), 1.1*max(ngBTR)), col=ifelse( TA < 0.5, "red", "blue"))
abline(a=0, b=0, col="gray", lty=2)
abline(v=0.5, col="gray", lty=2)
plot(TA, ngBT, xlab="Testing accuracy", ylab=expression(beta), xlim=c(0.35, 0.65), ylim=c(-1.1*max(ngBT), 1.1*max(ngBT)), col=ifelse( TA < 0.5, "red", "blue"))
idx=which(ngPT < pCut)
if(length(idx) > 0)
{
  points(TA[idx], ngBT[idx], pch=15, col="black")  
}
abline(a=0, b=0, col="gray", lty=2)
abline(v=0.5, col="gray", lty=2)
