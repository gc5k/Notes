source("MDRFun.R")
loci=10
Mlocus=2
frq=runif(loci, min=0.05, max=0.5)
ld=array(0, dim=length(frq))
NTr=1000
Nt=500
g=GenerateGeno(frq, ld, NTr)
gt=GenerateGeno(frq, ld, Nt)
yTr=rnorm(NTr)#+g[,1]*0.5
yTr=scale(yTr)
yT=rnorm(Nt)#+gt[,1]*0.5
yT=scale(yT)
comb=t(combn(loci, Mlocus))
TRA=array(0, dim=nrow(comb))
StatTr=array(0, dim=nrow(comb))
BTR=array(0, dim=nrow(comb))
PTR=array(0, dim=nrow(comb))
RsqTR=array(0, dim=nrow(comb))
StatT=array(0, dim=nrow(comb))
TA=array(0, dim=nrow(comb))
BT=array(0, dim=nrow(comb))
PT=array(0, dim=nrow(comb))
RsqT=array(0, dim=nrow(comb))
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
  modTr=lm(yTr~xTr)
  BTR[i]=modTr$coefficients[2]
  StatTr[i] = summary(modTr)$fstatistic[1]
  PTR[i]=summary(modTr)$coefficients[2,4]
  TRA[i] = MDRScoreTrain(scheme, MG, gCode, yTr)
  RsqTR[i]=summary(modTr)$r.squared
  Loss=LossFun(gCode, gCodeT)
  print(Loss)
  xT=MDR(scheme, MG, gCodeT)
  modT=lm(yT~xT)
  BT[i]=modT$coefficients[2]
  StatT[i] = summary(modT)$fstatistic[1]
  PT[i]=summary(modT)$coefficients[2,4]
  TA[i] = MDRScoreTest(scheme, MG, gCodeT, yT)
  RsqT[i] = summary(modT)$r.squared
}
layout(matrix(1:6, 3, 2, byrow=F))
plot(BTR, TRA)
plot(BTR, -log10(PTR))
plot(TRA, -log10(PTR))
plot(BT, TA)
plot(BT, -log10(PT))
plot(TA, -log10(PT))
plot(density(PTR))
layout(matrix(1:2, 1, 2))
hist(PTR, breaks=25, xlab="p-values", main="Distribtuion of the p-value \nfor the training set")
hist(PT, breaks=25, xlab="p-values", main="Distribution of the p-value \nfor the test set")
pCut=0.05/nrow(comb)
layout(matrix(1:2, 1, 2))
plot(TRA, RsqTR, xlab="Training accuracy", ylab=expression(R^{2}), xlim=c(0.35, 0.65), ylim=c(0, 0.05), col=ifelse( PTR < pCut, "red", "blue"))
plot(TA, RsqT, xlab="Test accuracy", ylab=expression(R^{2}), xlim=c(0.35, 0.65), ylim=c(0, 0.05), col=ifelse(PT < pCut, "red", "blue"))
layout(matrix(1:2, 1, 2))
plot(TRA, BTR, xlab="Training accuracy", ylab=expression(beta), xlim=c(0.35, 0.65), ylim=c(-0.6, 0.6), col=ifelse( TA < 0.5, "red", "blue"))
abline(a=0, b=0, col="gray", lty=2)
abline(v=0.5, col="gray", lty=2)
plot(TA, BT, xlab="Testing accuracy", ylab=expression(beta), xlim=c(0.35, 0.65), ylim=c(-0.6, 0.6), col=ifelse( TA < 0.5, "red", "blue"))
idx=which(PT < pCut)
points(TA[idx], BT[idx], pch=15, col="red")
abline(a=0, b=0, col="gray", lty=2)
abline(v=0.5, col="gray", lty=2)
