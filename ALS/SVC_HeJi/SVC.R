dat=read.csv("SVC_2020_0626.csv", as.is = T, header = T)

Smat=matrix(0, ncol(dat)-2, 4)
for(i in 3:ncol(dat)) {
  mod=glm(dat[,2]~dat[,i], family = "binomial")
  Smat[i-2,]=summary(mod)$coefficients[2,]
}
rownames(Smat)=colnames(dat)[3:ncol(dat)]

library(pROC)
aucMat=matrix(0, ncol(dat)-2, 1)
rownames(aucMat)=colnames(dat)[3:ncol(dat)]
for(i in 3:ncol(dat)) {
  rocObj=roc(dat[,2], dat[,i], quiet = T)
  aucMat[i-2,1]=rocObj$auc
}

layout(matrix(c(1,2,3,4), 2, 2, byrow = F))
par(las=2, cex=0.5, mai=c(1,0.5,0.5,0.5))
pcut=-log10(0.05/(ncol(dat)-2))
barplot(main="case-control",-log10(Smat[,4]), col=ifelse(-log10(Smat[,4])>pcut, "cyan", "grey"))
abline(h=-log10(0.05/ncol(dat)), col="red", lty=2)

barplot(main="case-control",aucMat, beside = T, col=ifelse(-log10(Smat[,4])>pcut, "cyan", "grey"))
plot(-log10(Smat[,4]), pch=16,aucMat,col=ifelse(-log10(Smat[,4])>pcut, "cyan", "grey"))
abline(v=pcut, col="red", lty=2)

###########
SVCmat=matrix(1, ncol(dat)-2, 4)
for(i in 3:ncol(dat)) {
  if(i!=68) {
    mod=lm(dat[,68]~dat[,i])
    if(nrow(summary(mod)$coefficients)>1) {
      SVCmat[i-2,]=summary(mod)$coefficients[2,]
    }
  }
}
rownames(SVCmat)=colnames(dat)[3:ncol(dat)]

par(las=2, cex=0.5, mai=c(1.5,0.5,0.5,0.5))
pcut=-log10(0.05/(ncol(dat)-2))
barplot(-log10(Smat[,4]), col=ifelse(-log10(Smat[,4])>pcut, "cyan", "grey"))
abline(h=-log10(0.05/ncol(dat)), col="red")
