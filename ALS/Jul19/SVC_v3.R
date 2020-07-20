dat=read.csv("20200718all_cgb.csv", as.is = T, header = T)

par(las=2, font=5, bty='l', ps=6)
boxplot(dat[,-2], cex=0.5, pch=16, col=runif(nrow(dat), 1, nrow(dat)))

##basic accessment of the data
###fivenum, mean, sd, missing, isNumeric
cutoff=5 #outlier pickup

SUMmat=matrix(0, ncol(dat), 9)
colnames(SUMmat)=c("minimum", "lower-hinge", "median", "upper-hinge", "maximum", "mean", "sd", "missingRate", "isNumeric")
rownames(SUMmat)=colnames(dat)
for(i in 1:nrow(SUMmat)) {
  if(is.numeric(dat[,i])) {
    SUMmat[i,1:5]=fivenum(dat[,i])
    SUMmat[i,6]=mean(dat[,i], na.rm = T)
    SUMmat[i,7]=sd(dat[,i], na.rm = T)
    SUMmat[i,8]=length(which(is.na(dat[,i])))/nrow(dat)
    SUMmat[i,9]=T
    
    idx=which(dat[,i] > SUMmat[i,6]+cutoff*SUMmat[i,7])
    if(length(idx)>0) {
      print(paste0("", colnames(dat)[i], ", sample id: ",dat[idx,2]))
      print(paste0("outlier value: ",  dat[idx, i]))

    }
  } else {
    SUMmat[i,1:8]=NA
    SUMmat[i,9]=F
  }
}


Smat=matrix(0, ncol(dat)-2, 4)

for(i in 3:ncol(dat)) {
  mod=glm(dat[,1]~dat[,i], family = "binomial")
  Smat[i-2,]=summary(mod)$coefficients[2,]
}
rownames(Smat)=colnames(dat)[3:ncol(dat)]

library(pROC)
aucMat=matrix(2, ncol(dat)-2, 1)
rownames(aucMat)=colnames(dat)[3:ncol(dat)]
for(i in 3:ncol(dat)) {
  rocObj=roc(dat[,1], dat[,i], quiet = T)
  aucMat[i-2,1]=rocObj$auc
}

layout(matrix(c(1,2,3,4), 2, 2, byrow = F))
par(las=2, cex=0.5, mai=c(1,0.5,0.5,0.5))
pcut=-log10(0.05/(ncol(dat)-2))
barplot(main="case-control",-log10(Smat[,4]), col=ifelse(-log10(Smat[,4])>pcut, "cyan", "grey"))
abline(h=-log10(0.05/ncol(dat)), col="red", lty=2)

barplot(main="case-control", aucMat, beside = T, col=ifelse(-log10(Smat[,4])>pcut, "cyan", "grey"))
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
