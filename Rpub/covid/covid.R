library(pROC) #Roc
library(MASS) #stepAIC
library(psych) #pairs, pairs.panel
library(glmnet) #lasso
library(lmvar) #cv.lm, does not allow missing data
library(boot) #cv.glm
library(randomForest) #
library(mice) #md.pattern
library(leaps) #best subset
source("~/R/MyLib/shotgun.R")

#read data & clean off missing data
dat=read.csv("covid-data-pph.csv", header = T)
dat=dat[dat[,22]!=0.5,] #delete missing value
md.pattern(dat) #check missing value pattern
rmNA=which(apply(is.na(dat), 1, any)) #remove

if(length(rmNA) > 0) {
  dt=data.frame(dat[-rmNA,])
  print(paste("Remove", length(rmNA), "records due to missing."))
  print(paste(nrow(dt), "samples remained."))
} else {
  dt=dat
}

#feature selection using AUC
aucMat=matrix(0, ncol(dt)-1, 1)
rownames(aucMat)=colnames(dt)[1:nrow(aucMat)]
for(i in 1:nrow(aucMat)) {
  rocObj=roc(dt$diag, dt[,i], quiet = T)
  aucMat[i,1]=rocObj$auc
}
layout(matrix(c(1,1,2,3), byrow=T, 2, 2))
barplot(main="COVID ROC", aucMat[,1], beside = T, border = F)
abline(h=c(0.8, 0.7), col=c("blue", "yellow"), lty=2)
plot(rocObj, main="ROC using clinical-score")

respose=dt[,c("score", "diag")]
fm=c("gaussian", "binomial")

sReg=array(0, dim=c(length(fm), ncol(dt)-1, 7))
for(i in 1:dim(sReg)[1]) {
  for(j in 1:dim(sReg)[2]) {
    lm1=glm(respose[,i]~dt[,j], family=fm[i], data=dt)
    sReg[i,j,1:4]=summary(lm1)$coefficient[2,]
    sReg[i,j,5]=lm1$null.deviance
    sReg[i,j,6]=lm1$deviance
  }
}
logP=t(matrix(-log10(sReg[,-21,4]), 20, 2, byrow = T))
colnames(logP)=colnames(dt)[1:20]

layout(matrix(c(1,1,3,2,2,3), byrow=T, 2, 3))
barplot(main="COVID ROC", aucMat[,1], beside = T, border = F, col="green")
abline(h=c(0.8, 0.7), col=c("blue", "yellow"), lty=2)

par(las=2, mar=c(6,4,2,2))
barplot(logP,angle = 45, beside = T, horiz = F, border = F, ylab="-log10(p)", col=c("blue", "green"), main="Simple linear regression")
abline(h=c(8,10), col=c("red", "blue"), lty=2)
legend("topleft", legend = c("ClinicScore", "KitDiag"), pch=15, col=c("black", "green"), bty = 'n')

plot(rocObj, main="ROC using clinical-score")

##data confusion matrix
DiagCutoff=c(3, 3.5, 4)
DiagTab=array(0, dim=c(length(DiagCutoff), 2, 2))
for(i in 1:length(DiagCutoff)) {
  KitDiag=ifelse(dt$diag ==0, "Yin", "Yang")
  ScoreDiag=ifelse(dt$score < DiagCutoff[i], "Yin", "Yang")
  DiagTab[i,,]=table(KitDiag, ScoreDiag)
  pairs(datT[,c(21,22)])
  pairs.panels(datT[, c(21,22)])
}

#make different score
pCut=c(0.05, 0.01, 0.001)
aucM=matrix(0, length(pCut), 1)
for(i in 1:length(pCut)) {
  pV=which(sReg[2,,4]< pCut[i])
  ns=apply(dt[,pV],1, sum)
  rocObj=roc(dt[,22], ns)
  aucM[i,1]=rocObj$auc
  print(rocObj)
  lines(rocObj, type="b", pch=21, col=i, by="grey")
}

##cross-validaton for lm, cv.lm in lmvar package
###!no missing data is allowed for cv.lm
diaMod_LM=lm(diag~as.matrix(dt)[,pV], x=TRUE, y=TRUE, data=dt)
cvLM=cv.lm(diaMod_LM, k = 5, seed=1000, max_cores = 2)

TrSet=sort(sample(nrow(dt), ceiling(0.9*nrow(dt))))
yTr=dt$diag[TrSet]
dtTr=as.matrix(dt[TrSet, c(2,3,8,9,17,20)])
dtTs=data.frame(as.matrix(dt[-TrSet, c(2,3,8,9,17,20)]))
LM_Tr=lm(yTr~x2+x3+x8+x9+x17+x20, data=as.data.frame(dtTr))
LM_Tt=predict(LM_Tr, dtTs) #it is very important to match the column names for training and test model

##cross-validation for glm in boot package

Ksize=5
KF=sample(KFold(Ksize, nrow(dt)), nrow(dt))
cvMod=array(0, dim=c(Ksize, 4, 4))
for(i in 1:Ksize) {
  TrSet=which(KF!=i)
  dtTrGLM=data.frame(dt[TrSet, c("Home_Work", "CT", "diag")])
  dtTestGLM=dt[-TrSet, c("Home_Work", "CT", "diag")]

  diaMod_GLM0=glm(diag~Home_Work+CT, data=dtTrGLM)
  cvMod[i,1:3,]=summary(diaMod_GLM0)$coefficient
  preGLM=predict.glm(diaMod_GLM0, dtTestGLM)
  roc_obj=roc(dtTestGLM$diag, preGLM, quiet = T)
  cvMod[i,4,1]=roc_obj$auc
#  cvGLM=cv.glm(dtTrGLM, diaMod_GLM0, K=5)
}
barplot(cvMod[,4,1], border = F, xlab="Cross-validation", ylab="ROC")
##score cutoff

yRes=dt[,22]
dtX=dt[,c(-21, -22)]
mod=glm(yRes~., family=binomial, data=dtX)
fitStep=stepAIC(mod, direction="both", trace = TRUE, scope = list(upper=~., lower=~1))
fitBack=stepAIC(mod, direction="backward", trace = TRUE, scope = list(upper=~., lower=~1))
fitFor=stepAIC(mod, direction="forward", trace = FALSE, scope = list(upper=~., lower=~1))

dtSelX=dt[,c(names(fitBack$coefficients)[-1], "diag")]
diaMod_GLM0=glm(diag~x1+x2+x3+x6+x7+x8+x9+x16+x17+x18+x19+x20, family=binomial, data=dtSelX)
cvGLM=cv.glm(dtSelX, diaMod_GLM0, K=5)

##lasso library(glmnet)
#glmnet.mod=glmnet(as.matrix(dt[,-ncol(dt)]), dt$diag, family="binomial")
dtm=as.matrix(dt[,-ncol(dt)])
dtmS=apply(dtm, 2, scale)
glmnet.mod1=glmnet(dtmS, dt$diag, family="binomial")
plot(glmnet.mod1)
cvglm=cv.glmnet(dtmS, dt$diag, family="binomial")
plot(cvglm)
