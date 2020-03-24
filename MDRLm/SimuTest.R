source("MDRFun.R")

###quantitative traits
M=10 #SNP feature
N=100 #sample size
frq=runif(M, 0.1, 0.9) #frequency
locID=c(1,3) #functional loci
gFun=c(12, 1, 13) #
intcp=0 #intercept of the model
gEff=1 #effect of the functional loci
cBeta=c(2) #
rsq= 0.2 #Rsq of the model

dat=SimuQt(frq, N, locIdx = locID, gFun=gFun, gbeta = gEff, cBeta = cBeta, rsq)
mod=lm(y~X-1, data = dat)
res=mod$residuals

###case-control
csN=100
ctrlN=100
datCC=SimuCC(frq, csN, ctrlN, locIdx=locID, gFun = gFun, intercept = -4.5, cBeta = cBeta)
modL=glm(y~X-1, data = datCC, family="binomial")
resCC=datCC$y-modL$fitted.values
modL_ft.value=exp(datCC$X%*%modL$coefficients)/(1+exp(datCC$X%*%modL$coefficients))

modL2=glm(y~1, data = datCC, family="binomial")
resCC2=datCC$y-modL2$fitted.values
modL2_ft.value=exp(datCC$X%*%modL2$coefficients)/(1+exp(datCC$X%*%modL2$coefficients))
