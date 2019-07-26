source("~/R/MyLib/shotgun.R")
M=2000 #marker
N=1000 #snp
h2=0.5 #snp-heritability
#fq=c(0.5, 0.3)

fq=runif(M, 0.1, 0.5)
#ld=c(0,0)
ld=rep(0, M)
ld=runif(M, -0.9, 0.9)

G=GenerateGenoDprime(fq, ld, N)
FRQ=1-colMeans(G)/2


REP=50
EXV=matrix(0, REP, 2)
JexpV=matrix(0, REP, 2)
expV=matrix(0, REP, 2)
for(i in 1:REP) {
#  bg=G %*% rnorm(M, 0, 1) #random
  bg=G %*% rnorm(M, 0, 1/sqrt(FRQ*(1-FRQ))) #inverse to 2pq
#  bg=G %*% rnorm(M, 0, 1/sqrt(FRQ*(1-FRQ))) #proportional to 2pq
#  bg=G %*% rnorm(M, 0, sqrt(FRQ*(1-FRQ))) #proportional to 2pq

  y=bg+rnorm(N, 0, sqrt(var(bg)/h2*(1-h2)))
  y=scale(y)  #standardization
  Y=matrix(y, nrow=length(y), ncol=1)%*%matrix(y, nrow=1, ncol=length(y))
  YY=matrix(Y[row(Y) < col(Y)], nrow=N*(N-1)/2, 1)
  
  sG=apply(G, 2, scale) #standardization snp
  
  idxB=which(FRQ<0.25)
  
  sGRM=sG %*% t(sG)/M
  sGRM1=sG[,idxB] %*% t(sG[,idxB])/length(idxB)
  sGRM2=sG[,-idxB] %*% t(sG[,-idxB])/(M-length(idxB))
  
  BG=sGRM[row(sGRM) < col(sGRM)]
  BG1=sGRM1[row(sGRM1) < col(sGRM1)]
  BG2=sGRM2[row(sGRM2) < col(sGRM2)]
  

  
  X2=matrix(0, nrow=length(BG1), 3)
  X2[,1]=1
  X2[,2]=BG1
  X2[,3]=BG2
  XtX2=t(X2)%*%X2
  
  beta=solve(XtX2)%*%t(X2)%*%YY
  
  me1=1/var(BG1)
  me2=1/var(BG2)
  me12=1/cov(BG1, BG2)
  
  NN=N*(N-1)/2
  detA=determinant(XtX2)
  int_X=solve(XtX2)
  
  XtX_xT=int_X %*%t(X2)
  xy=t(X2) %*% YY
  bbS=int_X %*% xy
  
  
#############
  m1=lm(YY~BG1)
  m2=lm(YY~BG2)
  expV[i, 1]=m1$coefficients[2]
  expV[i, 2]=m2$coefficients[2]
  JexpV[i,]=bbS[-1,1]
  EXV[i, 1]=cov(BG, BG1)/var(BG1)*h2
  EXV[i, 2]=cov(BG, BG2)/var(BG2)*h2
}
barplot(expV[,1])
abline(h=mean(EXV[,1]))
barplot(expV[,2])
abline(h=mean(EXV[,2]))
colMeans(EXV)
colMeans(expV)
