library(Rcpp)
sourceCpp("~/git/Notes/R/RLib/Shotgun.cpp")
source("~/git/Notes/R/RLib/shotgun.R")
M=10000
N=500
frq=runif(M, 0.1, 0.5)
Dp=sample(c(runif(M/2, 0, .5), runif(M/2, -0.5, 0)), M)
Dp=Dp[-1]

G=GenerateGenoDprimeRcpp(frq, Dp, N)
FQ=colMeans(G)/2

s=apply(G, 2, scale)
#Kcpp=CorMatrixRcpp(s)

#K
K=1/M * s %*% t(s)
KtK=K%*%t(K)
K2=K^2
diagK2=sum(diag(K2))
mVar=mean(1/(2*FQ*(1-FQ)))
ss=apply(s^2,1, sum)

me=var(K[row(K)<col(K)])
R2=(var(K[row(K)<col(K)])-1/M)*M/(M-1)
EK2=N*(N-1)/M+N+N/M*(mVar-2)+R2*N^2
SK2=sum(K2)

##################################
RmeMat=matrix(0, 30, 2)
for(i in 1:30) {
  Rme=MeVarK(ceiling(N/2), s)
  RmeMat[i,1]=Rme$me
  RmeMat[i,2]=Rme$me*N^2+N
}

bootS=10
EK2_=array(0, dim=c(bootS,2))
for(i in 1:bootS) {
  meFun=MeChi2(100, s)
  EK2_[i,1]=N+N/M*(mVar-2)+meFun$me*N^2
  EK2_[i,2]=meFun$me
}

K2S=array(0, dim=bootS)
for(i in 1:bootS) {
  K2S[i]=mean(MeSriram(100, s)$k2)
}

K2Ssub=array(0, dim=bootS)
sN=ceiling(N*0.75)
for(i in 1:bootS) {
  MeSsub=MeSriramSub(10, sN, G)
  K2Ssub[i]=(mean(MeSsub$k2[,1])-sN)/(sN^2)*N^2+N
}


vvM=matrix(0, bootS,1)
for(ii in 1:bootS){
  

SM=100
Sh2=matrix(0, SM, 1)
h2=0.01
Chi2Mat=matrix(0, SM, M)
for(i in 1:SM) {
  pp=0.5
  M1=ceiling(M*pp)
  beta=c(rnorm(M1, sd=sqrt(h2/M1)), rnorm(M-M1))
  bv=G%*%beta
  res=var(bv)*(1-h2)/h2
  Y=bv+rnorm(N, sd=sqrt(res))
  sY=scale(Y)
  sY=rnorm(N)
  yX=t(sY)%*%s
  Chi2Mat[i,]=(yX^2/N)
  yKy=t(sY)%*%K%*%sY
  yy=t(sY)%*%sY
  Sh2[i,1]=(yKy-N)/(SK2-N)
}
vvM[ii,1]=var(apply(Chi2Mat,1,sum))
}
