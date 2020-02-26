library(Rcpp)
sourceCpp("~/git/Notes/R/RLib/Shotgun.cpp")
M=10000
N=500
frq=runif(M, 0.5, 0.5)
Dp=sample(c(runif(M/2, 0, 0), runif(M/2, 0, 0)), M)
Dp=Dp[-1]


G=GenerateGenoDprimeRcpp(frq, Dp, N)
FQ=colMeans(G)/2
plot(colMeans(G)/2, frq)

#cr=cor(G)
#cr2=cr^2

#Ecr2=mean(cr2[row(cr2) < col(cr2)])
#Ocr2=mean(cr2[row(cr2) != col(cr2)])
#standardization
s=apply(G, 2, scale)
Kcpp=CorMatrixRcpp(s)

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


Sh2=matrix(0, 100, 1)
for(i in 1:100) {
  h2=0.5
  pp=0.5
  M1=ceiling(M*pp)
  beta=c(rnorm(M1, sd=sqrt(h2/M1)), rnorm(M-M1))
  bv=G%*%beta
  res=var(bv)*(1-h2)/h2
  Y=bv+rnorm(N, sd=sqrt(res))
  y=scale(Y)

  yX=t(y)%*%s

  yKy=t(y)%*%K%*%y
  yy=t(y)%*%y
  Sh2[i,1]=(yKy-N)/(SK2-N)
  print(h2)
}
plot(Sh2)



##################################
RmeMat=matrix(0, 100, 1)
for(i in 1:100) {
  Rme=RMe(100, s)
  RmeMat[i,1]=Rme$me
}

bootS=10
EK2_=array(0, dim=bootS)
for(i in 1:bootS) {
  meFun=MeFun(N/4, G)
  EK2_[i]=N*(N-1)/M+N+N/M*(mVar-2)+meFun$R2*N^2
}

MeFun <-function(sn, G) {
  n=nrow(G)
  m=ncol(G)
  SP=sample(n, sn)
  Gs_=G[SP,]
  s=apply(Gs_, 2, scale)
  K_=1/m * s %*% t(s)
  K2_=K_^2
  FQ_=colMeans(Gs_)/2
  mVar_=mean(1/(2*FQ_*(1-FQ_)))

  me_=var(K_[row(K_)<col(K_)])
  R2_=(var(K_[row(K_)<col(K_)])-1/m)*m/(m-1)
  EK2_=sn*(sn-1)/m+sn+sn/m*(mVar_-2)+R2_*sn^2
  SK2_=sum(K2_)
  
  pa=list("EK2"=EK2_, "SK2"=SK2_, "R2"=R2_, "me"=me_, "mVar"=mVar_)
}

RMe <-function(B, sG) {
  z=matrix(rnorm(B*nrow(sG)), B, nrow(sG))
  v1=(z%*%s)^2/N
  qqplot(v1[1,], rchisq(M, 1))
  abline(a=0, b=1, col="red")
  vv=apply(v1, 1, sum)
  v_vv=var(vv)
  me_v=v_vv/(2*M^2)-1/(N-1)
  pa=list("me"=me_v)
}

#simu
Gc=GenerateGenoDprimeRcpp(frq, Dp, N)
Gn=GenerateGenoDprimeRcpp(frq, Dp, N)
plot(colMeans(Gc)/2, frq)

sGc=apply(Gc, 2, scale)
sGn=apply(Gn, 2, scale)

grm_c=sGc%*%t(sGc)/M
GRM_C=CorMatrixRcpp(sGc)
GRM_N=CorMatrixRcpp(sGn)
GRM_CN=CorMatrix2Rcpp(sGc, sGn)


GRM_C=matrix(c(0.81, 0, 0, 0.81), 2, 2, byrow = T)
GRM_N=matrix(0, 3, 3)
diag(GRM_N)=c(.01, 1.61, 1.61)
GRM_N[2,3]=GRM_N[3,2]=-1.6
GRM_CN=matrix(0, 2, 3)
GRM_CN[,2]=0.8
GRM_CN[,3]=-0.8

IC=solve(GRM_C)
Pcn=IC %*% GRM_CN
Pnc=t(GRM_CN) %*% IC
Pcn_Gcn=t(Pcn) %*% GRM_CN
Mnn=diag(diag(GRM_N)-diag(Pcn_Gcn), nrow=nrow(GRM_N), ncol=nrow(GRM_N))
IMnn=solve(Mnn)

v1_1=rbind(IC, matrix(0, nrow=nrow(GRM_N), ncol=ncol(GRM_C)))
v1_2=rbind(-1*Pcn%*%IMnn, IMnn)
V1=cbind(v1_1, v1_2)

v2_1=rbind(diag(1, nrow=nrow(GRM_C), ncol=nrow(GRM_C)), -1*Pnc)
v2_2=rbind(matrix(0, nrow=nrow(GRM_C), ncol=nrow(GRM_N)), diag(1, nrow(GRM_N), ncol(GRM_N)))
V2=cbind(v2_1, v2_2)

IV=V1 %*% V2
IVI=solve(IV)

ld=Dprime2LDRcpp(frq, Dp)
c2=Dprime2CorRcpp(frq, Dp)

G1=GenerateGenoLDRcpp(frq, ld, N)
plot(colMeans(G1)/2, frq)

G2=GenerateHapLDRcpp(frq, ld, N)
plot(colMeans(G2), frq)

FQ=colMeans(G)/2
layout(matrix(1:2,1,2))
plot(colMeans(G)/2, colMeans(G2)/2)
plot(frq, colMeans(G)/2)
plot(frq, colMeans(G2)/2)

cr=cor(G)

cr2=cr^2

Ecr2=mean(cr2[row(cr2) < col(cr2)])
Ocr2=mean(cr2[row(cr2) != col(cr2)])
#standardization
s=apply(G, 2, scale)
Kcpp=CorMatrix(s)

#K
K=1/M * s %*% t(s)
K2=K^2
diagK2=sum(diag(K2))
mVar=mean(1/(2*FQ*(1-FQ)))
ss=apply(s^2,1, sum)

me=var(K[row(K)<col(K)])
R2=(var(K[row(K)<col(K)])-1/M)*M/(M-1)
EK2=N*(N-1)/M+N+N/M*(mVar-2)+R2*N^2
SK2=sum(K2)

z=rnorm(N)
Me=
r2=(diagK2*M^2/N-sVar-M*(M-1))/(M*(M-1))
#Moments
s2=s^2

#validate var & cov
s2_v=var(s2)

s2_cor=cor(s2)
hist(diag(s2_v), main=paste("var=", mean(diag(s2_v))))
plot(diag(s2_v)+1, 1/(2*FQ*(1-FQ)))
abline(a=0, b=1, col="red")

s2_cov=cov(s2)
s2_crE=matrix(0, M, M)
for(i in 1:M) {
  for(j in 1:i) {
    s2_crE[i,j]=cr[i,j]^2+cr[i,j]*(FQ[i]-(1-FQ[i]))*(FQ[j]-(1-FQ[j]))/(2*sqrt(frq[i]*(1-frq[i])*frq[j]*(1-frq[j])))
  }
}

l2=s2_crE[row(s2_crE) > (col(s2_crE))]
l1=s2_cov[row(s2_cov) > (col(s2_cov))]
idx=which(abs(l2) < 0.2)
plot(l1[idx], l2[idx])
abline(a=0, b=1, col="red")

plot(s2_cov[row(s2_cov) > (col(s2_cov))], s2_crE[row(s2_crE) > (col(s2_crE))], cex=0.3, pch=16)
abline(a=0, b=1, col="red", lwd=2)

Est_cr2=mean(s2_crE[row(s2_crE ) > col(s2_crE)])
mK2_diag=mean(diag(K2))
eK2_diag=
#validate s4
s4=s^4
plot(colMeans(s4), 1/(2*frq*(1-frq)))
abline(a=0, b=1, col="red", lwd=2)

Me <-function(s) {
  z=matrix(rnorm(nrow(s)), nrow=1, ncol=nrow(s))
  t2=(z%*%s)^2/nrow(s)
}
