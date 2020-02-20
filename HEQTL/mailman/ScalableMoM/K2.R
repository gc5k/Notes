source("~/R/MyLib/shotgun.R")
M=1000
N=500
frq=runif(M, 0.1, 0.5)
#frq=runif(M, 0.1, 0.5)
Dp=sample(c(runif(M/2, 0.5, 1), runif(M/2, -1, -0.5)), M)

#simu
G=GenerateGenoDprime(frq, Dp, N)

FQ=colMeans(G)/2
plot(frq, 1-colMeans(G)/2)
cr=cor(G)

cr2=cr^2

Ecr2=mean(cr2[row(cr2) < col(cr2)])
Ocr2=mean(cr2[row(cr2) != col(cr2)])
#standardization
s=apply(G, 2, scale)

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
