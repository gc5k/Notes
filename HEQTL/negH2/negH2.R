#see 10.1534/genetics.120.303161, equation 6
library(mnormt)
m=1000
n=500
frq=0.5
h2=0.01
he=1-h2
phi=h2/he
alpha=0.5 

G=matrix(rbinom(m*n, 2, frq), n, m)
sG=scale(G)
A=sG%*%t(sG)/m
G1=solve(phi*A+diag(n))
G2=alpha*(A-diag(n))
GG=solve(G1+G2)

REP=30

HEG=matrix(0, REP, 2)
for(s in 1:REP) {
  PARA=matrix(0, m, 4)
  y=rmnorm(1, mean=rep(0, n), GG)
  for(i in 1:m) {
    mod=lm(y~G[,i])
    PARA[i,]=summary(mod)$coefficient[2,]
  }
  qqplot(PARA[,3]^2, rchisq(m,1))
  abline(a=0, b=1, col="red")
  HEG[s,1]=median(PARA[,3]^2)
  HEG[s,2]=mean(PARA[,3]^2)
}
