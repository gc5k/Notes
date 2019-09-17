N=1000
Mc=1
fc=runif(Mc, 0.1, 0.5)
Mr=10

B=rnorm(Mc+Mr)
h2=0.5

gMat=matrix(0, N, Mc+Mr)
for(i in 1:Mc) {
  gMat[,i]=rbinom(N, 2, fc[i])
}
for(i in (Mc+1):(Mc+Mr)) {
  gMat[1,i]=1
}

sG=apply(gMat, 2, scale)
G=sG %*% t(sG)/(Mc+Mr)

gV=gMat %*% B

REP=100

REG=matrix(0, REP, 2)
for(i in 1:REP) {
  y=gV+rnorm(N, 0, sd=sqrt(var(gV)/h2*(1-h2)))
  y=scale(y)
  yM=y %*% t(y)
  
  yy=yM[row(yM) < col(yM)]
  gg=G[row(G) < col(G)]
  
  HE=lm(yy~gg)
  REG[i,] = coefficients(HE)
}
colMeans(REG)
