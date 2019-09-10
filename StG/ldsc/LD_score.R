source("~/R/MyLib/shotgun.R")
REP=50
N = 500
M = 1000
h2 = 0.5

h2E=array(0, dim=c(2,REP))

for(rep in 1:REP) {
  fq = runif(M, 0.05, 0.95)
  Dl = array(0, dim=M)
  Dl = runif(M, 0.8, 0.9)
  Dl[seq(10, M, 10)] = 0
  
  G = GenerateGenoDprime(fq, Dl[1:(length(Dl)-1)], N)
  b = rnorm(M, 0, sqrt(h2/M))
  bv= G%*%b
  Y = bv + rnorm(N, 0, 1)
  chi = array(0, dim=M)
  for(i in 1:M) {
    sm = summary(lm(Y~G[,i]))
    chi[i] = sm$coefficients[2,3]^2
  }
  plot(chi, b^2*fq*(1-fq), pch=16)

  lds = array(0, dim=M)
  for(i in 1:(M/10)) {
    cg = cor(G[,((i-1)*10+1):(i*10)])
    for(j in 1:10) {
      lds[(i-1)*10+j] = sum(cg[,j]^2)
    }
  }

  ldSc=lm(chi~lds)
  h2E[1,rep]=ldSc$coefficients[2]*M/N
  h2E[2,rep]=var(bv)/var(Y)
}
plot(h2E[1,], h2E[2,])
