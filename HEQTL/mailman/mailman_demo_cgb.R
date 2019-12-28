source("buildU.R")
set.seed(1000)
sigma = c(0, 1, 2)
sl = length(sigma)
m = 100
n = 800
logSN = ceiling(log(n, sl))
bs=c(sl^(seq(logSN-1,0)))

A = matrix(sample(sigma, m * n, replace = TRUE), nrow = m, ncol = n)

x1 = c(runif(n, -5, 5))
N=sl^logSN
X1=append(x1, rep(0, N-n))

if (m < n) {
  M=ifelse(m<=logSN, logSN, ceiling(m/logSN)*logSN)
}
BLK=ceiling(M/logSN)
A=rbind(A, matrix(0, nrow=M-m, ncol=n))
A=cbind(A, matrix(0, nrow=M, ncol=N-n))

#U = buildU(logSN, sigma)

for(Bk in 1:BLK) {
#  Pcompact=buildPF(A[((Bk-1)*logSN+1):(logSN*Bk),], bs)
  Pcompact=buildPF2(A[((Bk-1)*logSN+1):(logSN*Bk),], bs, n)

  Ax1=mailman_product2(logSN, N, Pcompact, sigma, X1)
  Ax2= A[((Bk-1)*logSN+1):(logSN*Bk),] %*% X1

  print(paste("Bk", Bk))
  print(t(Ax1))
  print(t(Ax2))
  print(paste("A %*% X1 == Ax1: ", all.equal(Ax2, Ax1)))
}
