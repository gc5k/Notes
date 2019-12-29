source("buildU.R")
set.seed(1000)
sigma = c(0, 1, 2)
sl = length(sigma)
m = 100
n = 800
logSN = ceiling(log(n, sl))
bs=c(sl^(seq(logSN-1,0)))

A = matrix(sample(sigma, m * n, replace = TRUE), nrow = m, ncol = n)

N=sl^logSN
xDim=10
x1 = matrix(c(runif(n*xDim, -5, 5), rep(0, xDim*(N-n))), xDim, N, byrow = F)

if (m < n) {
  M=ifelse(m<=logSN, logSN, ceiling(m/logSN)*logSN)
}
BLK=ceiling(M/logSN)
A=rbind(A, matrix(0, nrow=M-m, ncol=n))
A=cbind(A, matrix(0, nrow=M, ncol=N-n))

AX=matrix(0, nrow(A), nrow(x1))
for(Bk in 1:BLK) {
  Pcompact=buildPF2(A[((Bk-1)*logSN+1):(logSN*Bk),], bs, n)

  for(j in 1:xDim) {
    Ax1=mailman_product2(logSN, N, Pcompact, sigma, x1[j,])
    AX[((Bk-1)*logSN+1):(logSN*Bk), j] = Ax1
#    Ax2= A[((Bk-1)*logSN+1):(logSN*Bk),] %*% x1[j,]
#    print(paste("Bk", Bk))
#    print(t(Ax1))
#    print(t(Ax2))
#    print(paste("A %*% X1 == Ax1: ", all.equal(Ax2, Ax1)))
  }
}
BX=A %*% t(x1)
