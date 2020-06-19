M=100000 #marker
N=200 #sample size
P=runif(M, 0.2, 0.8) # ancestral
G=matrix(0, N, M)
for(i in 1:N) {
  G[i,] = rbinom(M, 2, P)
}
Y = rbinom(N, 1, 0.5)

Gs=apply(G, 2, scale)
GG=Gs %*% t(Gs) / M
EigenG=eigen(GG)
plot(EigenG$vectors[,1], EigenG$vectors[,2])

# get lambda GC
para = matrix(0, M, 2)
for (i in 1:M){
  model = lm(Y~Gs[,i])
  para[i,1]=summary(model)$coefficients[2,3]^2 #p value
  para[i,2]=summary(model)$coefficients[2,4] #p value
}
hist(para[,2])

pcut=c(0.001, 0.01, 0.05, 0.5, 0.95, 1)
layout(matrix(1:6, 3, 2))
for(i in 1:length(pcut)) {
  idx=which(para[,2] < pcut[i])
  Gs1=Gs[, idx]
  GG1=Gs1%*%t(Gs1)/ncol(Gs1)
  eG1=eigen(GG1)
  plot(eG1$vectors[,1], eG1$vectors[,2], col=(Y+1), main=paste("p-cut", pcut[i]))
}
