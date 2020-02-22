u=1
sg=2
y=rnorm(50, u, sg)

burnin=1000
Cyc=10000
Gibbs=matrix(0, Cyc, 2)
for(i in 1:Cyc) {
  sd0=sd(y)/sqrt(length(y))
  u0=mean(y)
  
  u1=rnorm(1,u0, sd0)
  sl=sum((y-u1)^2)
  sd1=sl*1/rchisq(1,df = length(y)-2)
  Gibbs[i,]=c(u1, sd1)
}
print(colMeans(Gibbs))
layout(matrix(1:2, 2, 1))
plot(Gibbs[,1], cex=0.3, pch=16)
plot(Gibbs[,2], cex=0.3, pch=16)

