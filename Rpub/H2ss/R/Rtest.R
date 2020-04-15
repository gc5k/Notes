n=200 #sample size
m=500 #marker
SM=50 #simulation
BS=25 #randomization factor

LK=0

fq=runif(m, 0.1, 0.5)
x=matrix(0, n, m)
for(i in 1:m) {
  x[,i]=rbinom(n, 2, fq[i])
}
sx=apply(x, 2, scale)
k=sx%*%t(sx)/m
ksq=k^2
k2=k%*%t(k)
k4=k2%*%t(k2)
#evaluate LB

LK=array(0, SM)
for(i in 1:SM) {
  Lb=array(0, BS)
  for(j in 1:BS) {
    z=matrix(rnorm(n), n, 1)
    x1=t(sx)%*%z
    x2=sx%*%x1
    Lb[j]=(t(x2)%*%x2)[1,1]
  }
  LK[i]=sum(Lb)/(BS*m^2)
}
print(paste("Obs mean", mean(LK), "vs expected mean", n^2/m+n)) 
print(paste("obs var", var(LK), " vs expected var", 2*sum(diag(k4))/BS))
print(paste("incorrected var in Wu's paper", (n^2/m+n)/BS))
