n1=1000
N=c(10000)#, 2000, 5000, 10000)
h2=0.01 #heritability
M=30000 #number of markers
Q=1000 #causal QTLs
n12=100

Simu=100
layout(matrix(1:(length(N)*2), 2, length(N)))
OBS=matrix(0, Simu, length(N))
LAM=matrix(0, Simu, length(N))

for(S in 1:Simu)
{
  for(I in 1:length(N))
  {
    N1=n1 #sample 1
    N2=N[I] #sample 2
    rho=n12/sqrt(N1*N2)
    
    E1=rnorm(M)
    E2=rho*E1 +sqrt(1-rho^2)*rnorm(M)

    B=rnorm(Q, 0, sqrt(h2/Q)) #sampling genetic effect
    e1=E1/sqrt(N1) #residual 1
    e2=E2/sqrt(N2) #residual 2

    b1=e1
    b1[1:Q] = e1[1:Q]+B
    
    b2=e2
    b2[1:Q] = e2[1:Q]+B
    
    Ts=(b1-b2)^2/(e1^2+e2^2)
    Lm=median(Ts)
    LAM[S, I] = (1-Lm)/(2*sqrt(N1*N2)/(N1+N2))
    Z1=sqrt(N1)*b1
    Z2=sqrt(N2)*b2
    P1=ifelse(Z1 < 0, 2*pnorm(Z1), 2*(1-pnorm(Z1)))
    P2=ifelse(Z2 < 0, 2*pnorm(Z2), 2*(1-pnorm(Z2)))
    
    obsCov=cor(Z1, Z2)
    OBS[S, I]=obsCov
    COV_h2=sqrt(N1*N2)/M*h2
    V1=N1/M*h2+1
    V2=N2/M*h2+1
    
    COV_n12=rho/(sqrt(N1/M*h2+1)*sqrt(N2/M*h2+1))
    preCov=COV_h2/sqrt(V1*V2)+COV_n12
    
    
    V=matrix(c(1, obsCov, obsCov, 1), 2, 2)
    IV=solve(V)
    
    X=matrix(0, length(Z1), 1)
    for(i in 1:length(Z1))
    {
      z=matrix(c(Z1[i], Z2[i]), 1, 2)
      zT=t(z)
      X[i,1] = z %*% IV %*% zT
    }
    hist(pchisq(X[1:Q,1], 2, lower.tail=F), xlab="P-value", main=paste("N=", N[I], ", M=", M, "\nrho=", obsCov))
    hist(pchisq(X[(Q+1):length(Z1),1], 2, lower.tail=F), xlab="P-value", main=paste("N=", N[I], ", M=", M, "\nrho=", obsCov))  
  }
}
print(paste(obsCov, " ", preCov))
print(paste(mean(OBS[,1]), " ", sd(OBS[,1])))
print(paste(mean(LAM[,1]), " ", sd(LAM[,1])))
