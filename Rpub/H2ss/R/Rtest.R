n=500 #sample size
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
k4Alt=k2^2
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
print(paste("sum(diag(K4))=", sum(diag(k4)), "sum[(kkT)^2]", sum(k4Alt)))


Sum_D_k2_sq=sum(diag(k2))^2

ZX=n^2*m^4+(2*n^3+2*n^2+8*n)*m^3+(n^4+2*n^3+21*n^2+20*n)*m^2+(8*n^3+20*n^2+20*n)*m

print(paste(sum(diag(k2)^2)*n, ZX/m^4, sum(diag(k2))^2))

EU=n*m^2+(n^2+n)*m
Sum_D_k2=sum(diag(k2))
print(paste(EU/m^2, Sum_D_k2))

k2_sq=n^2*m^4+(2*n^2+2*n+8)*m^3+(4*n+40)*m^2+48*m
k2_sq=n^2*m^4+(2*n^2+2*n)*m^3+(4*n)*m^2

sum(diag(k2)^2)
mean(diag(k2)^2)

ZX/m^4/n
(n^2*m^4+(2*n^2+6*n-8)*m^3+(24*n+20)*m^2+48*m)/(n*m^4)
((3*n^2+2*n-4)*m^4+(8*n+4)*m^3+(20*n+24)*m^2+48*m)/(n^2*m^4)

Uii=mean(diag(k))
Vii=var(diag(k))
Uij=mean(k[col(k)<row(k)])
Vij=var(k[col(k)<row(k)])

A1=Uii
A2=Uii^2+Vii
A3=Uii^3+3*Uii*Vii
A4=Uii^4+6*Uii^2*Vii+3*Vii^2

B1=Uij
B2=Uij^2+Vij
B3=Uij^3+3*Uij*Vij
B4=Uij^4+6*Uij^2*Vij+3*Vij^2

mean(diag(k4))

pK4D=A4+(n-1)*B4+2*(n-1)*A2*B2+(n-1)*(n-2)*B2*B2
print(pK4D)
mean(diag(k4Alt))

pK4O=2*A2*B2+(n-2)*B2*B2+2*A1*B2*A1+(n-2)*A1*B1^3+(n-2)*B1^3*A1+(n-1)*(n-2)*B1^4
print(pK4O)
mean(k4Alt[col(k4Alt)<row(k4Alt)])

pK4D*n+n*(n-1)*pK4O
sum(k4Alt)

