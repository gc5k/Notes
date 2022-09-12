
n=500 #sample size
m=1000 #marker
h2=0.05 #heritability
SM=100 #simulation
BS=5 #randomization factor

H2=matrix(0, SM, 2)
LK=0
TIMER=array(0, 2)

I=diag(1, n,n)
fq=runif(m, 0.5, 0.5)
x=matrix(0, n, m)
for(i in 1:m) {
  x[,i]=rbinom(n, 2, fq[i])
}
sx=apply(x, 2, scale)
K=sx%*%t(sx)/m
K2=K%*%K
K3=K2%*%K
K4=K3%*%K
k_3=m^3+8*m^2+m*(n^3+6*n^2+8*n-9)
#Randomized HE
Tstart=proc.time()

#evaluate LB
yKyM=array(0, dim=SM)
yKyV=array(0, dim=SM)
yKIyV=array(0, dim=SM)
yV=array(0, dim=SM)

for(i in 1:SM) {
  b=rnorm(m, 0, sqrt(h2/m)) #effect
  bv=sx%*%b
  
  y=2*(bv+rnorm(n, 0, sqrt(1-h2)))
  y=rnorm(n)
  y=scale(y)
  ys=scale(y)
  yV=var(y)
#  ys=scale(y)
  Lb=0
  for(j in 1:BS) {
    z=matrix(rnorm(n), n, 1)
    x1=t(sx)%*%z
    x2=sx%*%x1
    Lb=Lb+(t(x2)%*%x2)[1,1]
  }
  LK=Lb/(BS*m^2)

  m2=matrix(0, 2, 2)
  m2[1,1]=LK
  m2[1,2]=n
  m2[2,1]=n
  m2[2,2]=n
  wy=t(sx)%*%ys #yKy
  yIy=n-1
  yKy=mean(wy^2)
  yKyM[i]=yKy
  yKyV[i]=yKy
  yKIyV[i]=(t(y)%*%(K-I)%*%y)[1,1]
  yVec=matrix(c(yKy, yIy), 2, 1) #numerator
  B2=solve(m2)%*%yVec
  H2[i,2]=B2[1,1]/(B2[1,1]+B2[2,1])
}
print(paste(mean(yKyM), n*(1+n*h2/m)))
print(paste(var(yKyV), 2*mean(yV)^2*(n^2/m+n)))
#      , 2*(sum(diag(K4))*(h2/m)^2+2*h2/m*(1-h2/m)*sum(diag(K3))+(1-h2/m)^2*sum(diag(K2)))))
print(paste(var(yKIyV), 2*mean(yV)^2*(n^2/m)))

BV=matrix(bv, m, 1)
2*sum(diag(n^2/m+n))
2*(n^2/m)

Tend=proc.time()
TIMER[1]=Tend[3]-Tstart[3]

#old he
Tstart1=proc.time()

for(i in 1:SM) {
  y=sx%*%b+rnorm(n, 0, sqrt(1-h2))
  
  K=sx %*% t(sx)/m

  yM=y%*%t(y)
  
  yc=yM[col(yM) < row(yM)]
  Kc=K[col(K) < row(K)]
  md=lm(yc~Kc)
  H2[i, 1]=coefficients(md)[2]
}
Tend1=proc.time()
TIMER[2]=Tend1[3]-Tstart1[3]

layout(matrix(c(1,1,2,3), byrow=T, 2, 2))
barplot(TIMER,col=c("red", "blue"))
barplot(H2[,1],beside = T, col="red")
barplot(H2[,2],beside = T, col="blue")

abline(h=h2)
