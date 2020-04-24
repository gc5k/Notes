
n=500 #sample size
m=2000 #marker
h2=0.3 #heritability
SM=1 #simulation
BS=30 #randomization factor
b=rnorm(m, 0, sqrt(h2/m)) #effect

H2=matrix(0, SM, 2)
LK=0
TIMER=array(0, 2)

I=diag(1, n,n)
fq=runif(m, 0.1, 0.5)
x=matrix(0, n, m)
for(i in 1:m) {
  x[,i]=rbinom(n, 2, fq[i])
}
sx=apply(x, 2, scale)

#Randomized HE
Tstart=proc.time()

#evaluate LB

for(i in 1:SM) {
  y=sx%*%b+rnorm(n, 0, sqrt(1-h2))
  ys=scale(y)
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
  yVec=matrix(c(yKy, yIy), 2, 1) #numerator
  B2=solve(m2)%*%yVec
  H2[i,2]=B2[1,1]/(B2[1,1]+B2[2,1])
}
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
