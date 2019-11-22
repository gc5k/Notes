lk<-function(y, x, beta) {
  ll=t(y)%*%x%*%beta-matrix(1, 1, length(y)) %*% log(1+exp(x%*%beta))
  return(ll)
}

SIMU=1000
n=100
b=c(0.1,1,2)

epi=0.001
MAX_It=10

BETA=matrix(0, SIMU, length(b)+1)
for(ii in 1:SIMU) {
  x=matrix(1, n, length(b))
  x[,2]=rnorm(n)
  x[,3]=rbinom(n, 2, 0.5)
  y=rbinom(n, 1, exp(x%*%b)/(1+exp(x%*%b)))

  b0=c(1, 1, 1)
  b1=c(1, 1, 1)
  l1=1
  l0=2
  maxIt=0
  while(abs(l1-l0) > epi & maxIt<MAX_It) {
    l0=lk(y, x, b0)
    pb=exp(x%*%b0)/(1+exp(x%*%b0))
    m=diag(c(pb*(1-pb)), n, n)
    XtX=t(x)%*%m%*%x
    XtX_I=solve(XtX)
    yres=y-pb
    b1=b0+XtX_I%*%t(x)%*%yres
    l1=lk(y, x, b1)
    b0=b1
    maxIt=maxIt+1
  } 
  BETA[ii,]=c(b1, maxIt)
}
colMeans(BETA)
