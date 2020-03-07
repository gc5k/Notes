mean=0.5
sd=1
M=1/dnorm(0,mean,sd)
n=10000
a=-5
b=5
X=array(0, dim=n)
cnt=1
Qcnt=1
while(cnt<n) {
  p=runif(1)
  x=runif(1,a,b)
  if(p < (1/M*dnorm(x, mean, sd))) {
    X[cnt]=x
    cnt=cnt+1
  }
  Qcnt=Qcnt+1
}
print(cnt/Qcnt)
hist(X, breaks = 50)


##mcmc
mu=0.1
sig=1
n=10000
X=array(0, dim=n)
Xp=array(0, dim=n)
pp=array(0, dim=n)
X[1]=rnorm(1)
Xp[1]=1
pp[1]=1
cnt=2
while(cnt<n) {
  p=runif(1)
  xt=X[cnt-1]
  xt1=rnorm(1, xt)
  pcut=min(1, dnorm(xt1, mean=mu, sd=sig)/dnorm(xt, mean=mu, sd=sig))
  Xp[cnt]=pcut
  pp[cnt]=p
  if(p > pcut) {
    X[cnt]=xt
  } else {
    X[cnt]=xt1
  }
  cnt=cnt+1
}
hist(X, breaks = 50)
