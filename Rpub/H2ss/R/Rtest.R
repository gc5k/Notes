n=500 #sample size
m=1000 #marker
h2=0.3 #heritability
h2d=0.3
b=rnorm(m, 0, sqrt(h2/m)) #effect
d=rnorm(m, 0, sqrt(h2d/m))
SIMU=5

#simu g
fq=runif(m, 0.3, 0.5)
x=matrix(0, n, m)
for(i in 1:m) {
  x[,i]=rbinom(n, 2, fq[i])
}
FQ=colMeans(x)/2
sA=apply(x, 2, scale)

K=sA%*%t(sA)/m
me=var(K[col(K)<row(K)])


##dominance
xd=matrix(0, n, m)

for(i in 1:m) {
  cd=c(0, 2*FQ[i], 4*FQ[i]-2)
  xd[,i]=cd[x[,i]+1]
}
dCt=matrix(rep(2*FQ^2, n), n, m, byrow = T)
dV=matrix(rep(sqrt(4*FQ^2*(1-FQ)^2), n), n, m, byrow = T)
sD=(xd-dCt)/dV
Kd=sD%*%t(sD)/m
med=var(Kd[col(K)<row(K)])

Dm=ifelse(x==1, 1, 0)
#simu y
H2=matrix(0, SIMU, 4)
for(i in 1:SIMU) {
  y=x%*%b+Dm%*%d
  vy=var(y)
  y=y+rnorm(n, 0, sqrt(vy/(h2+h2d)*(1-h2-h2d)))
  y=scale(y)
  yy=y%*%t(y)
  h2Mod=lm(yy[col(yy)<row(yy)]~K[col(yy)<row(yy)]+Kd[col(yy)<row(yy)])
  H2[i,1]=summary(h2Mod)$coefficients[2,1]
  H2[i,2]=summary(h2Mod)$coefficients[3,1]
  ss=matrix(0, m, 5)
  ssd=matrix(0, m, 5)
  for(j in 1:m) {
    mod=lm(y~sA[,j]+sD[,j])
    ss[j,1:4]=summary(mod)$coefficient[2,]
    ssd[j,1:4]=summary(mod)$coefficient[3,]
  }
  ss[,5]=ss[,3]^2
  H2[i,3]=((mean(ss[,5])-1)*n)/(n*n*me)
  ssd[,5]=ssd[,3]^2
  H2[i,4]=((mean(ssd[,5])-1)*n)/(n*n*med)
}
barplot(t(H2), beside = T)
abline(h=c(h2, h2d))
