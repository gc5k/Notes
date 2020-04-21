n=500 #sample size
m=1000 #marker
h2=0.3 #heritability

b=rnorm(m, 0, sqrt(h2/m)) #effect
SIMU=50

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

sW=x-2*FQ
Kw=sW%*%t(sW)/sum(2*FQ*(1-FQ))
meW=var(Kw[col(Kw)<row(Kw)])

#simu y
H2=matrix(0, SIMU, 4)
H2w=matrix(0, SIMU, 2)
for(i in 1:SIMU) {
  y=x%*%b
  vy=var(y)
  y=y+rnorm(n, 0, sqrt(vy/(h2)*(1-h2)))
  y=scale(y)
  yy=y%*%t(y)
  h2Mod=lm(yy[col(yy)<row(yy)]~K[col(yy)<row(yy)])
  H2[i,1]=summary(h2Mod)$coefficients[2,1]
  ss=matrix(0, m, 5)
  for(j in 1:m) {
    mod=lm(y~sA[,j])
    ss[j,1:4]=summary(mod)$coefficient[2,]
  }
  ss[,5]=ss[,3]^2
  H2[i,2]=((mean(ss[,5])-1)*n)/(n*n*me)
  
  h2wMod=lm(yy[col(yy)<row(yy)]~Kw[col(yy)<row(yy)])
  H2[i,3]=summary(h2wMod)$coefficients[2,1]
  chiW=ss[,5]*2*FQ*(1-FQ)
  H2[i,4]=((mean(chiW)/mean(2*FQ*(1-FQ))-1)*n)/(n*n*meW)
}
barplot(t(H2), beside = T)
abline(h=c(h2))
legend("topleft", legend=c("h2", "ssh2", "h2w", "ssh2w"))

ME=1/me
##delta
s=sqrt(2*ME/n^2+2*ME^3/n^5)
