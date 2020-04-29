hT1=0.3
hT2=0.5
rg=0.5
re=0.3
n=1000 #sample size
m=500 #marker
frq=runif(m, 0.1, 0.5)

x=matrix(rbinom(m*n, 2, frq), n, m, byrow=T)

b1=rnorm(m)
b2=rg*b1+sqrt(1-rg^2)*rnorm(m)

e1=rnorm(n)
e2=re*e1+sqrt(1-re^2)*rnorm(n)

bv1=x%*%b1
vG1=var(bv1[,1])
y1=bv1+sqrt(vG1/hT1*(1-hT1))*e1
sy1=scale(y1)

bv2=x%*%b2
vG2=var(bv2[,1])
y2=bv2+sqrt(vG2/hT2*(1-hT2))*e2
sy2=scale(y2)

ss1=matrix(0, m, 4)
ss2=matrix(0, m, 4)
for(i in 1:m) {
  mod1=lm(sy1~x[,i])
  ss1[i,]=summary(mod1)$coefficients[2,]
  
  mod2=lm(sy2~x[,i])
  ss2[i,]=summary(mod2)$coefficients[2,]
}

rh2=(mean(ss1[,3]*ss2[,3])-mean(sy1*sy2))/(n/m)
hT1est=(mean(ss1[,3]^2)-mean(sy1*sy1))/(n/m)
hT2est=(mean(ss2[,3]^2)-mean(sy2*sy2))/(n/m)

rGEst=rh2/sqrt(hT1est*hT2est)
