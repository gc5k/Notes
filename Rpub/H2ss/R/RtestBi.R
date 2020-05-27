S=30

SIMU=array(0, dim=c(S, 4))
Ey=array(0, dim=c(S, 2))

for(s in 1:S) {

hT1=0.5
hT2=0.5
rg=1
re=0.3
nc=999
n1=1
n2=1

N1=nc+n1
N2=nc+n2

m=500 #marker
frq=runif(m, 0.1, 0.5)

x=matrix(rbinom(m*nc, 2, frq), nc, m, byrow=T)
x1=matrix(rbinom(m*n1, 2, frq), n1, m, byrow = T)
x2=matrix(rbinom(m*n2, 2, frq), n2, m, byrow = T)

X1=rbind(x, x1)
X2=rbind(x, x2)

b1=rnorm(m)
b2=rg*b1+sqrt(1-rg^2)*rnorm(m)

e0=rnorm(nc)
e1=rnorm(n1)
e2=rnorm(n2)
E1=c(e0, e1)
E2=c(re*e0+sqrt(1-re^2)*rnorm(nc), e2)

bv1=X1%*%b1
vG1=var(bv1[,1])
y1=bv1+sqrt(vG1/hT1*(1-hT1))*E1
sy1=scale(y1)
sy1c=sy1[1:nc]

bv2=X2%*%b2
vG2=var(bv2[,1])
y2=bv2+sqrt(vG2/hT2*(1-hT2))*E2
sy2=scale(y2)
sy2c=sy2[1:nc]

ss1=matrix(0, m, 4)
ss2=matrix(0, m, 4)
for(i in 1:m) {
  mod1=lm(sy1~X1[,i])
  ss1[i,]=summary(mod1)$coefficients[2,]
  
  mod2=lm(sy2~X2[,i])
  ss2[i,]=summary(mod2)$coefficients[2,]
}

Q=(mean(ss1[,3]*ss2[,3])-nc/sqrt(N1*N2)*mean(sy1c*sy2c)) / (sqrt(N1*N2)/m)

EQ=((sqrt(N1*N2)*rg*sqrt(hT1*hT2)/m+nc/sqrt(N1*N2)*re) - nc/sqrt(N1*N2)*re) / (sqrt(N1*N2)/m)
hT1est=(mean(ss1[,3]^2)-mean(sy1*sy1))/(N1/m)
hT2est=(mean(ss2[,3]^2)-mean(sy2*sy2))/(N2/m)
rG=Q/sqrt(hT1est*hT2est)
rGEst=EQ/sqrt(hT1est*hT2est)

SIMU[s, ]=c(hT1est, hT2est, rG, rGEst)

sX1=scale(X1)
sX2=scale(X2)
K=sX1%*%t(sX2)/m
I=diag(1, N1, N2)
Ey[s,1]=t(sy1)%*%(K-I)%*%sy2

}
trK2=N1*N2/m+N1
E_Ey=(trK2-N1)*rg*sqrt(hT1*hT2)

Sig1=sy1%*%t(sy1)
Sig2=sy2%*%t(sy2)
V_Ey=sum(diag((K-I)%*%(K-I)))*2
plot(Ey[,1])
abline(h=E_Ey)
colMeans(SIMU)

print(paste(mean(Ey[,1]), E_Ey))
print(paste(var(Ey[,1]), V_Ey))

