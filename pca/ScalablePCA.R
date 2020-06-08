#row dimension vs column dimension
RD=1

m=1000
n1=250
n2=250
n=n1+n2

fst=0.2
frq=runif(m, 0.2, 0.8)

frq1=rbeta(m, frq*(1-fst)/fst, (1-frq)*(1-fst)/fst)
G1=matrix(rbinom(n1*m, 2, frq1), n1, m, byrow = T)
frq2=rbeta(m, frq*(1-fst)/fst, (1-frq)*(1-fst)/fst)
G2=matrix(rbinom(n2*m, 2, frq2), n2, m, byrow = T)

G=rbind(G1, G2)
plot(colMeans(G)/2, frq)

if(RD==1) {
  Y=t(scale(G))
} else {
  Y=scale(G)
}
yy=t(Y)%*%Y/ncol(Y)

eY=eigen(yy)
svdY=svd(Y)
k=10

C0=matrix(rnorm(nrow(Y)*k), nrow(Y), k)
for(i in 1:5) {
  #E-step
  inv_CC=solve(t(C0)%*%C0)
  X=inv_CC%*%t(C0)%*%Y
  
  #M-step
  XXt=X%*%t(X)
  inv_XXt=solve(XXt)
  C1=Y%*%t(X)%*%inv_XXt
  print(cor(C1[,1], C0[,1]))
  C0=C1
}
Ye=C0%*%X
plot(Ye[,1], svdY$u[,1])
