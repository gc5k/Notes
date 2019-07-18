source("~/R/MyLib/shotgun.R")
M=2000 #marker
N=1000 #snp
h2=0.5 #snp-heritability
#fq=c(0.5, 0.3)
fq=runif(M, 0.1, 0.5)
#ld=c(0,0)
ld=rep(0, M)
ld=runif(M, -0.9, 0.9)

G=GenerateGenoDprime(fq, ld, N)

bg=G[,1]*1 #first mark is qtl
y=bg+rnorm(N, 0, sqrt(var(bg)/h2*(1-h2)))
y=scale(y)  #standardization
Y=matrix(y, nrow=length(y), ncol=1)%*%matrix(y, nrow=1, ncol=length(y))
YY=matrix(Y[row(Y) < col(Y)], nrow=N*(N-1)/2, 1)

sG=apply(G, 2, scale) #standardization snp

sGRM=sG %*% t(sG)/M
BG=sGRM[row(sGRM) < col(sGRM)]
X2=matrix(0, nrow=length(BG), 2)
X2[,1]=1
X2[,2]=BG
XtX2=t(X2)%*%X2

beta=solve(XtX2)%*%t(X2)%*%YY

me=1/var(BG)
NN=N*(N-1)/2
detA=XtX2[1,1]*XtX2[2,2]-XtX2[1,2]*XtX2[2,1]
xtx=XtX2
xtx[1,1]=XtX2[2,2]
xtx[2,2]=XtX2[1,1]
xtx[1,2]=-1*XtX2[2,1]
xtx[2,1]=-1*XtX2[1,2]
mInt_x=xtx
int_X=solve(XtX2)

XtX_xT=mInt_x %*%t(X2)
bb=XtX_xT %*% YY
bbF=bb/detA

xy=t(X2) %*% YY
bbS=int_X %*% xy
