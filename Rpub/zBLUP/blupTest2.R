M=200
N=200
ha=0.5
hd=0.3

frq=rep(0.5, M)
G=matrix(rbinom(M*N, 2, frq), N, M)
Gd=matrix(ifelse(G==1, 1, 0), N, M)
a=rnorm(M)
d=rnorm(M)
BVa=G%*%a
BVd=Gd%*%d

Beta=matrix(c(1, 2), 2, 1)
X=matrix(rbinom(2*N, 2, 0.5), N, 2)
vBVa=var(BVa)[1,1]
vBVd=var(BVd)[1,1]
ve=vBVa+vBVd
y=X%*%Beta+BVa+BVd+rnorm(N, 0, sqrt(ve))

#MME
C11=t(X)%*%X
C12=t(X)
C21=X
C13=t(X)
C31=X
C22=diag(1, M)+1.5*diag(1, M)
C33=diag(1, M)+3*diag(1, M)
MME_1=cbind(C11, C12, C13)
MME_2=cbind(C21, C22, diag(1, M))
MME_3=cbind(C31, diag(1, M), C33)

MME_mat=rbind(MME_1, MME_2, MME_3)
MME_y=matrix(c(t(X)%*%y, y, y), 2*M+2, 1)
MME_b=solve(MME_mat)%*%MME_y
plot(BVa, MME_b[3:(M+2),1])
abline(a=0, b=1)
plot(BVd, MME_b[(M+3):(nrow(MME_b)),1])

##GLMM
V=(diag(vBVa, N)+diag(vBVd, N)+diag(ve, N))
VI=solve(V)
bEst=solve(t(X)%*%VI%*%X)%*%t(X)%*%VI%*%y

uA=vBVa*VI%*%(y-X%*%bEst)
uD=vBVd*VI%*%(y-X%*%bEst)

plot(uA, MME_b[3:(M+2),1])
plot(uD, MME_b[(M+3):(nrow(MME_b)),1])

uD2=(vBVd/vBVa)*uA
