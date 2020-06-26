dat=read.table("hendersonMat.txt", as.is = T, header = F)
dat=as.matrix(dat)
Lt=chol(dat)
L=t(Lt)

D=diag(L[col(L)==row(L)], nrow(L), ncol(L))

Tm=as.matrix(read.table("hendersonT.txt", as.is = T, header = F))

Tmt=t(Tm)
V=Tm%*%D%*%t(D%*%Tm)
solve(V)
