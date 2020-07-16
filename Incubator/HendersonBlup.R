dat=read.table("hendersonMat.txt", as.is = T, header = F)
dat=as.matrix(dat)
Lt=chol(dat)
L=t(Lt)
tm=L
diag(tm)=1
D=diag(L[col(L)==row(L)], nrow(L), ncol(L))

#chol factorization
ch <- chol(dat)
dd <- diag(ch)
Lf <- t(ch/dd)
DD <- diag(dd, nrow(Lf), ncol(Lf))
Lf%*%DD%*%DD%*%t(Lf)

solve(t(Lf))%*%diag(1/dd^2, nrow(Lf), ncol(Lf))%*%solve(Lf)
