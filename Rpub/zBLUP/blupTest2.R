#walsh example
y=matrix(c(9, 12, 11, 6, 7, 14), 6, 1)
x=matrix(0, 6, 2)
x[,1]=c(1, 0, 1, 1, 1, 0)
x[,2]=c(0, 1, 0, 0, 0, 1)
z=matrix(0, 6, 3)
z[,1]=c(1,1,0,0,0,0)
z[,2]=c(0,0,1,1,0,0)
z[,3]=c(0,0,0,0,1,1)
Se=6
Sa=2
V=Sa*z%*%diag(1, 3, 3)%*%t(z)+Se*diag(1, 6, 6)
VI=solve(V)

bL=solve(t(x)%*%VI%*%x)%*%t(x)%*%VI%*%y
bS=Sa*diag(1, 3, 3)%*%t(z)%*%VI%*%(y-x%*%bL)

#MME
R=diag(Se, 6, 6)
RI=solve(R)
c11=t(x)%*%RI%*%x
c12=t(x)%*%RI%*%z
c21=t(c12)
c22=solve(Sa*diag(1, 3, 3))+t(z)%*%RI%*%z
MME_m=rbind(cbind(c11, c12), cbind(c21, c22))
MME_y=rbind(t(x)%*%RI%*%y, t(z)%*%RI%*%y)

MME_b=solve(MME_m)%*%MME_y
