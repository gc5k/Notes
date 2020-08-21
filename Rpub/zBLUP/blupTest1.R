#Robinson's example (1991)

y=c(110, 100, 110, 100, 100, 110, 110, 100, 100)
x=matrix(0, 9, 3)
x[,1]=c(1,1,0,0,0,0,0,0,0)
x[,2]=c(0,0,1,1,1,0,0,0,0)
x[,3]=c(0,0,0,0,0,1,1,1,1)

z=matrix(0, 9, 4)
z[,1]=c(1,0,0,0,0,0,0,0,0)
z[,2]=c(0,0,1,0,0,0,0,0,0)
z[,3]=c(0,0,0,0,0,1,1,0,0)
z[,4]=c(0,1,0,1,1,0,0,1,1)

R=diag(1, 9, 9)
RI=solve(R)
G=diag(0.1, 4, 4)
GI=solve(G)

c11=t(x)%*%RI%*%x
c12=t(x)%*%RI%*%z
c21=t(c12)
c22=t(z)%*%RI%*%z+GI

xy=t(x)%*%RI%*%y
zy=t(z)%*%RI%*%y
MME_y=rbind(xy, zy)

MME_M=rbind(cbind(c11, c12), cbind(c21, c22))
b=solve(MME_M)%*%MME_y

VI=solve(diag(0.1, 9, 9)+R)
BL=solve(t(x)%*%VI%*%x)%*%t(x)%*%VI%*%y
bs=0.1*t(z)%*%VI%*%(y-x%*%BL)
plot(b, c(BL, bs))
