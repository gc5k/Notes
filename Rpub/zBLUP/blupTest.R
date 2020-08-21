Ml=20
Ms=100
N=1000

Xl=apply(matrix(rbinom(Ml*N, 2, 0.5), N, Ml), 2, scale)
Gl=Xl%*%t(Xl)/Ml
Xs=apply(matrix(rbinom(Ms*N, 2, 0.5), N, Ms), 2, scale)
Gs=Xs%*%t(Xs)/Ms

hs=0.5
hsQ=hs/(2*0.5*0.5*Ms)
he=1
bl=rnorm(Ml, 1)
bs=rnorm(Ms, 0, sqrt(hs/(2*0.5*0.5*Ms)))

Fb=Xl%*%bl
Sb=Xs%*%bs
Ev=rnorm(N, 0, sqrt(he))
y=Fb+Sb+Ev

V=Gs*hs+diag(he, N, N)
V_Inv=solve(V)

layout(matrix(1:6, 2, 3, byrow = T))
bl_est=solve(t(Xl)%*%V_Inv%*%Xl)%*%t(Xl)%*%V_Inv%*%y
plot(bl, bl_est, pch=16, cex=1.5, col="green")
abline(a=0, b=1, col="red")
ypre=y-Xl%*%bl_est

bs_est=hs/Ms*t(Xs)%*%V_Inv%*%ypre
plot(bs, bs_est, pch=16, cex=1.5, col="green")
abline(a=0, b=1, col="red")

gs_est=hs/Ms*Xs%*%t(Xs)%*%V_Inv%*%ypre
plot(Sb, gs_est, pch=16, cex=0.5,col="green")
abline(a=0, b=1, col="red")

####MME
Ehe=var(Ev)
c11=t(Xl)%*%Xl
c12=t(Xl)%*%Xs
c21=t(c12)
c22_I=t(Xs)%*%Xs
c22=c22_I+solve(diag(1, Ms, Ms))*he/(hs/Ms)

Y1=t(Xl)%*%y
Y2=t(Xs)%*%y

M1=cbind(c11, c12)
M2=cbind(c21, c22)

MMe_mat=rbind(cbind(c11, c12), cbind(c21, c22))
MMe_y=rbind(Y1, Y2)
bMMe=solve(MMe_mat)%*%MMe_y

plot(bl, bMMe[1:Ml], pch=16, col="blue")
abline(a=0, b=1, col="red")
plot(bs, bMMe[(Ml+1):(Ms+Ml)], col="blue")
abline(a=0, b=1, col="red")
plot(bs_est, bMMe[(Ml+1):(Ms+Ml)], col="blue")
abline(a=0, b=1, col="red")
