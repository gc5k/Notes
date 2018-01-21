args = commandArgs(T)
if(length(args) >=6 )
{
  N1=as.numeric(args[1])
  N2=as.numeric(args[2])
  N12=as.numeric(args[3])
  h2=as.numeric(args[4])
  Mk=as.numeric(args[5])
  Q=as.numeric(args[6])  
} else {
  N1=100
  N2=100
  N12=50
  h2=0.01
  Mk=100
  Q=10
}

library("plotrix")
Exp=(sqrt(N1*N2)/Mk * h2 + N12/sqrt(N1*N2))/(sqrt(N1/Mk*h2+1)*sqrt(N2/Mk*h2+1))

if(Sys.info()[['sysname']] == "Darwin")
{
  source("~/R/MyLib/shotgun.R")
} else {
  source("~/bin/MyLib/shotgun.R")
}

beta=matrix(0, Mk, 1)
beta[sample(Mk, Q),1]=rnorm(Q)

freq=runif(Mk, 0.05, 0.5)
ld=runif(Mk, -1, 1)
ld=rep(0, Mk)

REP=16
CR=matrix(0, REP, 1)

png(width=1600, height=1600, filename=paste0("N1=", N1, "_N2=", N2, "_N12", N12, "_h2=", h2, "_Mk=", Mk, "_Q=", Q, ".png"))
layout(matrix(1:REP, 4, 4))
par(ps=20)
par(mai=c(0.8, 0.8, 0.6, 0.6))

##################
for(r in 1:REP)
{
  #generating data
  g1=GenerateGeno(freq,ld, N1)
  g2=GenerateGeno(freq,ld, N2)
  gv1=g1 %*% beta
  gv2=g2 %*% beta
  e1=sqrt(var(gv1)/h2 * (1-h2))
  E1=rnorm(N1, 0, e1)
  E2=rnorm(N2, 0, e1)
  y1=gv1+E1
  y2=gv2+E2

  #setting overlapping samples
  if(N12 > 0)
  {
    g1[1:N12,] = g2[1:N12,]
    y1[1:N12,1] =y2[1:N12,1]
  }
  
  ###########gwas1
  Z=matrix(Mk, Mk, 2)
  for(i in 1:Mk)
  {
    md1=lm(y1~g1[,i])
    md2=lm(y2~g2[,i])
    Z[i,1]=summary(md1)$coefficients[2,3]
    Z[i,2]=summary(md2)$coefficients[2,3]
  }
  plot(Z[,1], Z[,2], xlab="Z1", ylab="Z2")
  abline(a=0, b=Exp, col="red")
  draw.ellipse(x=mean(x1), y=mean(x2), a=min(max(x1), -1*min(x1)), b=min(max(x2), -1*min(x2))*(1-Exp), border = 'red', lwd = 1, angle=Exp*45)
  CR[r,1] = cor(Z[,1], Z[,2])
}
dev.off()
print(args)
print(paste("Exp rho=", Exp))
print(mean(CR[,1]))
print(sd(CR[,1]))
print(CR)

r=0.8
x1=rnorm(1000)
x2=r*x1+sqrt(1-r^2)*rnorm(1000)
plot(x1, x2)
abline(lm(x1~x2))
draw.ellipse(x=mean(x1), y=mean(x2), a=min(max(x1), -1*min(x1)), b=min(max(x2), -1*min(x2))*(1-r), border = 'red', lwd = 1, angle=tan(1/(1+r^2))/pi*180)

