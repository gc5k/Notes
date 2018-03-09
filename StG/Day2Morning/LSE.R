N=400
REP=5000

#######wrong model
layout(matrix(1:4, 2,2))
beta=matrix(0, REP, 5)
for(i in 1:REP)
{
  f2=rbinom(N, 2, 0.5)
  dcode=ifelse(f2 == 0, -0.5, ifelse(f2==1, 0.5, -0.5))
  dcode=dcode
  y=matrix(rnorm(N,1), N, 1)
  
  md=lm(y[,1]~dcode - 1) #no mean
  beta[i,1:4] = summary(md)$coefficients[1,]
  beta[i, 5] = summary(md)$fstatistic[1]
}
hist(beta[,1])
hist(beta[,4])

#####
beta=matrix(0, REP, 5)
for(i in 1:REP)
{
  f2=rbinom(N, 2, 0.5)
  dcode=ifelse(f2 == 0, -0.5, ifelse(f2==1, 0.5, -0.5))
  dcode=dcode+0.5
  y=matrix(rnorm(N,1), N, 1)
  
  md=lm(y[,1]~dcode-1)
  beta[i,1:4] = summary(md)$coefficients[1,]
  beta[i, 5] = summary(md)$fstatistic[1]
}
hist(beta[,1])
hist(beta[,4])

##########right model
layout(matrix(1:4, 2,2))
beta=matrix(0, REP, 5)
for(i in 1:REP)
{
  f2=rbinom(N, 2, 0.5)
  dcode=ifelse(f2 == 0, -0.5, ifelse(f2==1, 0.5, -0.5))
  dcode=dcode
  y=matrix(rnorm(N,1), N, 1)
  
  md=lm(y[,1]~dcode)
  beta[i,1:4] = summary(md)$coefficients[2,]
  beta[i, 5] = summary(md)$fstatistic[1]
}
hist(beta[,1])
hist(beta[,4])

beta=matrix(0, REP, 5)
for(i in 1:REP)
{
  f2=rbinom(N, 2, 0.5)
  dcode=ifelse(f2 == 0, -0.5, ifelse(f2==1, 0.5, -0.5))
  dcode=dcode+0.5
  y=matrix(rnorm(N,1), N, 1)
  
  md=lm(y[,1]~dcode)
  beta[i,1:4] = summary(md)$coefficients[2,]
  beta[i, 5] = summary(md)$fstatistic[1]
}
hist(beta[,1])
hist(beta[,4])
