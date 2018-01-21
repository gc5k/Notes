source("EpiNOIA.R")
N=400
REP=5000
M=500
freq=runif(M, 0.3, 0.5)

F2=matrix(0, N, M)
F2D=matrix(0, N, M)
Ha=matrix(0, N, M)
Hd=matrix(0, N, M)
for(i in 1:M)
{
  F2[,i]= rbinom(N, 2, freq[i])
  F2D[,i] = ifelse(F2[,i] == 1, 1, 0)
  Ha[,i] = NOIA_A(F2[,i])
  Hd[,i] = NOIA_D(F2[,i])
}

A=Ha %*% t(Ha)
GA=N*A/sum(diag(A))

D=Hd %*% t(Hd)
GD=N*D/sum(diag(D))

#######AD epi
##Vitezica
GAD_=GA*GD
GAD=N*GAD_/sum(diag(GAD_))

#CGB
GAD_t=matrix(0, nrow=nrow(GAD_), ncol=ncol(GAD_))
for(i in 1:nrow(GAD_t))
{
  for(j in 1:ncol(GAD_t))
  {
    for(k in 1:M)
    {
      GAD_t[i,j] =GAD_t[i,j] + Ha[i,k] * Ha[j,k] * Hd[i,k] * Hd[j,k]
    }
  }
}
gad=A*D-GAD_t
GAD_cgb=N*gad/sum(diag(gad))

####AA epi,
##Vitezica
GAA_=GA*GA
GAA=N*GAA_/sum(diag(GAA_))

##CGB
GAA_t=matrix(0, nrow=nrow(GAA_), ncol=ncol(GAA_))
for(i in 1:nrow(GAA_t))
{
  for(j in 1:i)
  {
    for(k in 1:M)
    {
      GAA_t[i,j] =GAA_t[i,j] + Ha[i,k] * Ha[j,k] * Ha[i,k] * Ha[j,k]
    }
    GAA_t[j,i] =GAA_t[i,j]
  }
}
gaa=A*A-GAA_t
GAA_cgb=N*gaa/sum(diag(gaa))

layout(matrix(1:4, 2,2))
hist(GAD_cgb[row(GAD_cgb) > col(GAD_cgb)], breaks = 50, main="AD")
hist(GAD[row(GAD) > col(GAD)], breaks = 50, main="AD (Vitezica)")

hist(GAA_cgb[row(GAA_cgb) > col(GAA_cgb)], breaks = 50, main="AA")
hist(GAA[row(GAA) > col(GAA)], breaks = 50, main="AA (Vitezica)")



#######HE
REP=100
VC=matrix(0, REP, 2)
for(rep in 1:REP)
{
  B=rnorm(M, 0, sqrt(1/(sum(2*freq*(1-freq)))))
  D=rnorm(M, 0, sqrt(1/(sum(4*freq^2*(1-freq)^2))))
  y=scale(F2 %*% B+F2D %*% D + rnorm(N, 0, 1))
  
  Y=array(0, dim=N*(N-1)/2)
  heX=array(0, dim=N*(N-1)/2)
  heXD=array(0, dim=N*(N-1)/2)
  cnt=1
  for(i in 2:nrow(y))
  {
    for(j in 1:(i-1))
    {
      Y[cnt]=(y[i]-y[j])^2
      heX[cnt]=GA[i,j]
      heXD[cnt]=GD[i,j]
      cnt=cnt+1
    }
  }
  
  md=lm(Y~heX+heXD)
  VC[rep, 1] = -md$coefficients[2]/md$coefficients[1]
  VC[rep, 2] = -md$coefficients[3]/md$coefficients[1]
}
