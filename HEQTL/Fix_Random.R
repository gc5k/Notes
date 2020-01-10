N=500 #sample size 
fMin=0.1 #freq min
fMax=0.5 #freq max

h2F=0.1 #h2 for fixed
QF=3 #loci number of fixed

h2R=0.4 #h2 for random
QR=100 #loci number of random

#simulate fixed 
fqF=runif(QF, fMin, fMax)
gF=matrix(0, N, QF)
for(i in 1:QF) gF[,i]=rbinom(N, 2, fqF[i])
bF=rnorm(QF)
gbF=gF%*%bF

###########
fqR=runif(QR, fMin, fMax)
gR=matrix(0, N, QR)
for(i in 1:QR) gR[,i]=rbinom(N, 2, fqR[i])

sVR=var(gbF)/h2F*h2R
sv=sum(apply(gR, 2, var))

s2=sqrt(sVR/sv)
bR=rnorm(QR, sd = s2)
gbR=gR%*%bR
vR=var(gbR)

vG=var(gbF+gbR)
se=sqrt(vG/(h2F+h2R)*(1-h2F-h2R))

y=gbF+gbR+rnorm(N,sd=se)
G=cbind(gF, gR)
print(paste0("variance fixed: ", var(gbF)))
print(paste0("variance random: ", var(gbR)))
print(paste0("variance total: ", var(y)))
