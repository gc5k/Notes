##########################
#D1-1
# sample size
##########################
N=1000
M=20000
#uniform
freq=runif(M)
#0.5
freq=rep(0.5, M)
#

gmat=matrix(0, N, M)
for(i in 1:M)
{
  gmat[,i] = rbinom(N, 2, freq[i])
}
layout(matrix(1:4, 2, 2))
freqDat=colMeans(gmat)/2
plot(freqDat, frame.plot = F, xlab="Locus", ylab="Freq", pch=16)
qqplot(freqDat, rnorm(M, freq, sqrt(freq*(1-freq)/(2*N))), xlab = "Observed", ylab="Expected", frame.plot=F, pch=16)
abline(a=0, b=1, col="red")
pvN=pnorm((freqDat-0.5)/sqrt(freqDat*(1-freqDat)/(2*N)))
pvB=pbinom(apply(gmat, 2, sum), 2*N, 0.5)
hist(pvN, main="Normal distribution")
hist(pvB, main="Binomial distribution")

hdtest=matrix(0, M, 2)
for(i in 1:M)
{
  gcount=table(gmat[,i])
  hdtest[i,1] = (gcount[1] - N/4)^2/(N/4) + (gcount[2] - N/2)^2/(N/2) + (gcount[3] - N/4)^2/(N/4)
}
layout(matrix(1:2, 1, 2))
plot(density(hdtest[,1]), main="density of chisq[2df]")
hist(hdtest[,2], main="Hardy p-value")

###############################
#D1-2
#Maize real inbred lines
###############################
layout(matrix(1:6, 2, 3))
FinaC=read.table("Fina_Genotype.ped", as.is = T, na.strings = "99")[,-c(1:6)]
GCnt=array(0, ncol(FinaC))
for(i in 1:ncol(FinaC))
{
  GCnt[i] = length(which(!is.na(FinaC[,i])))
}

FinaSC = matrix(NA, nrow = nrow(FinaC), ncol = ncol(FinaC))
for(i in 1:ncol(FinaC))
{
  idx=which(!is.na(FinaC[,i]))
  FinaSC[idx,i]=rbinom(length(idx), 1, 0.5)*2
}

FinaFreq=colMeans(FinaC, na.rm = T)/2
FinaSCFreq=colMeans(FinaSC, na.rm = T)/2
plot(FinaFreq, col="grey", bty="l")
qqplot(FinaFreq, FinaSCFreq, xlab = "Observed freq", ylab="Expected freq", bty="l")
abline(a=0, b=1, col="red")

#########t-test
#FinaFreq = FinaSCFreq
tTest=matrix(0, ncol(FinaC), 2)
tTest[,1]=(FinaFreq-0.5)/sqrt(FinaFreq*(1-FinaFreq)/(GCnt))
tTest[,2]=2*pt(-1*abs(tTest[,1]), GCnt-1)
hist(tTest[,2], breaks = 50, main="t test")

#######chi-sq
#FinaC=FinaSC
chiTest=matrix(0, ncol(FinaC), 2)
for(i in 1:ncol(FinaC))
{
  gcnt=table(FinaC[,i])
  gm=mean(gcnt)
  chiTest[i,1] = (gcnt[1]-gm)^2/gm + (gcnt[2]-gm)^2/gm
  chiTest[i,2] = pchisq(chiTest[i,1], 1, lower.tail = F)
}
hist(chiTest[,2], breaks = 50, main="Chisq test")

plot(FinaFreq, col=ifelse(chiTest[,2] < (0.05/length(FinaFreq)), "red", "grey"), pch=16, bty="l")

####################
#D1-3 Gold Standard
#p-value distribution
###################
layout(matrix(1:4, 2, 2))
N=c(500, 500, 500)

tS=rt(N[1], df=2)
chiS=rchisq(N[2], df=2)
normS=rnorm(N[3], mean=20)
hist(tS, xlim=c(-10, 30), main="3 distribution")
hist(chiS, add=T, col="red")
hist(normS, add=T, col="blue")
hist(c(pt(tS, df = 2), pchisq(chiS, df=2), pnorm(normS, mean=20)), main="p-value", xlab="P-value")

##########p-value, different df, mu
N=c(500, 500, 500)
dft=runif(N[1], 2, 10)
tS=rt(N[1], df=dft)

dfchi=runif(N[2], 1, 8)
chiS=rchisq(N[2], df=dfchi)

muN=runif(N[3], 10, 20)
normS=rnorm(N[3], mean=muN)
hist(tS, xlim=c(-10, 30), main="3 distribution (various df, mu)")
hist(chiS, add=T, col="red")
hist(normS, add=T, col="blue")
hist(c(pt(tS, df = dft), pchisq(chiS, df=dfchi), pnorm(normS, mean=muN)), main="p-value", xlab="P-value")
