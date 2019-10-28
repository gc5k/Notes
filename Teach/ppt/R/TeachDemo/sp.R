pop=1:100
rp=100
k=c(4,25,49)
mV=matrix(0, rp, length(k))
for(i in 1:rp) {
  for(j in 1:length(k)) {
    s=sample(pop,k[j])
    mV[i,j]=mean(s)
  }  
}
par(mfrow=c(2,3))
plot(main=paste("K=", k[1]), (mV[,1]), ylim=c(1, 100), ylab="mean", pch=16, col="red")
plot(main=paste("K=", k[2]),(mV[,2]), ylim=c(1, 100), ylab="mean", pch=16, col="gray")
plot(main=paste("K=", k[3]),(mV[,3]), ylim=c(1, 100), ylab="mean", pch=16, col="blue")

#plot(sort(mV[,1]), ylim=c(1, 100), pch=16, col="red")
#plot(sort(mV[,2]), ylim=c(1, 100), pch=16, col="gray")
#plot(sort(mV[,3]), ylim=c(1, 100), pch=16, col="blue")

qqplot(xlim=c(0, 100), ylim=c(0, 100), rnorm(rp, 50, sd(pop)/sqrt(k[1])), sort(mV[,1]), pch=16, col="red")
abline(a=0, b=1)
qqplot(xlim=c(0, 100), ylim=c(0, 100), rnorm(rp, 50, sd(pop)/sqrt(k[2])), sort(mV[,2]), pch=16, col="gray")
abline(a=0, b=1)
qqplot(xlim=c(0, 100), ylim=c(0, 100), rnorm(rp, 50, sd(pop)/sqrt(k[3])), sort(mV[,2]), pch=16, col="blue")
abline(a=0, b=1)


##coin
k=c(10, 50, 100)
rp=1000
coinMean=matrix(0, rp, length(k))
for(i in 1:length(k)) {
  for(j in 1:rp) {
    s=rbinom(1, k[i], 0.5)
    coinMean[j,i]=s
  }
}
par(mfrow=c(1,3))
hist(coinMean[,1], xlim=c(0, k[1]), breaks = 25, main=k[1], col="red")
abline(v=mean(coinMean[,1]), lwd=2, lty=2)
hist(coinMean[,2], xlim=c(0, k[2]), breaks = 25, main=k[2], col="green")
abline(v=mean(coinMean[,2]), lwd=2, lty=2)
hist(coinMean[,3], xlim=c(0, k[3]), breaks = 25, main=k[3], col="blue")
abline(v=mean(coinMean[,3]), lwd=2, lty=2)

sort(coinMean[,1])[c(25,975)]
sort(coinMean[,2])[c(25,975)]
sort(coinMean[,3])[c(25,975)]
