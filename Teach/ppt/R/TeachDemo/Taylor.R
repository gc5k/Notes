x=seq(0,20,length=100)
alpha=qchisq(0.95,1)

curve(dchisq(x,df=1), xlim=c(0, 40), ylab="density", lwd=3)
abline(v=alpha, col="gray", lty=2)
NCP=c(1, 5, 10, 20)
for(i in 1:length(NCP)) {
  curve(dchisq(x,df=1, ncp=NCP[i]), add = T, col=NCP[i])
  pchisq(alpha, df=1, ncp=1, lower.tail = F)
  print(pchisq(alpha, 1, ncp=NCP[i], lower.tail = F))
}
