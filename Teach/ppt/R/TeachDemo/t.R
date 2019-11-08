#sample pop
mt=matrix(c(0, 1, 1, 2:7), 3, 3)
layout(mt)
pop=rnorm(1000)
hist(pop, xlim=c(-8,8), freq = T, col="grey")
curve(dnorm(x)*500, seq(-5, 5, length=  1000), xlim=c(-8,8), lwd=2, col="red", bty='l', add=TRUE)
for(i in 1:6) {
  p1=sample(pop, 100)
  plot(density(p1), bty='l', main=paste("Sample", i))
#  hist(p1, xlim=c(-8,8), col="grey", add=T)
#  curve(dnorm(x)*50, p1, add=TRUE, col="red")
}

#inteval
sm=array(200)
for(i in 1:200) {
  sm[i]=mean(sample(pop, 100))
}
hist(main="Mean of samples", sm)
abline(v=c(-1.96/10, 1.96/10), col="red")
