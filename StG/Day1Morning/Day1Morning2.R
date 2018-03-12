##########################
#D1.2-1
# Shoot
##########################
layout(matrix(1:4, 2,2))
REP=c(1,2,4,256)
n=100
for(i in 1:length(REP))
{
  plot(main=paste(REP[i], "shoot"), x=NULL, y=NULL, xlim=c(-4,4), ylim=c(-4,4), axes = F, xlab="", ylab="")
  for(j in 1:n)
  {
    points(mean(rnorm(REP[i])), mean(rnorm(REP[i])), pch=16, cex=0.5, col="red")
    abline(h=0, v=0)
  }
  points(rep(0, 6),rep(0,6), cex=seq(5,30, 5))
}


##########################
#D1.2-2
#QQplot 1
##########################
layout(matrix(1:4, 2,2))
REP=c(1,2,4,256)
n=1000
for(i in 1:length(REP))
{
  xp=array(0, n)
  for(j in 1:n)
  {
    xp[j] = mean(rnorm(REP[i]))
  }
  qqplot(xp, rnorm(n, 0, sqrt(1/REP[i])), bty='l', xlab="Simulated", ylab="Theoretical")
  abline(a=0, b=1, col="red")
}

##########################
#D1.2-3
#QQplot 2
##########################
layout(matrix(1:4, 2,2))
REP=c(1,2,4,256)
n=1000
for(i in 1:length(REP))
{
  xp=array(0, n)
  for(j in 1:n)
  {
    xp[j] = mean(rnorm(REP[i]))
  }
  qqplot(xp, rnorm(n), bty='l', xlab="Simulated", ylab="Theoretical")
  abline(a=0, b=1, col="red")
}

