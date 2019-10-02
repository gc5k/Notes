#Example from
#http://faculty.washington.edu/eliezg/teaching/StatR503/CourseMaterials/Week9/BayesianMCMC.html

prior <- function(theta) dunif(theta)
likelihood <- function(theta)
  ifelse(theta >= 0 & theta <=1, theta^5, 0)
theta <- 0.3
n <- 1000
p.old <- prior(theta)*likelihood(theta)
thetas <- theta
while(length(thetas) <= n){
  theta.new <- theta + rnorm(1,0,0.05)
  p.new <- prior(theta.new)*likelihood(theta.new) #posterior
  if(p.new > p.old | runif(1) < p.new/p.old){
    theta <- theta.new
    thetas <- c(thetas,theta)
    p.old <- p.new
  }
}

plot(thetas, type="o", pch=21, cex=0.5, bg="darkgrey")
hist(thetas[-(1:100)], freq=FALSE, col="grey", bor="darkgrey")
curve(6*x^5, add=TRUE, col=2, lwd=3)
