#this example is adapted from http://ljwolf.org/teaching/gibbs/

library(MASS)
library(MCMCpack)
Beta = c(1,2)
sd = 1
N = 100
X = cbind(rnorm(N), rnorm(N))
y = 0.5 + X[,1]*1 + X[,2]*2 + rnorm(N)
n_obs = length(y)

X = cbind(rep(1, n_obs), X) # gonna include a constant
XtX = t(X) %*% X
n_params = 3 # two covariates and a constant

beta_hat = solve(XtX, t(X) %*% y) # compute this ahead of time cause we'll need it a lot
XtXi = solve(XtX)

beta = c(0,0,0) # starting value
sigma2 = 0 #starting value
n_iterations = 5000

beta_out = matrix(data=NA, nrow=n_iterations, ncol=n_params)
sigma_out = matrix(data = NA, nrow = n_iterations, ncol=1)
for (i in 1:n_iterations){
  beta = mvrnorm(n=1, beta_hat, sigma2 * XtXi) # need the multivariate cause we have three beta
  
  part = (y - X %*% beta)
  
  sigma2 = 1/rgamma(1, n_obs/2, t(part) %*% part * .5)
  #or sigma2 = rinvgamma(1, n_obs/2, t(part) %*% part * .5 ) from library (MCMCpark)
  
  beta_out[i,] = beta
  sigma_out[i,] = sigma2
}