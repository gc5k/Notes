#https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

SamplingMethod = 1 # 0 for gibbs; 1 for mcmc

beta=c(5, 0, 1)
sd=1
sampleSize=200
# create independent x-values 
X=matrix(1, sampleSize, length(beta))
for(i in 2:length(beta)) {
  X[,i]=rbinom(sampleSize, 2, 0.5)
}

# create dependent values according to ax + b + N(0,sd)
y <-  X%*%beta+rnorm(sampleSize, 0, sd)
XtX=t(X)%*%X
XtX.i=solve(XtX)
B=XtX.i%*%t(X)%*%y

par(mfrow=c(1,3))
plot(X[,1], y, main="Test Data")
plot(X[,2], y, main="Test Data")
plot(X[,3], y, main="Test Data")

log.likelihood <- function(param){
  Beta = param[1:(length(param)-1)]
  sd0 = param[length(param)]

  pred = X%*%Beta
  singlelikelihoods = dnorm(y, mean = pred, sd = sd0, log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

# Example: plot the likelihood profile of the slope a

# Prior distribution
log.prior <- function(param){
  return(sum(dnorm(param, sd=1, log=T)))
}

# Posterior distribution
log.posterior <- function(param){
  return (log.likelihood(param) + log.prior(param))
}

######## Metropolis algorithm ################

proposalfunction <- function(param) {
  
  return(c(rnorm(length(param)-1, mean = param, sd= 0.3), runif(1, 0.1, 3)))
}

library(MASS)
proposalfunctionGibbs <- function(param) {
  b=mvrnorm(1, B, XtX.i)

  rs=y-X%*%B
  sg=1/rgamma(1, length(y)/2, t(rs)%*%rs/2)
  return(c(b, sg))
}

run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    if(SamplingMethod == 1) {
      proposal = proposalfunction(chain[i,])
    } else {
      proposal = proposalfunctionGibbs(chain[i,])
    }
    probab = exp(log.posterior(proposal) - log.posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

startvalue = c(4,0,1,3)
chain = run_metropolis_MCMC(startvalue, 50000)

burnIn = 2000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

### Summary: #######################

par(mfrow = c(2,length(beta)+1))
hist(chain[-(1:burnIn),1],nclass=30, , main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]), col="green")
abline(v = beta[1], col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b1", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]), col="green")
abline(v = beta[2], col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of b2", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]), col="green")
abline(v = beta[3], col="red" )
hist(chain[-(1:burnIn),4],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),4]), col="green")
abline(v = sd, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = beta[1], col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = beta[2], col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = beta[3], col="red" )
plot(chain[-(1:burnIn),4], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = beta[4], col="red" )

# for comparison:
summary(lm(y~X-1))