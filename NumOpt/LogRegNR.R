#Example from
#https://statacumen.com/teach/SC1/SC1_11_LogisticRegression.pdf
f.lr.p <- function(X, beta) {
# compute vector p of probabilities for logistic regression with logit link
X <- as.matrix(X)
beta <- as.vector(beta)
p <- exp(X %*% beta) / (1 + exp(X %*% beta))
return(p)
}

f.lr.l <- function(y, m, p) {
# binomial log likelihood function
# input: vectors: y = counts; m = sample sizes; p = probabilities
# output: log-likelihood l, a scalar
l <- t(y) %*% log(p) + t(m - y) %*% log(1 - p)
return(l)
}



f.lr.FS <- function(X, y, m, beta.1
                    , eps1 = 1e-6, eps2 = 1e-7, maxit = 50) {
# Fisher's scoring routine for estimation of LR model (with line search)
# Input:
# X = n-by-(r+1) design matrix
# y = n-by-1 vector of success counts
# m = n-by-1 vector of sample sizes
# beta.1 = (r+1)-by-1 vector of starting values for regression est
# Iteration controlled by:
# eps1 = absolute convergence criterion for beta
# eps2 = absolute convergence criterion for log-likelihood
# maxit = maximum allowable number of iterations
# Output:
# out = list containing:
# beta.MLE = beta MLE
# NR.hist = iteration history of convergence differences
# beta.hist = iteration history of beta
# beta.cov = beta covariance matrix (inverse Fisher's information matrix at MLE)
# note = convergence note
beta.2 <- rep(-Inf, length(beta.1)) # init beta.2
diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
llike.1 <- f.lr.l(y, m, f.lr.p(X, beta.1)) # update loglikelihood
llike.2 <- f.lr.l(y, m, f.lr.p(X, beta.2)) # update loglikelihood
diff.like <- abs(llike.1 - llike.2) # diff
if (is.nan(diff.like)) { diff.like <- 1e9 }
i <- 1 # initial iteration index
alpha.step <- seq(-1, 2, by = 0.1)[-11] # line search step sizes, excluding 0
NR.hist <- data.frame(i, diff.beta, diff.like, llike.1, step.size = 1) # iteration history
beta.hist <- matrix(beta.1, nrow = 1)
while ((i <= maxit) & (diff.beta > eps1) & (diff.like > eps2)) {
  i <- i + 1 # increment iteration
  # update beta
  beta.2 <- beta.1 # old guess is current guess
  mu.2 <- m * f.lr.p(X, beta.2) # m * p is mean
  # variance matrix
  v.2 <- diag(as.vector(m * f.lr.p(X, beta.2) * (1 - f.lr.p(X, beta.2))))
  score.2 <- t(X) %*% (y - mu.2) # score function
  # this increment version inverts the information matrix
  # Iinv.2 <- solve(t(X) %*% v.2 %*% X) # Inverse information matrix
  # increm <- Iinv.2 %*% score.2 # increment, solve() is inverse
  # this increment version solves for (beta.2-beta.1) without inverting Information
  increm <- solve(t(X) %*% v.2 %*% X, score.2) # solve for increment
  # line search for improved step size
  llike.alpha.step <- rep(NA, length(alpha.step)) # init llike for line search
  for (i.alpha.step in 1:length(alpha.step)) {
  llike.alpha.step[i.alpha.step] <- f.lr.l(y, m
                                           , f.lr.p(X, beta.2 + alpha.step[i.alpha.step] * increm))
  }
  # step size index for max increase in log-likelihood (if tie, [1] takes first)
  ind.max.alpha.step <- which(llike.alpha.step == max(llike.alpha.step))[1]
  beta.1 <- beta.2 + alpha.step[ind.max.alpha.step] * increm # update beta
  diff.beta <- sqrt(sum((beta.1 - beta.2)^2)) # Euclidean distance
  llike.2 <- llike.1 # age likelihood value
  llike.1 <- f.lr.l(y, m, f.lr.p(X, beta.1)) # update loglikelihood
  diff.like <- abs(llike.1 - llike.2) # diff
  # iteration history
  NR.hist <- rbind(NR.hist, c(i, diff.beta, diff.like, llike.1, alpha.step[ind.max.alpha.step]))
  beta.hist <- rbind(beta.hist, matrix(beta.1, nrow = 1))
}

# prepare output
out <- list()
out$beta.MLE <- beta.1
out$iter <- i - 1
out$NR.hist <- NR.hist
out$beta.hist <- beta.hist
v.1 <- diag(as.vector(m * f.lr.p(X, beta.1) * (1 - f.lr.p(X, beta.1))))
Iinv.1 <- solve(t(X) %*% v.1 %*% X) # Inverse information matrix
out$beta.cov <- Iinv.1
if (!(diff.beta > eps1) & !(diff.like > eps2)) {
out$note <- paste("Absolute convergence of", eps1, "for betas and"
                  , eps2, "for log-likelihood satisfied")
}
if (i > maxit) {
out$note <- paste("Exceeded max iterations of ", maxit)
}
return(out)
}



#beet <- read.table("http://statacumen.com/teach/SC1/SC1_11_beetles.dat", header = TRUE)
beet=read.table("beet.txt", as.is = T, header = T)
beet$rep <- factor(beet$rep)
# create data variables: m, y, X
n <- nrow(beet)
m <- beet$n
y <- beet$y
X.temp <- beet$conc
# quadratic model
X <- matrix(c(rep(1,n), X.temp, X.temp^2), nrow = n)
colnames(X) <- c("Int", "conc", "conc2")
r <- ncol(X) - 1 # number of regression coefficients - 1
# initial beta vector
beta.1 <- c(log(sum(y) / sum(m - y)), rep(0, r))
# fit betas using our Fisher Scoring function
out <- f.lr.FS(X, y, m, beta.1)
out

X0 <- matrix(c(X.temp, X.temp^2), nrow = n)

glm.beetles2=glm(formula = cbind(y, n - y) ~ X0, family = binomial,data = beet)
