AUCmax <- function (auc, k)
{
  t = qnorm(1-k)
  i = dnorm(t)/k
  v = -i * k /(1-k)
  
  Q=pnorm(auc)
  
  hl= 2 * Q^2 / ( (v-i)^2 + Q^2 * i * (i-T) + v * (v-T) )
  hl
}
