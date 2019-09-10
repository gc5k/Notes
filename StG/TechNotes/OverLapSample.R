R=100
V_delta_y=array(0, R)
N=1000
M=1000

for(r in 1:R)
{
  y=rnorm(1000)
  y_H=matrix(0, N, M)
  for(i in 1:M)
  {
    x=rnorm(1000)
    mod=lm(y~x)
    y_H[,i] = x*mod$coefficient[2]
  }
  y_pre=apply(y_H, 1, sum)
  
  delta_y=y-y_pre
  V_delta_y[r]=var(delta_y)
}
prediction=(1-M/N)^2+M/N
plot(V_delta_y)

abline(a=prediction, b=0, col="red")

lambda=seq(0.01, 5, length.out=100)
plot(lambda, (1-lambda)^2+lambda)
