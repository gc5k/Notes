rep=100
b=matrix(0, rep, 1)
r=0.3
for(i in 1:rep)
{
  z1=rnorm(5000)
  z2=z1*r+sqrt(1-r^2)*rnorm(5000)
  cor(z1, z2)
  cov(z1, z2)
  c1=z1^2
  c2=z2^2
  mod=lm(c1~c2)
  b[i,1]=mod$coefficient[2]
}
plot(b[,1])
abline(h=r^2)
