u1_old=1
u1_new=1.05

sigma_old=1
sigma_new=1

sigma=sigma_old

k=1

alpha=0.05
beta=0.2

z_a=qnorm(1-alpha/2)
z_b=qnorm(1-beta)

n1=(1+1/k)*sigma^2*(z_a+z_b)^2/(u1_old-u1_new)^2
n2=k*n1

simu=1000
s_diff=array(0, dim=c(simu, 3))
for(i in 1:simu) {
  s1=rnorm(n1, u1_old, sigma_old)
  s2=rnorm(n2, u1_new, sigma_new)
  s_diff[i,1]=mean(s2)-mean(s1)
  tt=t.test(s1, s2)
  s_diff[i,2]=tt$p.value
  s_diff[i,3]=tt$statistic
}
hist(s_diff[,1])
pw=length(which(s_diff[,2]< alpha))/simu
print(paste("Power is [Diff]: " , pw))
