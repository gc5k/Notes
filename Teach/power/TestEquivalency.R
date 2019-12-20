delta=0.1
u1_old=1
u1_new=1

sigma_old=1
sigma_new=1

sigma=sigma_old

k=1

alpha=0.05
beta=0.2

z_a=qnorm(1-alpha/2) #note two tails
z_b=qnorm(1-beta/2) #note

n1=(1+1/k)*sigma^2*(z_a+z_b)^2/(delta-abs(u1_old-u1_new))^2
n2=k*n1

simu=1000
s_diff=array(0, dim=c(simu, 5))
for(i in 1:simu) {
  s1=rnorm(n1, u1_old, sigma_old)
  s2=rnorm(n2, u1_new, sigma_new)
  s_diff[i,1]=mean(s2)-mean(s1)
  tt=t.test(s1, s2)
  s_diff[i,2]=tt$p.value
  s_diff[i,3]=tt$statistic
  s_diff[i,4]=(mean(s1)-mean(s2)-delta)/(sigma*sqrt(1/n1+1/n2))
  s_diff[i,5]=(mean(s1)-mean(s2)+delta)/(sigma*sqrt(1/n1+1/n2))
}
hist(s_diff[,1])
pw=length(which(abs(s_diff[,1]) > delta))/simu
print(paste("Power is [EB]: " , pw))


RS=10000
SS1=rnorm(RS, mean = u1_new-u1_old, sd=sqrt(1/n1+1/n2))

Sd1=rnorm(RS, mean = -delta, sd=sqrt(1/n1+1/n2))
Sd2=rnorm(RS, mean = delta, sd=sqrt(1/n1+1/n2))

hist(Sd1, col="grey", breaks = 50,xlim=c(-delta*3, delta*3))
hist(Sd2, add=TRUE, col="grey", breaks = 50)
hist(SS1, add=TRUE, breaks = 50, col="blue")

abline(v=0, col="red")
T1=sort(Sd1)[RS*(1-alpha/2)]

abline(v=0, col="red")
T2=sort(Sd2)[RS*alpha/2]
abline(v=c(T1, T2, u1_new-u1_old), lty=2, col="green")

length(which(SS1 > T1 & SS1 < T2))/RS
