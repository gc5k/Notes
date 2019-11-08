zb=qnorm(0.8)
za=qnorm(0.975)
k=1
ua=0.3
ub=0
sa=1.1
sb=1.5

na=(zb*sb+za*sa*sqrt(k))^2/(k*(ua-ub)^2)
nb=k*na
na=(sa^2+sb^2/k)*(za+zb)^2/(ua-ub)^2
nb=(k*sa^2+sb^2)*(za+zb)^2/(ua-ub)^2


pv=array(0, 100)
for(i in 1:100) {
  sp1=rnorm(na, ua, sa)
  sp2=rnorm(nb, ub, sb)
  
  
  tt=t.test(sp1, sp2)
  pv[i]=tt$p.value  
}
length(which(pv < 0.05))