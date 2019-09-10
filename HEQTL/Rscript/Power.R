Mk=1000000
N=1000
h2=0.01
rsq=0.9
alpha=0.05/Mk
Power=0.5


falPos=1
power=0
while ( (falPos > alpha) || (power <Power) )
{
  N=N+1
  NCP=N*h2*rsq
  power=pchisq(NCP, 1, ncp=h2*rsq)
  falPos=pchisq(NCP, 1, lower.tail=F)
  print(paste("N=", N, "NCP=", NCP, " power=", power, " type I=", falPos))

}
