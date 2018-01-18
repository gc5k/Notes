N=100 #sample size
p1=0.4 #freq 1 (the proportion of cases)
p2=0.7 #freq 2 (the proportion of high-risk genotypes)

d1=rbinom(N, 1, p1) #simulate case-control
d2=rbinom(N, 1, p2) #simulate risk-cells
x=table(d1, d2) #contigency table

#as defined in wiki
a=x[1,1]
b=x[1,2]
c=x[2,1]
d=x[2,2]

ab=a+b
cd=c+d
ac=a+c
bd=b+d
n=sum(x)

num=log10(factorial(ab))+log10(factorial(cd))+log10(factorial(ac))+log10(factorial(bd))
dem=log10(factorial(a))+log10(factorial(b))+log10(factorial(c)) +log10(factorial(d)) +log10(factorial(n))
pobs=num-dem

#rearrange the contigency table 
if(ab > cd)
{
  tt=ab
  ab=cd
  cd=tt
}
if(ac > bd)
{
  tt=ac
  ac=bd
  bd=tt
}

Ptable=matrix(0:min(ab, ac), nrow=min(ab, ac)+1, ncol=2) #full probability table
pv=0 #pvalue
for(i in 0:min(ab, ac))
{
  a=i
  b=ab-a
  c=ac-a
  d=bd-b
  num=log10(factorial(ab)) +log10(factorial(cd)) +log10(factorial(ac)) + log10(factorial(bd))
  dem=log10(factorial(a)) +log10(factorial(b)) +log10(factorial(c)) + log10(factorial(d)) + log10(factorial(n))
  p=num-dem
  Ptable[i+1,2]=10^p

  if(10^p <= 10^pobs)
  {
    pv=pv+10^p
  }
}
plot(Ptable[,1], Ptable[,2], xlab="Observation", ylab="Probability")
print(paste("Pvalue is: ", pv))

fisher.test(x)
