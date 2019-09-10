REP=1
K=10
N=1000
TM=matrix(0, REP, K)
PM=matrix(0, REP, K)
LM=array(0, REP)
for(rep in 1:REP) {
  GM=matrix(0, K, K)
  GM[K,]=1:K
  gN=K
  for(i in K:2) {
    TM[rep,i]=rexp(1,  rate = (i*(i-1)/2)/N)
#    PM[rep,i]=qexp(TM[rep,i], ncol(combn(gN,2)))
    r2=sort(sample(which(GM[i,]>0), 2))
#    print(r2)
    GM[i-1,]=GM[i,]
    GM[i-1, r2[2]]=-GM[i, r2[1]]
    gN=gN-1
  }
  print(sum(PM[rep,]))
  LM[rep]= sum(PM[rep,])
}
