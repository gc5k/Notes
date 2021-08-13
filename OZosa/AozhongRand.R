seed=20210812 #2021812
set.seed(seed)
n=60
Unit=20
k=n/Unit

tSeq=matrix(0, n, 4)
for(i in 1:Unit) {
  tSeq[((i-1)*k+1):(i*k),1]=i
  tSeq[((i-1)*k+1):(i*k),2]=seq(1,k)
}

colnames(tSeq)=c("Unit", "index", "A", "B")
for(i in 1:n) {
  tSeq[i,3:4]=sample(c(1,0),2)
}

write.csv(tSeq, file="AzRand.csv", row.names = T)
