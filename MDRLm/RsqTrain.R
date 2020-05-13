simu=10000
N=c(100, 200, 500)
rsq=c(0.05, 0.1) #you can realized rsq with your own linear model
M=c(10, 20, 50)

PCIM=array(0, dim=c(length(N), length(rsq), length(M)))
for(i in 1:length(N)) {
  for(j in 1:length(rsq)) {
    for(k in 1:length(M)) {
      topP=pchisq(rchisq(simu, df = 1, ncp = N[i]*rsq[j]), 1, lower.tail = F)
      hist(topP)
      cb=length(combn(M[k],2))
      #proportion of correted identifed model
      PCIM[i,j,k]=length(which(topP< 0.05/cb))/simu
      print(paste(N[i], rsq[j], M[k], PCIM[i,j,k]))
    }
  }
}
barplot(PCIM[i,,], beside = T)
