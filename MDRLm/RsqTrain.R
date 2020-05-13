simu=10000
N=c(100, 200, 500, 1000)
rsq=c(0.05, 0.1) #you can realized rsq with your own linear model
M=c(10, 20, 50)


PCIM=array(0, dim=c(length(N), length(rsq), length(M)))
PCIM2=array(0, dim=c(length(N), length(rsq), length(M)))

layout(matrix(1:4, 2, 2))
for(i in 1:length(N)) {
  for(j in 1:length(rsq)) {
    for(k in 1:length(M)) {
      topP=pchisq(rchisq(simu, df = 1, ncp = N[i]*rsq[j]), 1, lower.tail = F)
      #      hist(topP)
      cb=length(combn(M[k],2))
      null_topP=rbeta(simu, 1, cb)
      
      #proportion of correted identifed model
      PCIM[i,j,k]=length(which(topP< 0.05/cb))/simu
      PCIM2[i,j,k]=length(which(topP<null_topP))/simu
      print(paste(N[i], rsq[j], M[k], PCIM[i,j,k]))
    }
  }
  barplot(PCIM2[i,,], beside = T)
}
