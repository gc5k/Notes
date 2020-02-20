source("shotgun.R")
M=100
N=100
frq=runif(M, 0.05, 0.5)
LewontinLD=runif(N, -1, 1) #LewontinLD is between -1 to 1.
LD=CalLD(frq, LewontinLD) #convert LewontinLD to LD parameter 
geno=GenerateGeno(frq, LD, N)
plot(1-colMeans(geno)/2, frq)

#under special case when frq=0.5
