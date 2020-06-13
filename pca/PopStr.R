library(Rcpp)
sourceCpp("~/git/Notes/R/RLib/Shotgun.cpp")
#you may find this file at
#https://github.com/gc5k/Notes/blob/master/R/RLib/Shotgun.cpp


N=200
M=100
frq=runif(M, 0.1, 0.9) #frequency
dprime=runif(M-1, -1, 1) #normalized correlation
# dprime=runif(M-1, 0,0) # no ld
 dprime=runif(M-1, 0.5, 0.5) # normalized r=0.5
Dprime2Cor=Dprime2CorRcpp(frq, dprime) #convert Dprime to pearson's correlation
print(paste(frq[1], frq[2], dprime[1], Dprime2Cor[1])) #check

m=GenerateGenoDprimeRcpp(frq, dprime, N) #generate data
cm=cor(m)
corL=cm[row(cm)==(col(cm)-1)]
plot(corL, Dprime2Cor)

