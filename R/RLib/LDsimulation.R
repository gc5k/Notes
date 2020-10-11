sourceCpp("Shotgun.cpp")
M=100
N=100
frq=runif(M, 0.4, 0.5)
LewontinLD=runif(M-1, -1, 1) #LewontinLD is between -1 to 1. 
LewontinLD=runif(M-1, 0.9, 0.9) #high ld 

genoDprime=GenerateGenoDprimeRcpp(frq, LewontinLD, N)
plot(1-colMeans(genoDprime)/2, frq)
gCorDprime=cor(genoDprime)

LD=Dprime2LDRcpp(frq, LewontinLD) #convert LewontinLD to LD parameter 
genoLD=GenerateGenoLDRcpp(frq, LD, N)
plot(1-colMeans(genoLD)/2, frq)
gCorLD=cor(genoLD)

#
plot(xlab="Dprime", ylab="LD", gCorDprime[col(gCorDprime)==(row(gCorDprime)-1)]^2,
     gCorLD[col(gCorLD)==(row(gCorLD)-1)]^2)
abline(a=0, b=1, col="red")

rec=0.01 #recombination
gn=100 #generation
LDdecay=LD*(1-rec)^gn #ld decay equation

genoLDdecay=GenerateGenoLDRcpp(frq, LDdecay, N) #a ld decayed population
plot(1-colMeans(genoLDdecay)/2, frq)
gCorLDdecay=cor(genoLDdecay)

ldsq1=LD2CorRcpp(frq, LD) #convert ld to correlation
ldsq2=LD2CorRcpp(frq, LDdecay)

plot(xlim=c(0, 1), ylim=c(0, 1), xlab="original ld", ylab="ld after100 generation", ldsq1^2, ldsq2^2, col="red", pch=16) #original vs decayed
points(gCorDprime[col(gCorDprime)==(row(gCorDprime)-1)]^2, gCorLDdecay[col(gCorLDdecay)==(row(gCorLDdecay)-1)]^2) #simulation
