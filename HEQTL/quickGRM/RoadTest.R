library(Rcpp)
sourceCpp("cormatrix.cpp")

sp=c(15)#, 20, 25)
sigma=c(0, 1)
s=length(sigma)
FD=10
for(i in 1:length(sp)) {
  n1=sp[i]
  m=s^sp
  loci=n1*s^n1*FD
  G=matrix(rbinom(loci, max(sigma), 0.5), n1*FD, m) #it generate a matrix
  Gt=t(G)
  Tstart1=proc.time()
  for(j in 1:FD) { #generate off-diagonal elements of the relationship matrix using mailman
    g1=G[((j-1)*n1+1):(j*n1),]
    g2=Gt[,1:(j*n1)]
    Pc=MailmanProductL(g1, g2, sigma)
  }
  Tend1=proc.time()
  print(Tend1[3]-Tstart1[3])
}

for(i in 1:length(sp)) {
  n1=sp[i]
  m=s^sp
  loci=n1*s^n1*FD
  G=matrix(rbinom(loci, max(sigma), 0.5), n1*FD, m)
  Gt=t(G)
  Tstart2=proc.time()
  for(j in 1:FD) { #us common matrix multiplication
    g1=G[((j-1)*n1+1):(j*n1),]
    g2=Gt[,1:(j*n1)]
    Pc=GRM(g1, g2)
  }
  Tend2=proc.time()
  print(Tend2[3]-Tstart2[3])
}


Tstart3=proc.time()
gg=G%*%Gt #buildin function
Tend3=proc.time()
print(Tend3[3]-Tstart3[3])
