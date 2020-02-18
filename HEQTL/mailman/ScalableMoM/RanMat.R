REP=1000
vz4=matrix(0, REP, 1)
for(i in 1:REP) {
  z1=rnorm(5000)
  z2=rnorm(5000)
  z3=z1*z2
  z4=z3^2
  vz4[i,1]=var(z4)
}
hist(vz4)
