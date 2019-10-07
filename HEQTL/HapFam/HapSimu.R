m = 15
n = 100
nNew = 100
gmat = matrix(rbinom(m*n, 1, 0.5), n, m)
gd = matrix(0.1, 1, m-1)

for(k in 1:nNew) {
  s = nrow(gmat)

  gIdx=ceiling(runif(1, 0, nrow(gmat)))
  gVec=matrix(0, 1, m)
  gVec[1,1]=gmat[gIdx, 1]

  for(l in 1:(m-1)) {
    r1 = (1-exp(-gd[l]/s))/s
    r2 = exp(-gd[l]/s)
    gS=ceiling(runif(1, 0, nrow(gmat)))
    if(gVec[l] == gmat[gS,l+1]) {
      r=r1+r2
    } else {
      r=r1
    }
    
    r0=runif(1, 0, 1)
    if(r < r0) {
      gIdx = gS
    }
    gVec[l+1]=gmat[gIdx,l+1]
  }
  ss=1/sum(seq(1, nrow(gmat)))
  theta=ss/(2*(n+k-1 + ss))
  gmat=rbind(gmat, gVec)
}
