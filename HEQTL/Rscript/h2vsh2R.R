B[1,]=c(1,1)
B[2,]=c(-1,1)

B[3,]=c(1,2)
B[4,]=c(-2,1)

B[5,]=c(1,3)
B[6,]=c(-3,1)
rho=seq(-1, 1, 0.01)
mat=matrix(1:6, 2, 3)
layout(mat)
for(i in 1:nrow(B)) {
  b=B[i,]
  vy=b[1]^2+b[2]^2+2*rho*b[1]*b[2]
  hsnp=0.5*(b[1]^2+b[2]^2)+2*rho*b[1]*b[2]/(1+rho^2)
  h2_f=hsnp/vy
  plot(rho, h2_f, ylim=c(0, 0.6), main=b, pch=16, xlab=expression(rho), ylab=expression(h[SNP]^2))
  abline(h=0.5)
}


#############

B[1,]=c(1, 1)
B[2,]=c(1, -1)
B[3,]=c(1, 3)
B[4,]=c(-1, 3)
B[5,]=c(3, 1)
B[6,]=c(-3, 1)

mat=matrix(1:6, 2, 3)
layout(mat)
for(i in 1:nrow(B)) {
  b=B[i,]
  vg=1/2*b[1]^2+3/8*b[2]^2+sqrt(3)/2*RHo*b[1]*b[2]
  vy=2*vg

  hsnp=(1/2*b[1]^2 + 3/8*b[2]^2)+sqrt(3)*RHo*b[1]*b[2]/(1+RHo^2)
  h2_nw=hsnp/vy
  plot(RHo, h2_nw, xlab=expression(rho), ylab=expression(h[snp]^2), cex=0.5, ylim=c(0, 0.7), main=b)

  hsnpW=(sqrt(0.5)+sqrt(3/8))*((sqrt(1/2)*(3/8*RHo^2*b[2]^2 +sqrt(3)/2 * RHo *b[1]*b[2]+0.5*b[1]^2))
                               + (sqrt(3/8)*(1/2*RHo^2*b[1]^2 +sqrt(3)/2 * RHo *b[1]*b[2] +3/8*b[2]^2)))/(7/8+sqrt(3)/2*RHo^2)
  h2_w=hsnpW/vy
  points(RHo, h2_w, col="red", cex=0.5)
  abline(h=c(0.5))
}
