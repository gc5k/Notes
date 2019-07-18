gear='java -jar /Users/gc5k/Documents/workspace/FromSVN/GEAR/gear.jar'
gear='java -jar ~/bin/gear.jar'

n=5000
m=2
hsq=0.5
rep=100

simu=30
result=matrix(0, simu, 4)
FQ=matrix(0, simu, 4)
for(s in 1:simu) {
  Dprime=runif(1, -1, 1)
  simuPop=paste(gear, "simuqt --n ", n,
                " --m ", m,
                " --unif-freq --hsq ", hsq,
                " --poly-effect --rep", rep,
                " --make-bed --out test",
                " --ld", Dprime,
                " --seed ", s*100)
  system(simuPop)

  locus=paste(gear, "locus --bfile test --out test")
  system(locus)

  grm=paste(gear, "wgrm --bfile test --out test")
  system(grm)

  wgrm=paste(gear, "wgrm --vandem --bfile test --out test_w")
  system(wgrm)

  for(i in 1:rep) {
    he=paste0(gear, " he --grm test --pheno test.phe --mpheno ", i,
              " --out t", i)
    system(he)
    heW=paste0(gear, " he --grm test_w --pheno test.phe --mpheno ", i, 
               " --out wt", i)
    system(heW)
  }

  system("grep Beta1 t*.he > t.txt")
  system("grep Beta1 wt*.he > wt.txt")

  TT=read.table("t.txt", as.is = T)
  WT=read.table("wt.txt", as.is = T)

###########
  frqF=read.table("test.locus", as.is = T, header = T)
  fq=frqF$Freq
  
  W=diag(2*fq*(1-fq), length(fq), length(fq))
  U=diag(1, length(fq), length(fq))

  
  w=matrix(2*fq*(1-fq), 1, length(fq))
  u=matrix(1,1,length(fq))

  I=diag(1, length(fq), length(fq))
  H=sqrt(W)

  bF=read.table("test.rnd", as.is = T)
  beta=matrix(bF$V3, length(fq), 1)

  maf=array(0, dim=length(fq))
  for(k in 1:length(fq)) {
    if(fq[k] < 0.5) {
      maf[k]=fq[k]
    } else {
      maf[k]=1-fq[k]
    }
  }

  rhoCap = min(maf[1]*(1-maf[2]), maf[2]*(1-maf[1]))

  FQ[s,1] = fq[1]
  FQ[s,2] = fq[2]
  FQ[s,3] = Dprime
  rhoD = rhoCap*Dprime
  FQ[s,4] = rhoD/sqrt(maf[1]*(1-maf[1])*maf[2]*(1-maf[2]))

  V = matrix(0, 2, 2)
  V[,1] = c(1, FQ[s,4])
  V[,2] = c(FQ[s,4],1)

  Pm = matrix(c(1, FQ[s,4]^2, FQ[s,4]^2, 1), length(fq), length(fq))

  La = array(0, dim=c(2,2,2))
  LaW = array(0, dim=c(2,2,2))
  for(v in 1:ncol(V)) {
    La[v,,] = u[1,v]*V[,v]%*% t(V[,v])
    LaW[v,,] =w[1,v]*V[,v]%*% t(V[,v])
  }

  h2=sum(diag(U))*t(beta)%*%H%*%(La[1,,]+La[2,,])%*%H%*%beta/(u%*%Pm%*%t(u))
  h2w=sum(diag(W))*t(beta)%*%H%*%(LaW[1,,]+LaW[2,,])%*%H%*%beta/(w%*%Pm%*%t(w))

  result[s,1] = mean(TT$V2)/-2
  result[s,2] = mean(WT$V2)/-2
  result[s,3] = h2
  result[s,4] = h2w
}
write.table(result, "res.txt", row.names = F, col.names = F, quote = F)
write.table(FQ, "FQ.txt", row.names=F, col.names=F, quote=F)

res=read.table("res.txt", as.is = T)

layout(matrix(1:4, 2, 2))
plot(main="No weight", res[,3], res[,1], bty="l", col=ifelse(res[,1]<res[,3], "red", "blue"), 
     xlab=expression(paste(h^2, " (inner mass)")), ylab=expression(paste(hat(h^2)," (observed)")), pch=16)
abline(a=0, b=1, lty=2)

plot(main="Weight", res[,4], res[,2], bty="l", col=ifelse(res[,2]<res[,4], "red", "blue"),
     xlab=expression(paste(h^2, " (inner mass)")), ylab=expression(paste(hat(h^2)," (observed)")), pch=16)
abline(a=0, b=1, lty=2)

plot(main="Theory", res[,1], res[,2], bty="l", col=ifelse(res[,1]<res[,2], "red", "blue"), xlab="Observed", ylab="Theory", pch=16)
abline(a=0, b=1, lty=2)

plot(main="Oberved", res[,3], res[,4], bty="l", col=ifelse(res[,3]<res[,4], "red", "blue"), xlab="No weight", ylab="Weight", pch=16)
abline(a=0, b=1, lty=2)
