source("~/R/MyLib/shotgun.R")
gear='java -jar /Users/gc5k/Documents/workspace/FromSVN/GEAR/gear.jar'
REP=100
LOCI=10
effMat=matrix(0, LOCI, 1)
effMat[,1]=rep(1, LOCI)
h2Mat=matrix(0, REP, LOCI)
write.table(effMat, "effect.txt", row.names = F, col.names = F, quote=F, append=F)
result=matrix(0, REP, LOCI)

LD = c(-0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75)
REC = c(0.01, 0.1, 0.25, 0.5)

LD = rep(0.75, LOCI-1)
REC = rep(0.1, LOCI-1)
LD=c(0.75)
REC=c(0.1)

mat=matrix(0, length(REC), length(LD))
smat=matrix(0, length(REC), length(LD))
mat1=matrix(0, length(REC), length(LD))
smat1=matrix(0, length(REC), length(LD))

for(r in 1:length(REC))
{
  for(ld in 1:length(LD))
  {
    for(rep in 1:REP)
    {
      cmd=paste(gear, "simufam --num-fam 1000 --num-marker ", LOCI , " --effect-file effect.txt --out test --hsq 0.5 --unif-freq")
      cmd=paste(cmd, "--ld ", LD[ld], " --rec ", REC[r])
      cmd=paste(cmd, "--make-bed --seed ", rep+1000)
      print(cmd)
      system(cmd)
      fam=read.table("test.fam", as.is=T)
      idx=which(fam$V3 == 0)
      write.table(fam[-idx, c(1,2)], "keep.txt", row.names = F, col.names = F, quote = F)

      he=paste0(gear, " hefam --keep keep.txt --ibd test --pheno test.phe --out test", rep)
      print(he)
      system(he)
    }
    system("grep Beta1 *.he > he.txt")
    he=read.table("he.txt", as.is = T)
    mat[r, ld] = mean(he$V2)
    smat[r, ld] = sd(he$V2)

    system("grep Mean *.he > heM.txt")
    heM=read.table("heM.txt", as.is = T)

    mat1[r, ld] = mean(heM$V2)
    smat1[r, ld] = sd(heM$V2)
  }
}
write.table(mat, "heBeta.txt", row.names = F, col.names = F, quote = F)
write.table(smat, "heBetas.txt", row.names = F, col.names = F, quote = F)

write.table(mat1, "heBetaM.txt", row.names = F, col.names = F, quote = F)
write.table(smat1, "heBetaMs.txt", row.names = F, col.names = F, quote = F)
