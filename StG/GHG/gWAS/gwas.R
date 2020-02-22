rsub='Rscript /clusterdata/gc5k/bin/rsub.R'
plink='/clusterdata/gc5k/bin/plink-1.07-x86_64/plink'
gear='java -Xmx10G -jar /clusterdata/gc5k/bin/gear.jar'

args=commandArgs(TRUE)
geno=args[1]
pheno=args[2]
pheIdx=as.numeric(args[3])
trt=args[4]
trtIdx=as.numeric(args[5])
REP=args[6]
type=args[7]

print(paste("Geno:", geno))
print(paste("Pheno:", pheno))
print(paste("Pheno index:", pheIdx))
print(paste("Trt:", trt))
print(paste("Trt index:", trtIdx))
print(paste("REP:", REP))
print(paste("type:", type))

phe=read.table(pheno, as.is=T)
TRT=read.table(trt, as.is=T)
keep=which(TRT[, trtIdx]==1)
phe=phe[keep, c(1, 2, pheIdx)]

root=unlist(strsplit(geno, "/"))
out=paste0(root[length(root)], ".", type, ".", REP)

write.table(phe, out, row.names=F, col.names=F, quote=F)

#GWAS
if(type == "Q")
{
  PLINKQ=paste(plink, " --bfile ", geno, " --allow-no-sex --pheno ", out, " --linear --out ", out, " --noweb")
} else if (type == "B")
{
  PLINKQ=paste(plink, " --bfile ", geno, " --allow-no-sex --pheno ", out, " --fisher --out ", out, " --noweb")        
}
print(PLINKQ)
system(PLINKQ)

#extract result
scoreF=paste0(out, ".score")
if(type == "B")
{
  scoreQ=read.table(paste0(out, ".assoc.fisher"), as.is=T, header=T)
  write.table(scoreQ[,c(2,4,9)], scoreF, row.names=F, col.names=F, quote=F)
  preQ=paste(gear, " profile --bfile ", geno, " --logit --score ", scoreF, " --no-score-header --no-weight --keep-atgc --out ", out)
} else if(type == "Q" )
{
  scoreQ=read.table(paste0(out, ".assoc.linear"), as.is=T, header=T)
  write.table(scoreQ[,c(2,4,7)], scoreF, row.names=F, col.names=F, quote=F)
  preQ=paste(gear, " profile --bfile ", geno, " --score ", scoreF, " --no-score-header --no-weight --keep-atgc --out ", out)   
}

print(preQ)
system(preQ)
