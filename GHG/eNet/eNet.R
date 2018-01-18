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
print(paste("Type:", type))

phe=read.table(pheno, as.is=T)
TRT=read.table(trt, as.is=T)
keep=which(TRT[, trtIdx]==1)
phe=phe[keep, c(1, 2, pheIdx)]

root=unlist(strsplit(geno, "/"))
out=paste0(root[length(root)], ".", type, ".", REP)

write.table(phe, out, row.names=F, col.names=F, quote=F)
keepF=paste0(out, ".keep")
write.table(phe[keep, c(1,2)], keepF, row.names=F, col.names=F, quote=F)

#GWAS
if (type == "Q")
{
  GWAS=paste(plink, " --bfile ", geno, " --allow-no-sex --pheno ", out, " --linear --out ", out, " --noweb")
} else if (type == "B")
{
  GWAS=paste(plink, " --bfile ", geno, " --allow-no-sex --pheno ", out, " --fisher --out ", out, " --noweb")        
}
print(GWAS)
system(GWAS)

#step2 filter
snpF=paste0(out, ".snp")
scoreF=paste0(out, ".score")
if (type == "B")
{
  scoreB=read.table(paste0(out, ".assoc.fisher"), as.is=T, header=T)
  odB=order(-1*log10(scoreB$P), decreasing=T)
  scoreB8000=scoreB[odB[1:8000],]
  m=merge(scoreB, scoreB8000, by.x="SNP", by.y="SNP", sort=F)
  write.table(m[,c(1,4,8)], snpF, col.names=F, row.names=F, quote=F)
} else if(type == "Q")
{
  scoreQ=read.table(paste0(out, ".assoc.linear"), as.is=T, header=T)
  odQ=order(-1*log10(scoreQ$P), decreasing=T)
  scoreQ8000=scoreQ[odQ[1:8000],]
  m=merge(scoreQ, scoreQ8000, by.x="SNP", by.y="SNP", sort=F)
  write.table(m[,c(1,4,9)], snpF, col.names=F, row.names=F, quote=F)
}

#step3 prepare data for elasticnet
exF=paste0(out, ".ex")
extract=paste(plink, "--bfile", geno, "--keep", keepF, "--extract", snpF, "--make-bed --out", exF, "--noweb")
print(extract)
system(extract)

#step4
if (type == "B")
{
  library(MultiPhen)
  datB=read.plink(exF)
  yB=read.table(out, as.is=T)
  famB=read.table(paste0(exF, ".fam"), as.is=T)
  mpB=merge(famB, yB, by.x="V2", by.y="V2", sort=F)
  y_B=mpB$V3.y
  library(glmnet)
  cvobB=cv.glmnet(datB, y_B, family="binomial", type.measure="auc")
  coB=coef(cvobB, s="lambda.1se")
  
  idxB=which(coB[,1]!=0)
  idxB=idxB[-1]
  bimB.ex=read.table(paste0(exF, ".bim"), as.is=T)
  bimB.ex=bimB.ex[idxB-1,]
  
  bimB.ex$eff=coB[idxB,]
  write.table(bimB.ex[,c(2,5,7)], scoreF, row.names=F, col.names=F, quote=F)  
} else if(type == "Q")
{
  library(MultiPhen)
  datQ=read.plink(exF)
  yQ=read.table(out, as.is=T)
  famQ=read.table(paste0(exF, ".fam"), as.is=T)
  mpQ=merge(famQ, yQ, by.x="V2", by.y="V2", sort=F)
  y_Q=mpQ$V3.y
  library(glmnet)
  cvobQ=cv.glmnet(datQ, y_Q, family="gaussian", type.measure="deviance")
  coQ=coef(cvobQ, s="lambda.1se")
  
  idxQ=which(coQ[,1]!=0)
  idxQ=idxQ[-1]
  bimQ.ex=read.table(paste0(exF, ".bim"), as.is=T)
  bimQ.ex=bimQ.ex[idxQ-1,]

  bimQ.ex$eff=coQ[idxQ,]
  write.table(bimQ.ex[,c(2,5,7)], scoreF, row.names=F, col.names=F, quote=F)
}

#step5 prediction
pre=paste(gear, "profile --bfile", geno, "--score", scoreF, "--no-score-header --no-weight --keep-atgc --out", out)  
print(pre)
system(pre)
