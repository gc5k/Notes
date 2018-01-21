plink='/clusterdata/gc5k/bin/plink-1.07-x86_64/plink'
gear='java -Xmx20G -jar /clusterdata/gc5k/bin/gear.jar'

args=commandArgs(TRUE)
root=args[1]
print(root)
src=paste0("../", root)
IDX=as.numeric(args[2])

#step0 keep
phe=read.table(paste0("../", root, ".phe.glm.rep"))
keepInd=which(phe[,4+IDX]==1)
keep=paste0(root, ".Tr", IDX, ".keep")
write.table(phe[keepInd,c(1,2)], keep, row.names=F, col.names=F, quote=F)

pheQ=paste0(root, ".Tr", IDX, ".Q")
pheB=paste0(root, ".Tr", IDX, ".B")
write.table(phe[keepInd, c(1,2,3)], pheB, row.names=F, col.names=F, quote=F)
write.table(phe[keepInd, c(1,2,4)], pheQ, row.names=F, col.names=F, quote=F)

#step1 GWAS
PLINKQ=paste(plink, "--bfile", src, "--pheno", pheQ, "--linear --out", pheQ, "--noweb")
print(PLINKQ)
system(PLINKQ)

PLINKB=paste(plink, "--bfile", src, "--pheno", pheB, " --1 --fisher --out", pheB, "--noweb")
print(PLINKB)
system(PLINKB)

#step2 filter
scoreQ=read.table(paste0(pheQ, ".assoc.linear"), as.is=T, header=T)
idxQ=which(scoreQ$TEST=="ADD")
scoreQ=scoreQ[idxQ,]
odQ=order(-1*log10(scoreQ$P), decreasing=T)
scoreQ8000=scoreQ[odQ[1:8000],]
m=merge(scoreQ, scoreQ8000, by.x="SNP", by.y="SNP", sort=F)
write.table(m[,c(1,4,9)], paste0(pheQ,".snp"), col.names=F, row.names=F, quote=F)

scoreB=read.table(paste0(pheB, ".assoc.fisher"), as.is=T, header=T)
#idxB=which(scoreB$TEST=="ADD")
#scoreB=scoreB[idxB,]
odB=order(-1*log10(scoreB$P), decreasing=T)
scoreB8000=scoreB[odB[1:8000],]
m=merge(scoreB, scoreB8000, by.x="SNP", by.y="SNP", sort=F)
write.table(m[,c(1,4,8)], paste0(pheB, ".snp"), col.names=F, row.names=F, quote=F)

#step3 prepare data for elasticnet
extractQ=paste(plink, "--bfile", src, "--keep", keep, "--extract", paste0(pheQ, ".snp"), "--make-bed --out", paste0(pheQ, ".extract"), "--noweb")
print(extractQ)
system(extractQ)

imputationQ=paste(gear, " --bfile", paste0(pheQ, ".extract"), "--naive-imputation --make-bed --out", paste0(pheQ, ".extract.imput"))
print(imputationQ)
system(imputationQ)

extractB=paste(plink, "--bfile", src, "--keep", keep, "--extract", paste0(pheB, ".snp"), "--make-bed --out", paste0(pheB, ".extract"), "--noweb")
print(extractB)
system(extractB)

imputationB=paste(gear, " --bfile", paste0(pheB, ".extract"), "--naive-imputation --make-bed --out", paste0(pheB, ".extract.imput"))
print(imputationB)
system(imputationB)

#step4
library(MultiPhen)
datQ=read.plink(paste0(pheQ,".extract.imput"))
yQ=read.table(pheQ, as.is=T)
famQ=read.table(paste0(pheQ, ".extract.imput.fam"), as.is=T)
mpQ=merge(famQ, yQ, by.x="V2", by.y="V2", sort=F)
y_Q=mpQ$V3.y
library(glmnet)
cvobQ=cv.glmnet(datQ, y_Q, family="gaussian", type.measure="deviance")
coQ=coef(cvobQ, s="lambda.1se")

idxQ=which(coQ[,1]!=0)
idxQ=idxQ[-1]
bimQ.ex=read.table(paste0(pheQ, ".extract.imput.bim"), as.is=T)
bimQ.ex=bimQ.ex[idxQ-1,]

bimQ.ex$eff=coQ[idxQ,]
write.table(bimQ.ex[,c(2,5,7)], paste0(pheQ, ".score"), row.names=F, col.names=F, quote=F)

library(MultiPhen)
datB=read.plink(paste0(pheB, ".extract.imput"))
yB=read.table(pheB, as.is=T)
famB=read.table(paste0(pheB, ".extract.imput.fam"), as.is=T)
mpB=merge(famB, yB, by.x="V2", by.y="V2", sort=F)
y_B=mpB$V3.y
library(glmnet)
cvobB=cv.glmnet(datB, y_B, family="binomial", type.measure="auc")
coB=coef(cvobB, s="lambda.1se")

idxB=which(coB[,1]!=0)
idxB=idxB[-1]
bimB.ex=read.table(paste0(pheB, ".extract.imput.bim"), as.is=T)
bimB.ex=bimB.ex[idxB-1,]

bimB.ex$eff=coB[idxB,]
write.table(bimB.ex[,c(2,5,7)], paste0(pheB, ".score"), row.names=F, col.names=F, quote=F)

#step5 prediction
preQ=paste(gear, "profile --bfile", src, "--score", paste0(pheQ, ".score"), "--no-score-header --no-weight --keep-atgc --out", pheQ)
print(preQ)
system(preQ)

preB=paste(gear, "profile --bfile", src, "--score", paste0(pheB, ".score"), "--no-score-header --no-weight --keep-atgc --out", pheB)
print(preB)
system(preB)
