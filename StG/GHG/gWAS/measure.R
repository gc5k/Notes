library(caTools)

args=commandArgs(TRUE)
geno=args[1]
pheno=args[2]
trt=args[3]
rep=as.numeric(args[4])
type=args[5]

COR=matrix(0, rep, 2)
B=matrix(0, rep, 2)
AUC=matrix(0, rep, 2)

Str=unlist(strsplit(geno, "/"))
root=Str[length(Str)]

PHE=read.table(pheno, as.is=T)
TRT=read.table(trt, as.is=T)
for(i in 1:rep)
{
  TestIdx=which(TRT[,(3+i)] == 0)
  ##make phe
  dat=PHE[TestIdx, c(1,2,3)]    

  ##make profile
  profileQ=read.table(paste0(root, ".", type, ".", i, ".profile"), as.is=T, header=T)
  profileQ=profileQ[TestIdx,]
  Q_r=cor(dat[,3], profileQ$SCORE.1)
  Q_b=cov(dat[,3], profileQ$SCORE.1)/var(profileQ$SCORE.1)
  Q_auc=colAUC(profileQ$SCORE.1, dat[,3], plotROC=FALSE)[1]

  COR[i,1] = Q_r

  B[i,1] = Q_b
  
  AUC[i,1] = Q_auc
  print(paste(Q_r))
  print(paste(Q_b))
  print(paste(Q_auc))
}
print("========")
print(paste(mean(COR[,1]), sd(COR[,1])/sqrt(rep), sd(COR[,1])))
print(paste(mean(B[,1]), sd(B[,1])/sqrt(rep), sd(B[,1])))
print(paste(mean(AUC[,1]), sd(AUC[,1])/sqrt(rep), sd(AUC[,1])))
