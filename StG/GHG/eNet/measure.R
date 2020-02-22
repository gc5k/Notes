library(caTools)

args=commandArgs(TRUE)
root=args[1]
src=paste0("../", root, ".phe.glm.rep")

rep=as.numeric(args[2])
COR=matrix(0, rep, 2)
B=matrix(0, rep, 2)
AUC=matrix(0, rep, 2)
dat=read.table(src, as.is=T)
for(i in 1:rep)
{
  idx=which(dat[,(4+i)] == 0)
  ##make phe
  datQ=dat[idx,c(1,2,4)]
  datB=dat[idx,c(1,2,3)]
  
  ##make profile
  profileQ=read.table(paste0(root, ".Tr", i, ".Q.profile"), as.is=T, header=T)
  profileQ=profileQ[idx,]
  Q_r=cor(datQ[,3], profileQ$SCORE.1)
  Q_b=cov(datQ[,3], profileQ$SCORE.1)/var(profileQ$SCORE.1)
  Q_auc=colAUC(profileQ$SCORE.1, datB[,3], plotROC=FALSE)[1]

  profileB=read.table(paste0(root, ".Tr", i, ".B.profile"), as.is=T, header=T)
  profileB=profileB[idx,]
  B_r=cor(datB[,3], profileB$SCORE.1)
  B_b=cov(datB[,3], profileB$SCORE.1)/var(profileB$SCORE.1)
  B_auc=colAUC(profileB$SCORE.1, datB[,3], plotROC=FALSE)[1]
  COR[i,1] = Q_r
  COR[i,2] = B_r

  B[i,1] = Q_b
  B[i,2] = B_b
  
  AUC[i,1] = Q_auc
  AUC[i,2] = B_auc
  print(paste(Q_r, B_r))
  print(paste(Q_b, B_b))
  print(paste(Q_auc, B_auc))
}
print("========")
print(paste(mean(COR[,1]), sd(COR[,1])/sqrt(rep)))
print(paste(mean(COR[,2]), sd(COR[,2])/sqrt(rep)))
print(paste(mean(B[,1]), sd(B[,1])/sqrt(rep)))
print(paste(mean(B[,2]), sd(B[,2])/sqrt(rep)))
print(paste(mean(AUC[,1]), sd(AUC[,1])/sqrt(rep)))
print(paste(mean(AUC[,2]), sd(AUC[,2])/sqrt(rep)))
