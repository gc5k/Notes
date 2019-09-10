source("qcFun.R")
library(vioplot)
##read data
proFile="protein_no_normal_impute.txt"
prot=read.table(proFile, as.is = T, header = T)

PPP0=read.csv("PPPB1_lineup.csv", as.is = T, header = T)
PPP1=PPP0
PPP1$PPPA_ID=tolower(PPP1$PPPA_ID)

dat=t(prot[,c(2:ncol(prot))])

misD=matrix(-1, nrow(dat), 4)
for(i in 1:nrow(misD))
{
  naIdx=which(is.na(dat[i,]))
  if(length(naIdx) < ncol(dat) )
  {
    misD[i, 1] = length(naIdx)
    misD[i, 2] = min(dat[i,-naIdx])
    misD[i, 3] = median(dat[i, -naIdx])
    misD[i, 4] = mean(dat[i, -naIdx])
  }
}
blankInd = which(misD[,1] == -1)
pdf("Abudence_Missing.pdf")
#layout(matrix(1:3, 1, 3))
plot(misD[-blankInd,1], misD[-blankInd,2], frame.plot = F, ylab="Abundence (min)", xlab="Missing proteins", pch=16, cex=0.5)
md=lm(misD[-blankInd,2] ~ misD[-blankInd,1])
abline(md, col="red", lwd=2, lty=2)
legend("topright", legend = paste("Rsq=", format(summary(md)$r.squared, digit=3)), bty='n')

# plot(misD[-blankInd,1], misD[-blankInd,3], frame.plot = F, ylab="Abundence (median)", xlab="Missing proteins", pch=16, cex=0.5)
# md=lm(misD[-blankInd,3] ~ misD[-blankInd,1])
# abline(md, col="red", lwd=2, lty=2)
# legend("topright", legend = paste("Rsq=", format(summary(md)$r.squared, digit=3)), bty='n')
# 
# plot(misD[-blankInd,1], misD[-blankInd,4], frame.plot = F, ylab="Abundence (mean)", xlab="Missing proteins", pch=16, cex=0.5)
# md=lm(misD[-blankInd,4] ~ misD[-blankInd,1])
# abline(md, col="red", lwd=2, lty=2)
# legend("topright", legend = paste("Rsq=", format(summary(md)$r.squared, digit=3)), bty='n')

dev.off()

###quantile normalization
###https://en.wikipedia.org/wiki/Quantile_normalization
library(preprocessCore)

qdat=t(normalize.quantiles(t(dat)))
#qdat=normalize.quantiles(dat)

dat=qdat

colnames(dat)=prot$protein_group
rownames(dat)=rownames(t(prot[,c(2:ncol(prot))]))

#line up data
exDat=matrix(0, nrow(dat), ncol(dat))
rID=matrix(0, nrow(dat), 1)
for (i in 1:nrow(PPP1))
{
  idx=which(rownames(dat) == PPP1$PPPA_ID[i])
  rID[i,1]=idx
  exDat[i,] = dat[idx,]
}
rownames(exDat) = PPP1$PPPA_ID[rID[,1]]
colnames(exDat)=colnames(dat)

#write.table(t(exDat), "QN_PPP1.txt", row.names = T, col.names = T, quote = F)
matDat=matrix(NA, nrow = nrow(t(exDat))+1, ncol= ncol(t(exDat)) + 1)
matDat[1,]= colnames(prot)
matDat[c(2:nrow(matDat)),1] = colnames(dat)
matDat[c(2:nrow(matDat)), c(2:ncol(matDat))] = t(exDat)
for(i in 2:nrow(matDat))
{
  naI=which(!is.na(matDat[i, c(2:ncol(matDat))]))
  matDat[i, naI+1] = ceiling(as.numeric(matDat[i, naI+1]))
}
write.table(matDat, "QN_PPP1.txt", row.names = F, col.names = F, quote=F)

#################missing data
ms=matrix(1, nrow = nrow(exDat), ncol=ncol(exDat))
for(i in 1:nrow(ms))
{
  idx=which(is.na(exDat[i,]))
  ms[i,idx]=0
}
TuIdx=which(PPP1$Tissue == "T")
NoIdx=which(PPP1$Tissue == "N")
pdf("Missing_Tissue.pdf")
layout(matrix(1:2, 1, 2))
plot(apply(ms[TuIdx,], 2, sum), apply(ms[NoIdx,], 2, sum), frame.plot = F, xlab = "Observed Tumor samples", ylab="Observed Normal samples", pch=16, cex=0.5)
abline(a=0, b=1, col="red", lwd=2, lty=2)

freq=colMeans(ms[c(TuIdx, NoIdx),])
p1=rbinom(ncol(exDat), length(TuIdx), freq)
p2=rbinom(ncol(exDat), length(NoIdx), freq)
plot(p1, p2, xlab="Expected Observed Tumor samples", ylab="Expected Normal samples", pch=16, cex=0.5, frame.plot = F)
dev.off()

pdf("Protein_Missing.pdf")
colM=colMeans(exDat, na.rm = T)
colMiss=array(0, dim=ncol(exDat))
for(i in 1:ncol(exDat))
{
  colMiss[i] = length(which(is.na(exDat[,i])))
}
layout(matrix(1:2, 1, 2))
idxM=which(colM > 0)
plot(xlab=expression(paste("Protein", log[10], "(TIC)")), ylab="Missing per protein", x=log10(colM[idxM]), y=colMiss[idxM], col="red", pch=16, bty='n', cex=0.3)
plot(xlab=expression(paste("Protein", log[10], "(-1*TIC)")), ylab="Missing per protein", x=log10(-1*colM[-idxM]), y=colMiss[-idxM], col="blue", pch=16, bty='n', cex=0.3)
dev.off()

# pdf("ind_Missing.pdf")
# rowM=rowMeans(exDat, na.rm = T)
# rowMiss=array(0, dim=nrow(exDat))
# for(i in 1:nrow(exDat))
# {
#   rowMiss[i] = length(which(is.na(exDat[i,])))
# }
# idxM=which(!is.na(rowM))
# plot(xlab=expression(paste(log[10], "(ind-TIC)")), ylab="Missing per protein", x=log10(rowM[idxM]), y=rowMiss[idxM], col="red", pch=16, bty='n', cex=0.3)
# dev.off()
# 
# pdf("ind_Missing_MS.pdf")
# layout(matrix(1:6, 2, 3, byrow = T))
# rowM=rowMeans(exDat, na.rm = T)
# rowMiss=array(0, dim=nrow(exDat))
# for(i in 1:nrow(exDat))
# {
#   rowMiss[i] = length(which(is.na(exDat[i,])))
# }
# idxM=which(!is.na(rowM))
# PC=names(table(PPP1$MS_ID))
# for(i in 1:length(PC))
# {
#   plot(main=paste("MS", i), xlim=range(log10(rowM[idxM])), ylim=c(400, 1600), xlab=expression(paste(log[10], "(ind-TIC)")), ylab="Missing per protein", x=log10(rowM[PPP1$MS_ID[idxM] == PC[i]]), y=rowMiss[PPP1$MS_ID[idxM] == PC[i]], col="red", pch=16, bty='n', cex=1)
# }
# dev.off()


jpeg("MissingData.jpeg", width = 800, height = 600)
#pdf("MissingData.pdf")
layout(matrix(c(1,1,2), 1, 3))
plot(x=NULL, y=NULL, xlim=c(1, ncol(ms)), ylim=c(1, nrow(ms)), axes = F, xlab="Protein", ylab = "Sample")
for(i in 1:nrow(ms))
{
  points(x=seq(1, ncol(ms)), y=rep(i, ncol(ms)), col=ms[i,], pch=15, cex=0.1)
}
hist(1-rowMeans(ms), breaks = 50, main="Missing rate per individual", xlab="Missing rate")
dev.off()

####remove samples completely missing
idxMS=which(rowMeans(ms) == 0)
write.table(PPP1[idxMS,], "Blank_ind.txt", row.names=F, col.names=T, quote=F)

exDat=exDat[-idxMS,]
PPP1=PPP1[-idxMS,]

write.table(exDat, "QN_PPP1.txt", row.names = T, col.names = T, quote = F)
write.table(PPP1, "ID_PPP1.txt", row.names = F, col.names = T, quote = F)

####pca analysis
#PCA(exDat, "prot", PPP1$PPPA_ID)

####common 
TSnames=c("AQUA", "CTRL", "Normal", "Tumor")
Aidx=which(PPP1$Tissue == "A")
Cidx=which(PPP1$Tissue == "C")
Nidx=which(PPP1$Tissue == "N")
Tidx=which(PPP1$Tissue == "T")
IDX=c(length(Aidx), length(Cidx), length(Nidx), length(Tidx))

TSidx=list(Aidx)
TSidx[2]=list(Cidx)
TSidx[3]=list(Nidx)
TSidx[4]=list(Tidx)

eve=read.table("prot.eigenvec", as.is = T)

#######grm plot
grm=read.table(gzfile("prot.grm.gz"), as.is = T)
eva=read.table("prot.eigenval", as.is = T)

#jpeg("EigenVal.jpeg", width = 800, height = 600)
pdf("EigenVal.pdf")
layout(matrix(1:2, 1, 2))
hist(grm[,4], breaks=50, main = "Protein relationship matrix", xlab="Score")
barplot(eva[1:20,1], main=paste(nrow(eva), "individuals"), xlab="Eigenvalue", col=ifelse(row(eva) < 6, "red", "grey"), border = "NA")
dev.off()

######################pcA for tissue

#jpeg("PCA_tissue1.jpeg", width = 800, height = 800)
pdf("PCA_tissue1.pdf")
mat=matrix(21, 5, 5)
mat[row(mat) >= col(mat)] = 1:15
mat[1,4:5]=16:17
mat[2,4:5]=18:19
mat[3,5]=20

layout(mat)
par(mai=c(0.3,0.3,0,0))
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i == j)
    {
      plot(x=NULL, y=NULL, xlim=range(eve[,i+2])*1.1, ylim=range(eve[,i+2])*1.1, bty='n', axes = F)
      text(median(eve[,i+2]), median(eve[,i+2]), labels = paste("PC",i))
    } else {
      plot(eve[,i+2], eve[,j+2], xlim=range(eve[,i+2])*1.1, ylim=range(eve[,j+2])*1.1, bty='n', axes = F, col=as.numeric(as.factor(PPP1$Tissue)), pch=16, cex=0.7)
#      points(eve[NAidx, i+2], eve[NAidx, j+2], col="pink", pch=15, cex=1)
      if (i==1)
      {
        axis(side=2)
      }
      if(j == 5)
      {
        axis(side=1)
      }
    }
  }
}

par(mai=c(0,0,0.3,0))

for(i in 1:length(TSidx))
{
  plot(eve[TSidx[[i]],3], eve[TSidx[[i]],4], bty='n', axes = T, xlim=range(eve[,3], na.rm = TRUE)*1.1, ylim=1.1*range(eve[,4], na.rm = TRUE), col=i, pch=16, cex=1)
  points(mean(eve[TSidx[[i]],3]), mean(eve[TSidx[[i]],4]), pch=1, cex=4, lwd=3, col="gold")
  legend("bottomleft", legend = TSnames[i], bty = 'n')
}
dev.off()

#####################################MS vs PC

#jpeg("MZ_tissue.jpeg", width = 800, height = 800)
pdf("MZ_tissue.pdf")
MS=names(table(PPP1$MS_ID))
PTest=array(0, dim=c(5,length(MS),length(MS))) #t-test

layout(matrix(1:30, 5, length(MS), byrow = T))
par(mai=c(0.3, 0.3, 0, 0))

for(i in 1:length(TSidx))
{
  dat=eve[TSidx[[i]],]
  pp=PPP1[TSidx[[i]],]
  for(j in 1:length(MS))
  {
    plot(xlim=range(eve[,3])*1.1, ylim=1.1*range(eve[,4]), dat[pp$MS_ID==MS[j],3], dat[pp$MS_ID==MS[j],4], bty='n', col=i, pch=16)
  }
  
  for(j1 in 2:length(MS))
  {
    for(j2 in 1:(j1-1))
    {
      PTest[i,j1,j2]=t.test(dat[pp$MS_ID==MS[j1], 3], dat[pp$MS_ID==MS[j2], 3])$p.value
    }
  }
}
dev.off()

########################PC for MS
#jpeg("PCA_vs_MS.jpeg", width = 800, height = 800)
pdf("PCA_vs_MS.pdf")
MS=names(table(PPP1$MS_ID))
mat=matrix(22, 5, 5)
mat[row(mat) >= col(mat)] = 1:15
mat[1,4:5]=16:17
mat[2,4:5]=18:19
mat[3,4:5]=20:21

layout(mat)
par(mai=c(0.3,0.3,0,0))

for(i in 1:5)
{
  for(j in i:5)
  {
    if(i == j)
    {
      plot(x=NULL, y=NULL, xlim=c(0, 1), ylim=c(0,1), bty='n', axes = F)
      text(0.1, 0.5, labels = paste("PC",i))
    } else {
      plot(eve[,i+2], eve[,j+2], bty='n', axes = F, col=as.numeric(as.factor(PPP1$MS_ID)), pch=16, cex=1)
      if(i==1)
      {
        axis(side=2)
      }
      if(j == 5)
      {
        axis(side=1)
      }
    }
  }
}

par(mai=c(0.3,0.3,0,0))

for(i in 1:length(MS))
{
  plot(eve[PPP1$MS_ID==MS[i],3], eve[PPP1$MS_ID==MS[i],4], bty='n', axes = T, xlim=1.1*range(eve[,3]), ylim=1.1*range(eve[,4]), col=i, pch=16, cex=1)
  points(mean(eve[PPP1$MS_ID==MS[i],3]), mean(eve[PPP1$MS_ID==MS[i],4]), pch=1, cex=4, col="gold", lwd=3)
  legend("bottomleft", legend = MS[i], bty = 'n')
}

dev.off()


#########mis matrix plot
#jpeg("MisMatrix.jpeg", width = 800, height = 800)
pdf("MisMatrix.pdf")

MisP=array(0, nrow(exDat))
for(i in 1:nrow(exDat))
{
  MisP[i]=length(which(is.na(exDat[i,])))
}

mat=matrix(24, 5, 5)
mat[row(mat) >= col(mat)] = 1:15
mat[1,3:5]=16:18
mat[2,4:5]=19:20
mat[3,4:5]=21:22
mat[4,5]=23
layout(mat)
par(mai=c(0.3,0.3,0,0))

MS=names(table(PPP1$MS_ID))
for(i in 1:5)
{
  for(j in i:5)
  {
    if(i == j)
    {
      plot(x=NULL, y=NULL, xlim=c(0, 1), ylim=c(0,1), bty='n', axes = F)
      text(0.1, 0.5, labels = paste("PC",i))
    } else {
      colM=max(MisP/ncol(exDat))
      colS=min(MisP/ncol(exDat))
      plot(eve[,i+2], eve[,j+2], bty='n', axes = F, col=rgb((MisP/ncol(exDat)-colS)/colM, 1-(MisP/ncol(exDat)-colS)/colM, (MisP/ncol(exDat)-colS)/colM), pch=1, cex=1-MisP/ncol(exDat))# pch=ifelse(MisP/(ncol(exDat)) > 0.5, 15, 2), cex=0.5)
      if(i==1)
      {
        axis(side=2)
      }
      if(j == 5)
      {
        axis(side=1)
      }
    }
  }
}

par(mai=c(0.2,0.2,0.3,0.1))
for(i in 1:length(TSidx))
{
  hist(MisP[TSidx[[i]]], main=TSnames[i], xlim = c(0, 2500))
  abline(v=mean(MisP[TSidx[[i]]], na.rm = T), col="red", lty=2)
}

dev.off()

pdf("violin_tissue.pdf")
vioplot(names=c("AQUA", "CTRL", "Tumor", "Normal"),
        MisP[Aidx], MisP[Cidx], MisP[Tidx], MisP[Nidx])
dev.off()

#jpeg("Miss_vs_eigen.jpeg", width = 800, height = 600)
pdf("Miss_vs_eigen.pdf")
layout(matrix(1:6, 2, 3))
par(mai=c(0.6,0.5,0.3,0.1))
for(i in 1:6)
{
  plot(MisP, eve[,i+2], col="grey", main=paste0("Missing vs E", i), xlab="Missing rate", ylab=paste("PC",i), bty='n', pch=16)
  m1=lm(eve[,i+2]~MisP)
  abline(m1)
  legend("topright", legend = paste("Rsq = ", format(summary(m1)$r.squared, digits=3)), bty='n')
}
dev.off()

pdf("Tissue_vs_eigen.pdf")
layout(matrix(1:6, 2, 3))
par(mai=c(0.6,0.5,0.3,0.1))
TN=which(PPP1$Tissue == "T" | PPP1$Tissue == "N")
for(i in 1:6)
{
  plot(as.numeric(as.factor(PPP1$Tissue[TN])), eve[TN,i+2], main=paste0("Tumor&Normal vs PC", i), xlab="Normal(blue) vs Tumor (red)", ylab=paste("Eigen",i), bty='n', pch=16, col=ifelse(PPP1$Tissue[TN] == "T", "red", "blue"))
  m1=lm(eve[TN,i+2]~PPP1$Tissue[TN])
  abline(m1)
  legend("topright", legend = paste("Rsq = ", format(summary(m1)$r.squared, digits=3)), bty='n')
}
dev.off()



pdf("QQ_MS.pdf")
layout(matrix(1:4, 2, 2))
qqplot(pch=16, xlim=1.1*range(eve[,3]), ylim=1.1*range(eve[,4]), eve[PPP1$MS_ID=="ProCan1",3],  eve[PPP1$MS_ID=="ProCan4", 3], xlab="MS 1", ylab="MS 4", bty='n')
abline(a=0,b=1, lwd=3, lty=2, col="grey")
qqplot(pch=16, xlim=1.1*range(eve[,3]), ylim=1.1*range(eve[,4]), eve[PPP1$MS_ID=="ProCan2",3],  eve[PPP1$MS_ID=="ProCan3", 3], xlab="MS 2", ylab="MS 3", bty='n')
abline(a=0,b=1, lwd=3, lty=2, col="grey")
qqplot(pch=16, xlim=1.1*range(eve[,3]), ylim=1.1*range(eve[,4]), eve[PPP1$MS_ID=="ProCan5",3],  eve[PPP1$MS_ID=="ProCan6", 3], xlab="MS 5", ylab="MS 6", bty='n')
abline(a=0,b=1, lwd=3, lty=2, col="grey")
dev.off()

pdf("Miss_vs_E1_2.pdf")
PC=table(PPP1$MS_ID)
layout(matrix(1:12, 3, 4, byrow = F))
par(mai=c(0.3, 0.5, 0.3, 0.1))
for(i in 1:2)
{
  for(j in 1:length(PC))
  {
    idx=which(PPP1$MS_ID == names(PC)[j])
    plot(col=ifelse(MisP[idx] < 500, "red", ifelse(MisP[idx] < 600 & MisP[idx] > 500, "cyan", "grey")), pch=ifelse(PPP1$Tissue[idx] == "T", 15, ifelse(PPP1$Tissue[idx] == "N", 16, 1)), main=paste0("PC", i, " ProCan ", j), eve[idx, 2+i], MisP[idx], xlim = 1.1*range(eve[,2+i]), ylim=1.1*range(MisP), xlab=paste("PC", i), ylab=paste("Missing count"), bty="n")
    mod=lm(MisP[idx]~eve[idx, 2+i])
    points(col="gold", mean(eve[idx, 2+i]), mean(MisP[idx]), pch=1, lwd=3, cex=3)
    abline(mod, lty=2, lwd=2, col="blue")
    legend(ifelse(i == 1, "topright", "topleft"), bty='n', legend = paste("Rsq=", format(summary(mod)$r.squared, digits = 3)), pch=1, col="white")
  }
}
dev.off()

pdf("Violin_MS.pdf")
vioplot(names=names(PC), MisP[PPP1$MS_ID == names(PC)[1]], MisP[PPP1$MS_ID == names(PC)[2]], MisP[PPP1$MS_ID == names(PC)[3]], MisP[PPP1$MS_ID == names(PC)[4]], MisP[PPP1$MS_ID == names(PC)[5]], MisP[PPP1$MS_ID == names(PC)[6]])
dev.off()

pdf("Miss_vs_AQUA.pdf")
layout(matrix(1:2, 2, 1))
plot(main="AQUA", eve[Aidx,3], MisP[Aidx], xlab="Eigen 1", ylab="Missing count", pch=16, bty='n')
md=lm(MisP[Aidx]~eve[Aidx,3])
abline(md, lwd=3, col="grey")
legend("topright", legend = paste("Rsq = ", format(summary(md)$r.squared, digits=3)), bty='n')

plot(main="CTRL", eve[Cidx,3], MisP[Cidx], xlab="Eigen 1", ylab="Missing count", pch=16, bty='n')
md=lm(MisP[Cidx]~eve[Cidx,3])
abline(md, lwd=3, col="grey")
legend("topright", legend = paste("Rsq = ", format(summary(md)$r.squared, digits=3)), bty='n')
dev.off()


#########Batch
#jpeg("pc_vs_Batch.jpeg", width = 600, height = 600)
pdf("pc_vs_Batch.pdf")
tb=names(table(PPP1$Batch))
or=order(as.numeric(tb))
mat=matrix(1:32, 4, 8)
layout(mat)
par(mai=c(0.15,0.15,0.15,0.15))
for(i in 1:length(tb))
{
  plot(xlim=range(eve[,3]), ylim=range(eve[,4]), eve[PPP1$Batch == tb[or[i]], 3], bty='n', eve[PPP1$Batch == tb[or[i]], 4], col=as.numeric(as.factor(PPP1$Tissue)), pch=16)
  legend("bottomleft", legend=tb[or[i]], bty='n')
}
dev.off()

############Date

#jpeg("PC_vs_Date.jpeg", width = 800, height = 800)
pdf("PC_vs_Date.pdf")

tDat=names(table(PPP1$Date))
mat=matrix(1:36, 6, 6)
layout(mat)
par(mai=c(0.15,0.15,0.15,0.15))
for(i in 1:length(tDat))
{
  plot(xlim=range(eve[,3]), ylim=range(eve[,4]), eve[PPP1$Date == tDat[i], 3], bty='n', eve[PPP1$Date == tDat[i], 4], pch=16, col=i)
  legend("bottomleft", legend=tDat[i], bty='n')
}
dev.off()

##############linear model cut for PCA
idx=which(PPP1$Tissue == "N" | PPP1$Tissue == "T")

exDat1=exDat[idx,]
PPP1_TS=PPP1$Tissue[idx]
PPP1_TS[PPP1_TS=="N"] = 0
PPP1_TS[PPP1_TS=="T"] = 1
pEx=matrix(2, ncol(exDat1))
for(i in 1:ncol(exDat1))
{
  Tc=length(which(PPP1_TS == 1 & !is.na(exDat1[,i])))
  Nc=length(which(PPP1_TS == 0 & !is.na(exDat1[,i])))
  if (Tc > 10 && Nc > 10)
  {
    mod=lm(PPP1_TS~exDat1[,i])
    pEx[i,1]=summary(mod)$coefficients[2,4]
  }
}

ACvar=matrix(0, 4, ncol(exDat))
for(i in 1:ncol(exDat))
{
  ACvar[1, i] = log10(sd(exDat[Aidx, i], na.rm = T))
  ACvar[2, i] = log10(sd(exDat[Cidx, i], na.rm = T))
  ACvar[3, i] = log10(sd(exDat[Nidx, i], na.rm = T))
  ACvar[4, i] = log10(sd(exDat[Tidx, i], na.rm = T))
}
rownames(ACvar)=c("SD(AQUA)", "SD(CTRL)", "SD(Normal)", "SD(Tumor)")
ACmean=matrix(0, 4, ncol(exDat))
for(i in 1:ncol(exDat))
{
  ACmean[1, i] = mean(log10(exDat[Aidx, i]), na.rm = T)
  ACmean[2, i] = mean(log10(exDat[Cidx, i]), na.rm = T)
  ACmean[3, i] = mean(log10(exDat[Nidx, i]), na.rm = T)
  ACmean[4, i] = mean(log10(exDat[Tidx, i]), na.rm = T)
}
rownames(ACmean)=c("Ave (AQUA)", "Ave (CTRL)", "Ave (Normal)", "Ave (Tumor)")

##AQUA, CTRL, Normal, Disease
pdf("AVE_SD.pdf")
TIT=c("AQUA", "CTRL", "Normal", "Tumor")
SDLab=c("SD (AQUA)", "SD (CTRL)", "SD (Normal)", "SD (Tumor)")
MLab=c("Ave (AQUA)", "Ave (CTRL)", "Ave (Normal)", "Ave (Tumor)")

mat=matrix(11, 4, 4)
mat[row(mat) >= col(mat)] = 1:10
mat[row(mat) < col(mat)] = 11:16
layout(mat)
par(mai=c(0.5,0.5,0.5,0.5))
for(i in 1:4)
{
  for(j in i:4)
  {
    if(i == j)
    {
      plot(x=ACvar[i,], y=ACmean[j,], xlab=SDLab[i], ylab=MLab[j], col="pink", xlim=c(0, 8), ylim=c(0, 8))
    } else {
      plot( xlab=SDLab[i], ylab=SDLab[j], 
            ACvar[i,], ACvar[j,], pch=16, cex=0.5, frame.plot = F, xlim=c(0, 8), ylim=c(0, 8))
      mod=lm(ACvar[i,]~ACvar[j,])
      abline(mod, lty=2, lwd=2, col="red")
      abline(a=0, b=1, lty=2, lwd=2, col="grey")
    }
  }
}

for(i in 1:3)
{
  for(j in (i+1):4)
  {
    plot(ylab= MLab[i], xlab= MLab[j],
         ACmean[j,], ACmean[i,], pch=16, cex=0.5, frame.plot = F, xlim=c(0, 8), ylim=c(0, 8), col="blue")
    mod=lm(ACmean[i,]~ACmean[j,])
    abline(mod, lty=2, lwd=2, col="red")
    abline(a=0, b=1, lty=2, lwd=2, col="grey")
  }
}

dev.off()

#####################F-value
#https://en.wikipedia.org/wiki/F-test_of_equality_of_variances
pdf("volcano.pdf")
F_ACvar=matrix(0, 6, ncol(exDat))
F_df1=matrix(0, 6, ncol(exDat))
F_df1[1,] = length(Cidx)
F_df1[2,] = F_df1[4,] = length(Nidx)
F_df1[3,] = F_df1[5,] = F_df1[6,] = length(Tidx)
F_df2=matrix(0, 6, ncol(exDat))
F_df2[1,] =F_df2[2,] = F_df2[3,] = length(Cidx)
F_df2[4,] = F_df2[5,] = length(Nidx)
F_df2[6,] = length(Nidx)

p_ACvar=matrix(0, 6, ncol(exDat))
for(i in 1:ncol(exDat))
{
  F_ACvar[1, i] = var(exDat[Cidx, i], na.rm = T)/var(exDat[Aidx, i], na.rm = T)
  F_ACvar[2, i] = var(exDat[Nidx, i], na.rm = T)/var(exDat[Aidx, i], na.rm = T)
  F_ACvar[3, i] = var(exDat[Tidx, i], na.rm = T)/var(exDat[Aidx, i], na.rm = T)
  F_ACvar[4, i] = var(exDat[Nidx, i], na.rm = T)/var(exDat[Cidx, i], na.rm = T)
  F_ACvar[5, i] = var(exDat[Tidx, i], na.rm = T)/var(exDat[Cidx, i], na.rm = T)
  F_ACvar[6, i] = var(exDat[Tidx, i], na.rm = T)/var(exDat[Nidx, i], na.rm = T)
  if(length(which(is.na(exDat[Aidx,i]))) > 0)
  {
    F_df2[1,i] = F_df2[1,i] - length(which(is.na(exDat[Aidx,i]))) - 1
    F_df2[2,i] = F_df2[2,i] - length(which(is.na(exDat[Aidx,i]))) - 1
    F_df2[3,i] = F_df2[3,i] - length(which(is.na(exDat[Aidx,i]))) - 1
  }
  
  if(length(which(is.na(exDat[Cidx,i]))) > 0)
  {
    F_df1[1,i] = F_df1[1,i] - length(which(is.na(exDat[Cidx,i]))) - 1
    F_df2[4,i] = F_df2[4,i] - length(which(is.na(exDat[Cidx,i]))) - 1 
    F_df2[5,i] = F_df2[5,i] - length(which(is.na(exDat[Cidx,i]))) - 1
  }
  
  if(length(which(is.na(exDat[Nidx,i]))) > 0)
  {
    F_df1[2,i] = F_df1[2,i] - length(which(is.na(exDat[Nidx,i]))) - 1
    F_df1[4,i] = F_df1[4,i] - length(which(is.na(exDat[Nidx,i]))) - 1
    F_df2[6,i] = F_df2[6,i] - length(which(is.na(exDat[Nidx,i]))) - 1
  }
  
  if(length(which(is.na(exDat[Tidx,i]))) > 0)
  {
    F_df1[3,i] = F_df1[3,i] - length(which(is.na(exDat[Tidx,i]))) - 1
    F_df1[5,i] = F_df1[5,i] - length(which(is.na(exDat[Tidx,i]))) - 1
    F_df1[6,i] = F_df1[6,i] - length(which(is.na(exDat[Tidx,i]))) - 1
  }
  
}

for(i in 1:6)
{
  p_ACvar[i,]=pf(F_ACvar[i,], F_df1[i,], F_df2[i, ], lower.tail = F)
}

mat=matrix(7, 3, 3)
mat[row(mat) >= col(mat)] = 1:6
mat[1,2]=9
mat[2,3]=8
layout(mat)
par(mai=c(0.5,0.5,0.5,0.5))
MN=c("V(C)/V(A)", "V(N)/V(A)", "V(T)/V(A)", "V(N)/V(C)", "V(T)/V(C)", "V(T)/V(N)")

for(i in 1:6)
{
  hist(p_ACvar[i,], main=MN[i], xlab="F-test, P-value", breaks = 25, col=ifelse(i == 6, "red", "grey"))
  abline(h=500, col="red", lty=2)
}

dev.off()

##############linear model
idx=which(PPP1$Tissue == "N" | PPP1$Tissue == "T")

exDat1=exDat[idx,]
PPP1_TS=PPP1$Tissue[idx]
PPP1_TS[PPP1_TS=="N"] = 0
PPP1_TS[PPP1_TS=="T"] = 1
pEx=matrix(NA, ncol(exDat1))
fchange=matrix(NA, ncol(exDat1))
for(i in 1:ncol(exDat1))
{
  Tm=which(PPP1_TS == 1 & !is.na(exDat1[,i]))
  Nm=which(PPP1_TS == 0 & !is.na(exDat1[,i]))
  Tc=length(Tm)
  Nc=length(Nm)
  if(Tc > 10 && Nc > 10)
  {
    mod=lm(PPP1_TS~exDat1[,i])
    pEx[i,1]=summary(mod)$coefficients[2,4]
    fchange[i,1]=mean(exDat1[Tm,i])/mean(exDat1[Nm,i])
  }
}

idxNA= unique(c(which(is.na(pEx)), which(fchange[,1] < 0)))
ProNames=colnames(exDat1)[-idxNA]
m1=fchange[-idxNA,]
pEx=pEx[-idxNA,]


pdf("Volcano_1.pdf")
layout(matrix(1:2, 2, 1))
plot(main="Bonferroni 0.01", x=log2(m1), xlim=c(-max(abs(log2(m1))), max(abs(log2(m1)))), -log10(pEx), cex=0.5, pch=16, bty='n', xlab="Fold change log2 scale", ylab=expression(paste(log[10](p))), col="grey")
abline(v=c(-1,1), lty=2)
idxG1=intersect(which(log2(m1) > 1), which(-log10(pEx) > -log10(0.01/length(pEx))))
points(log2(m1[idxG1]), -log10(pEx[idxG1]), col="red", cex=0.5, pch=16)
idxG2=intersect(which(log2(m1) < -1), which(-log10(pEx) > -log10(0.01/length(pEx))))
points(log2(m1[idxG2]), -log10(pEx[idxG2]), col="blue", cex=0.5, pch=16)

legend("topleft", legend = length(idxG2) + length(idxG1), bty='n')
legend("topright", legend = c("Tumor > Normal", "Tumor < Normal"), pch=16, col=c("red", "blue"), bty='n')

padj=p.adjust(pEx)

plot(main = "FDR 0.01", x=log2(m1), xlim=c(-max(abs(log2(m1))), max(abs(log2(m1)))), -log10(pEx), cex=0.5, pch=16, bty='n', xlab="Fold change log2 scale", ylab=expression(paste(log[10](p))), col="grey")
abline(v=c(-1,1), lty=2)
idxG1=intersect(which(log2(m1) > 1), which(-log10(padj) > -log10(0.01/length(pEx))))
points(log2(m1[idxG1]), -log10(pEx[idxG1]), col="red", cex=0.5, pch=16)
idxG2=intersect(which(log2(m1) < -1), which(padj < 0.01))
points(log2(m1[idxG2]), -log10(pEx[idxG2]), col="blue", cex=0.5, pch=16)
legend("topleft", legend = length(idxG2) + length(idxG1), bty='n')
dev.off()

idxG1=intersect(which(log2(m1) > 1), which(-log10(padj) > -log10(0.01/length(pEx))))
idxG2=intersect(which(log2(m1) < -1), which(-log10(padj) > -log10(0.01/length(pEx))))
write.table(ProNames[idxG1], "TumorHighExpProtein.txt", row.names = F, col.names = F, quote = F)
write.table(ProNames[idxG2], "TumorLowExpProtein.txt", row.names = F, col.names = F, quote = F)


######
eve=read.table("prot.eigenvec", as.is = T)
eve$V1=toupper(eve$V1)
eve$V2=toupper(eve$V2)
eval=read.table("prot.eigenval", as.is = T)
FT=10
Wt=array(eval[1:FT, 1]/sum(eval[1:FT, 1]), FT)

PPP0=read.csv("PPPB1_lineup.csv", as.is = T, header = T)
idxM=c()
for(i in 1:nrow(eve))
{
  idxM=c(idxM, which(PPP0$PPPA_ID == eve$V1[i]))
}
if(length(idxM) > 0)
{
  PPP1=PPP0[idxM,]
} else {
  PPP1=PPP0
}

pdf("PVCA.pdf")
layout(matrix(1:6, 2, 3))
TN=c("Tissue", "MS", "Batch")
wRsq=array(0, 3)
for(ii in 1:3)
{
  if (ii == 1)
  {
    idxT=which(PPP1$Tissue == "N" | PPP1$Tissue == "T")
    y=PPP1$Tissue[idxT]
  } else if (ii == 2) {
    idxT=which(!is.na(PPP1$MS_ID))
    y=PPP1$MS_ID[idxT]
  } else {
    idxT=which(!is.na(PPP1$Batch))
    y=PPP1$Batch[idxT]
  }
  
  Rsq=array(0, FT)
  names(Rsq)=seq(1,10)
  for (i in 1:FT)
  {
    an=aov(eve[idxT,2+i]~y)
    Rsq[i]=summary(an)[[1]][1,2]/sum(summary(an)[[1]][,2])
  }

  barplot(main=TN[ii], Rsq, ylim=c(0, 0.3), ylab=expression(R^2), border = F)
  barplot(Rsq*Wt, ylim=c(0, 0.15), ylab=expression(paste("Weighted ", R^2)), xlab="PC", border = F)
  wRsq[ii] = sum(Rsq*Wt)
}
dev.off()

