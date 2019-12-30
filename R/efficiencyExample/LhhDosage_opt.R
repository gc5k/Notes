conn <- gzfile("lr_10000.plk.ped.gz", "rt") #repalce the filename with your own ped
ped=read.table(conn, as.is = T)[,-c(1:6)]#remove the first 6 cols
hped=ped[,seq(1, ncol(ped), 2)] #inbred, one haploid is enough

##QC stats
freq=colMeans(hped, na.rm = T) #freq
sMiss=array(0, dim=nrow(hped)) #sample-level missing
lMiss=array(0, dim=ncol(hped)) #locus-leve missing
hped_T=t(hped)
for(i in 1:length(sMiss)) {
  sMiss[i]=length(which(is.na(hped_T[,i])))
}
rm(hped_T)

for(i in 1:length(lMiss)) {
  lMiss[i]=length(which(is.na(hped[,i])))
}

layout(matrix(1:3, ncol=3))
plot(main="Frequency", freq, pch=16, cex=0.5)
plot(main="Individual missing rate", sMiss/ncol(hped), pch=16, cex=0.5)
plot(main="Locus missing rate", lMiss/nrow(hped), pch=16, cex=0.5)

#########clean
#QC 1 for MAF
QCcut=list("fqMin"=0.03, "fqMax"=0.97, "IndMiss"=0.4)
fQC=which(freq > QCcut$fqMin & freq < QCcut$fqMax) #remove freq <0.03 and > 0.97
#QC 2 for ind miss
iQC=which(sMiss/ncol(hped) < QCcut$IndMiss) #individual genotyping missing rate

hpedQC=hped[iQC,fQC]
hpedQC_T=t(hpedQC)
freqQC=colMeans(hpedQC, na.rm = T)
vQC=freqQC*(1-freqQC)

##test sample size, the maxima of testInd=nrow(hpedQC)
#testInd=100
testInd=nrow(hpedQC)

#make grm
G=matrix(0, nrow = testInd, ncol = testInd)
for(i in 1:testInd) {
  print(paste0(i, "/", testInd))

##slow  x1=ifelse(!is.na(as.numeric(hpedQC[i,])), as.numeric(hpedQC[i,]), freqQC) 
  ##fast
  x1=ifelse(!is.na(as.numeric(hpedQC_T[,i])), as.numeric(hpedQC_T[,i]), freqQC)
  s1=(x1-freqQC)/sqrt(vQC)
  for(j in i:testInd) {
##slow    x2=ifelse(!is.na(as.numeric(hpedQC[j,])), as.numeric(hpedQC[j,]), freqQC)
    ##fast
    x2=ifelse(!is.na(as.numeric(hpedQC_T[,j])), as.numeric(hpedQC_T[,j]), freqQC)
    
    s2=(x2-freqQC)/sqrt(vQC)
    effM=length(freqQC)-length(which(x1==0 | x2==0))
    G[j,i]=G[i,j]=sum(s1*s2)/effM
  }
}
rm(hpedQC_T)
write.table(G, "G.txt", row.names = F, col.names = F, quote = F)
g1=read.table("G.txt", as.is = T)
eG=eigen(g1)
write.table(eG$values, "Gvalue.txt", row.names = F, col.names = F, quote = F)
write.table(eG$vectors, "Gvec.txt", row.names = F, col.names = F, quote = F)

eVec=read.table("Gvec.txt", as.is = T)
eVal=read.table("Gvalue.txt", as.is = T)
layout(matrix(1:2, 1, 2))
barplot(eVal[1:5,1])
plot(eVec$V1, eVec$V2, pch=16, cex=0.5)

MOD=matrix(2, nrow=ncol(hpedQC), 4)

for(i in 1:nrow(MOD)) {
  mod=lm(eVec$V1~hpedQC[1:testInd,i])
  sm=summary(mod)
  if(nrow(sm$coefficients)>1) {
    MOD[i,]=sm$coefficients[2,]
  }
}
MOD1=MOD[MOD[,4]!=2,]
colnames(MOD1)=colnames(sm$coefficients)
hist(MOD1[,4])
gc=qchisq(median(MOD1[,4]),1,lower.tail = F)/0.455
