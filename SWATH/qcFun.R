#############function pca
PCA <- function(exDat, fileN, ID)
{
  exDatS=matrix(0, nrow(exDat), ncol(exDat))
  
  for(i in 1:ncol(exDat))
  {
    exDatS[,i] = scale(exDat[,i])
  }
  
  exCor=matrix(0, nrow(exDat), nrow(exDat))
  datM=matrix(0, nrow=nrow(exCor) * (nrow(exCor)+1)/2, 4)
  cnt=1
  for(i in 1:nrow(exDatS))
  {
    for(j in 1:i)
    {
      idxD=which(!is.na(exDatS[i,]) & !is.na(exDatS[j,]))
      exCor[i,j]=exCor[j,i]=mean(exDatS[i,idxD] * exDatS[j,idxD])
      datM[cnt,3] = length(idxD)
      cnt=cnt+1
    }
  }
  
  cnt=1
  for(i in 1:nrow(exCor))
  {
    for(j in 1:i)
    {
      datM[cnt,1] = i
      datM[cnt,2] = j
      datM[cnt,4] = exCor[i,j]
      cnt=cnt+1
    }
  }
  
  matS=matrix(0, nrow(exCor), 2)
  for(i in 1:nrow(exCor))
  {
    matS[i,1] = rownames(exDat)[i]
    matS[i,2] = rownames(exDat)[i]
  }
  
  write.table(matS, paste0(fileN, ".grm.id"), row.names = F, col.names = F, quote = F)
  write.table(datM, paste0(fileN, ".grm"), row.names = F, col.names = F, quote = F)
  system(paste0("gzip -f ", fileN, ".grm"))
  ei = eigen(exCor)
  write.table(ei$values, paste0(fileN, ".eigenval"), row.names = F, col.names = F, quote = F)
  evec=matrix(0, nrow=nrow(exCor), ncol = 12)
  evec[,1]=ID
  evec[,2]=ID
  evec[,3:12]=ei$vectors[, 1:10]
  write.table(evec, paste0(fileN, ".eigenvec"), row.names = F, col.names = F, quote=F)
}
