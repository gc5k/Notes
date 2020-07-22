dat=read.csv("20200721all_cgb.csv", as.is = T, header = T)


##basic accessment of the data
###fivenum, mean, sd, missing, isNumeric
cutoff=5 #outlier pickup

SUMmat=matrix(0, ncol(dat), 9)
colnames(SUMmat)=c("minimum", "lower-hinge", "median", "upper-hinge", "maximum", "mean", "sd", "missingRate", "isNumeric")
rownames(SUMmat)=colnames(dat)
for(i in 1:nrow(SUMmat)) {
  if(is.numeric(dat[,i])) {
    SUMmat[i,1:5]=fivenum(dat[,i])
    SUMmat[i,6]=mean(dat[,i], na.rm = T)
    SUMmat[i,7]=sd(dat[,i], na.rm = T)
    SUMmat[i,8]=length(which(is.na(dat[,i])))/nrow(dat)
    SUMmat[i,9]=T
    
    idx=which(dat[,i] > SUMmat[i,6]+cutoff*SUMmat[i,7])
    if(length(idx)>0) {
      for(j in 1:length(idx)) {
        print(paste0("sample id: ",dat[idx[j],2], ", col:", colnames(dat)[i], ", ", "outlier value: ",  dat[idx[j], i]))
      }
      
    }
  } else {
    SUMmat[i,1:8]=NA
    SUMmat[i,9]=F
  }
}

par(las=2, font=5, bty='l', ps=6)
boxplot(dat[,-2], cex=0.5, pch=16, col=runif(nrow(dat), 1, nrow(dat)))
barplot(SUMmat[,7]/SUMmat[,6])
