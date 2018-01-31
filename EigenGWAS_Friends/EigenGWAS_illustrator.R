
FN="Arab295"
PC=10

source("EigenGWAS_Friends.R")

####GRM stats
grm=grmReader(FN)
grmS=grm[col(grm) > row(grm)]
ne=-1/mean(grmS)
me=1/var(grmS)
hist(grmS, xlab="GRM scores", main ="GRM distribution", breaks=25)
legend("topright", legend = c(paste0("ne=", ne), paste0("me=", me)), bty='n')

############PC plot
Evec=read.table(paste0(FN, ".eigenvec"), as.is = T)
pdf(paste0(FN, ".pdf"))
plot(Evec[,3], Evec[,4], xlab="PC 1", ylab="PC 2", frame.plot = F)
dev.off()

####EigenValue vs Lambda_GC
Evev=read.table(paste0(FN, ".eigenval"), as.is = T)
GC=array(0, dim=PC)

for(i in 1:PC)
{
  eg = read.table(paste0(FN, ".", i, ".egwas"), as.is = T, header = T)

  GC[i] = qchisq(median(eg$P), 1, lower.tail = F)/qchisq(0.5, 1)

  pdf(paste0(FN, "_", i, ".pdf"))
  layout(matrix(1:3, 3, 1))
  eg=eg[,-which(colnames(eg)=="P")]
  colnames(eg)[which(colnames(eg)=="PGC")]="P"
  manhattan(eg, pch=16, cex=0.5)
  FstPlot(eg, pch=16, cex=0.5)
  plot(eg$Chi, eg$Fst, xlab=expression(chi^2), ylab="Fst", pch=16, cex=0.5)
  dev.off()
}

egc=matrix(c(Evev[1:PC,1], GC), PC, 2, byrow = F)
rownames(egc)=seq(1, PC)
barplot(t(egc), beside = T, border = F)
legend("topright", legend = c("Eigenvalue", expression(paste(lambda[gc]))), pch=15, col=c("black", "grey"), bty='n')

########manhattan plot

