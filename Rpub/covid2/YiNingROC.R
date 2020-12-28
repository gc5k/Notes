dat=read.csv("Book1.csv")
library(pROC)

layout(matrix(1:6, 2, 3))
for(i in 1:3) {
  sc=dat[,i]
  plot(main=paste0("Model ", i), density(sc[dat$Diag==0]), col="green", xlim=range(sc), xlab="Predicted score")
  lines(density(sc[dat$Diag==1]), col="blue")
  rug(sc[dat$Diag==0], col="green", line=0)
  rug(sc[dat$Diag==1], col="blue", line=-0.5)
  legend("topright", legend = c("NCP", "Non-NCP"), col=c("blue","green"), lty=1, bty = 'n')
  a1=roc(dat$Diag, sc, ci=TRUE)
  print(a1)
  ci1=ci.auc(a1)
  ci_se=ci.se(a1, boot.n = 100, specificities = seq(0, 1, 0.01))
  plot.roc(main=paste0("AUC (95% CI)~", format(a1$auc, digits = 3), " (", format(ci1[1], digits=3), ", ", format(ci1[3], digits = 3),")"),
           dat$Diag, sc, print.thres="best", print.thres.best.method="youden", col="red")
  #plot(smooth(a1))#smooth
  #lines(1-ci_se[,1], ci_se[,2])
#  plot(ci_se)
  lines(as.numeric(rownames(ci_se)), ci_se[,1], col="grey", lty=2)
  lines(as.numeric(rownames(ci_se)), ci_se[,3], col="grey", lty=2)
}

MD=c("WBC", "X-Ray")
layout(matrix(1:6, 2, 3))
for(i in 1:2) {
  sc=dat[,i+3]
  plot(main=MD[i], density(sc[dat$Diag==0]), col="green", xlim=range(sc), xlab="Predicted score")
  lines(density(sc[dat$Diag==1]), col="blue")
  rug(sc[dat$Diag==0], col="green", line=0)
  rug(sc[dat$Diag==1], col="blue", line=-0.5)
  legend("topright", legend = c("NCP", "Non-NCP"), col=c("blue","green"), lty=1, bty = 'n')
  a1=roc(dat$Diag, sc, ci=TRUE)
  print(a1)
  ci1=ci.auc(a1)
  ci_se=ci.se(a1, boot.n = 100, specificities = seq(0, 1, 0.01))
  plot.roc(xlim=c(1,0), ylim=c(0,1),main=paste0("AUC (95% CI)~", format(a1$auc, digits = 3), " (", format(ci1[1], digits=3), ", ", format(ci1[3], digits = 3),")"),
           dat$Diag, sc, print.thres="best", print.thres.best.method="youden", col="red")
  #plot(smooth(a1))#smooth
  #lines(1-ci_se[,1], ci_se[,2])
  #  plot(ci_se)
  lines(as.numeric(rownames(ci_se)), ci_se[,1], col="grey", lty=2)
  lines(as.numeric(rownames(ci_se)), ci_se[,3], col="grey", lty=2)
}
