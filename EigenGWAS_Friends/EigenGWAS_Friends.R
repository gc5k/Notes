grmStats <- function(root, pdf="n")
{
  grm=grmReader(root)
  grmS=grm[col(grm) > row(grm)]
  ne=-1/mean(grmS)
  me=1/var(grmS)
  if (pdf == "pdf")
  {
    pdf(paste0(root, "_grm.pdf"))
  }
  hist(grmS, xlab="GRM scores", main ="GRM distribution", breaks=25)
  legend("topright", legend = c(paste0("ne=", format(ne, digits=3, nsmall=2)), paste0("me=", format(me, digits=3, nsmall=2))), bty='n')

  if (pdf == "pdf")
  {
    dev.off()
  }
}

pcPlot <- function(root, pdf="n")
{
  Evec=read.table(paste0(root, ".eigenvec"), as.is = T)
  if(pdf == "pdf")
  {
    pdf(paste0(root, "_pc.pdf"))
  }
  plot(Evec[,3], Evec[,4], pch=16, xlab="PC 1", ylab="PC 2", frame.plot = F)
  if(pdf == "pdf")
  {
    dev.off()
  }
}

EigenValuePlot <- function(root, PC, pdf="n")
{
  Evev=read.table(paste0(root, ".eigenval"), as.is = T)
  GC=array(0, dim=PC)

  for(i in 1:PC)
  {
    eg = read.table(paste0(root, ".", i, ".egwas"), as.is = T, header = T)
    GC[i] = qchisq(median(eg$P), 1, lower.tail = F)/qchisq(0.5, 1)
  }

  egc=matrix(c(Evev[1:PC,1], GC), PC, 2, byrow = F)
  rownames(egc)=seq(1, PC)
  if(pdf == "pdf")
  {
    pdf(paste0(FN, "_EV.pdf"))
  }
  barplot(t(egc), beside = T, border = F)
  legend("topright", legend = c("Eigenvalue", expression(paste(lambda[gc]))), pch=15, col=c("black", "grey"), bty='n')
  if(pdf == "pdf")
  {
    dev.off()
  }
}

DeepEigenValuePlot <- function(root, PC, cutV, pdf="n")
{
  Evev=read.table(paste0(root, ".eigenval"), as.is = T)
  GC=matrix(0, nrow=1, ncol=length(cutV)+1)
  GC[1]=Evev[PC,1]
  eg = read.table(paste0(root, ".", PC, ".egwas"), as.is = T, header = T)

  for(i in 1:length(cutV))
  {
    GC[1,1+i] = qchisq(sort(eg$P, T)[ceiling(cutV[i]*nrow(eg))], 1, lower.tail = F)/qchisq(cutV[i], 1, lower.tail = T)
  }
  colnames(GC)=c("EV", cutV)
  par(las=2)
  barplot(GC, ylim=c(0, max(GC)*1.3), beside = T, border = F, col=c("black", rep("grey", length(cutV))))
  legend("topright", legend = c("Eigenvalue", expression(paste(lambda[gc]))), pch=15, col=c("black", "grey"), bty='n')
  if(pdf == "pdf")
  {
    dev.off()
  }
}


EigenGWASPlot <- function(root, pc, pdf='n')
{
  eg=read.table(paste0(root, ".", pc, ".egwas"), as.is = T, header = T)

  if(pdf == "pdf")
  {
    pdf(paste0(root, "_E", pc, ".pdf"))
  }
  layout(matrix(1:3, 3, 1))
  
  eg=eg[,-which(colnames(eg)=="P")]
  colnames(eg)[which(colnames(eg)=="PGC")]="P"
  manhattan(eg, pch=16, cex=0.5)
  FstPlot(eg, pch=16, cex=0.5)
  plot(eg$Chi, eg$Fst, xlab=expression(chi[1]^2), ylab=expression(paste(italic("F")[italic("ST")])), pch=16, cex=0.5, frame.plot = F)
  rsq=cor(eg$Chi, eg$Fst, use="pairwise.complete.obs")^2
  legend("topright", legend = c(paste("Rsq=", format(rsq, digits = 4, nsmall = 3))), bty='n')

  if(pdf == "pdf")
  {
    dev.off()
  }
}

SWEigenGWASPlot <- function(root, pc, kb=5, pdf='n')
{
  eg=read.table(paste0(root, ".", pc, ".egwas"), as.is = T, header = T)
  lgc=qchisq(median(eg$P), 1, lower.tail = F)/qchisq(0.5, 1)
  CHR=names(table(eg$CHR))
  for(i in 1:length(CHR))
  {
    egs=eg[which(eg$CHR == CHR[i]),]
    rg=range(egs$BP)
    meg=matrix(0, nrow(egs), 5)
    
    for(j in 1:nrow(meg))
    {
      idx=which((egs$BP >= egs$BP[j] - kb*1000) & (egs$BP <= egs$BP[j] + kb * 1000))
      meg[j, 1] = egs$CHR[j]
      meg[j, 2] = egs$BP[j]
      meg[j, 3] = mean(egs$Fst[idx])
      meg[j, 4] = pchisq(mean(egs$Chi[idx]/lgc), 1, lower.tail = F)
      meg[j, 5] = mean(egs$Chi[idx])
    }

    if(i == 1)
    {
      Meg = meg
    }
    else {
      Meg = rbind(Meg, meg)
    }
  }
  DMeg=as.data.frame(Meg)
  colnames(DMeg)=c("CHR", "BP", "Fst", "P", "Chi")

  
  if(pdf == "pdf")
  {
    pdf(paste0(root, "_SE", pc, ".pdf"))
  }
  layout(matrix(1:3, 3, 1))

  manhattan(DMeg, cex=0.5, pch=16)
  FstPlot(DMeg, pch=16, cex=0.5)
  plot(DMeg$Chi, DMeg$Fst, xlab=expression(chi[1]^2), ylab=expression(paste(italic("F")[italic("ST")])), pch=16, cex=0.5, frame.plot = F)
  rsq=cor(DMeg$Chi, DMeg$Fst, use="pairwise.complete.obs")^2
  legend("topright", legend = c(paste("Rsq=", format(rsq, digits = 4, nsmall = 3))), bty='n')
  
  if(pdf == "pdf")
  {
    dev.off()
  }
}


grmReader <- function(root)
{
  grmFile = gzfile(paste0(root, ".grm.gz"))
  grm = read.table(grmFile)
  
  idFile = read.table(paste0(root, ".grm.id"), as.is=T)
  
  mat=matrix(0, nrow(idFile), nrow(idFile))
  for(i in 1:nrow(idFile))
  {
    idx1 = i * (i-1)/2 + 1
    idx2 = i * (i+1)/2
    mat[i, 1:i] = grm[idx1:idx2, 4]
    mat[1:i, i] = grm[idx1:idx2, 4]
  }
  return(mat)
}

manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  if (numchroms==1) {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), frame.plot = F, ...))
  }  else {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", frame.plot = F, ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...))
  }
  #  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (!is.null(genomewideline)) {
    abline(h=genomewideline, col="gray")
  } else {
    abline(h=-log10(0.05/nrow(d)), col="gray")    
  }
}

FstPlot <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=NULL, genomewideline=NULL, annotate=NULL, title="", ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "Fst" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (Fst>0 & Fst<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = d$Fst
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-max(d$logp)*1.1
  #  if (ymax<8) ymax<-8
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  if (numchroms==1) {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(italic("F")[italic("ST")])), xlab=paste("Chromosome",unique(d$CHR),"position"), main=title, frame.plot=F, ...))
  }  else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(italic("F")[italic("ST")])), xlab="Chromosome", xaxt="n", type="n", main=title, frame.plot = F, ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...))
  }
  #  if (suggestiveline) abline(h=suggestiveline, col="blue")
  #  if (genomewideline) abline(h=genomewideline, col="red")
}
