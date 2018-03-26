manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=NULL, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (!is.null(limitchromosomes)) {
    d=d[d$CHR %in% limitchromosomes, ]
  }
  
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0 #  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  colors <- rep(colors,max(d$CHR))[1:length(unique(d$CHR))]
  
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    Uchr=unique(d$CHR)
    for (i in 1:length(Uchr)) {
      if (i==1) {
        d[d$CHR==Uchr[i], ]$pos=d[d$CHR==Uchr[i], ]$BP
      } else {
        lastbase=lastbase+tail(subset(d, CHR==Uchr[i-1])$BP, 1)
        d[d$CHR==Uchr[i], ]$pos=d[d$CHR==Uchr[i], ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==Uchr[i], ]$pos[floor(length(d[d$CHR==Uchr[i], ]$pos)/2)+1])
    }
  }
  if (numchroms==1) {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  } else {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    Uchr=unique(d$CHR)
    for (i in 1:length(Uchr)) {
      with(d[d$CHR==Uchr[i], ], points(pos, logp, col=colors[icol], ...))
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
