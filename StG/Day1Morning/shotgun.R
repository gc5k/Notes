
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
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }  else {
    with(d, plot(main=title, pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
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

Dprime2Correlation <- function(freq, dprime)
{
  ld=Dprime2LD(freq, dprime)
  corl=array(0, dim=length(ld))
  for(i in 1:length(corl))
  {
    corl[i] = ld[i]/sqrt(freq[i]*(1-freq[i])*freq[i+1]*(1-freq[i+1]))
  }
  return(corl)
}

LD2Correlation <- function(freq, ld)
{
  corl=array(0, dim=length(ld))
  for(i in 1:length(corl))
  {
    corl[i] = ld[i]/sqrt(freq[i]*(1-freq[i])*freq[i+1]*(1-freq[i+1]))
  }
  return(corl)
}

Dprime2LD <- function(freq, dprime)
{
  if(length(freq) == length(dprime))
  {
    ld=array(0, dim=length(dprime)-1)
  }
  else if(length(dprime) == (length(freq) -1) )
  {
    ld=array(0, dim=length(dprime))
  }
  
  for(i in 1:length(ld))
  {
    if(dprime[i] > 0)
    {
      ld[i] = dprime[i] * min(freq[i] * (1-freq[i+1]), (1-freq[i])*freq[i+1])
    }
    else
    {
      ld[i] = dprime[i] * min(freq[i] * (freq[i+1]), (1-freq[i]) * (1-freq[i+1]))
    }
  }
  return(ld)
}

GenerateGenoDprime <- function(freq, Dprime, N)
{
  g = matrix(0, nrow=N, ncol=length(freq))
  
  ld=Dprime2LD(freq, Dprime)
  for(h in 1:N)
  {
    gMat = matrix(0, nrow=length(freq), ncol=2)
    for(i in 1:length(freq))
    {
      for(j in 1:2)
      {
        idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1)
        if(i == 1)
        {
          gMat[i,j] = idx
        }
        else
        {
          d = runif(1, 0, 1)
          a = gMat[i-1,j]
          f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
          f2 = ifelse(a == 0, freq[i], 1-freq[i])
          gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j])
        }
      }
    }
    g[h,] = gMat[,1] + gMat[,2]
  }
  return(g)
}


GenerateGeno <- function(freq, ld, N)
{
  g = matrix(0, nrow=N, ncol=length(freq))

  for(h in 1:N)
  {
    gMat = matrix(0, nrow=length(freq), ncol=2)
    for(i in 1:length(freq))
    {
      for(j in 1:2)
      {
        idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1)
        if(i == 1)
        {
          gMat[i,j] = idx
        }
        else
        {
          d = runif(1, 0, 1)
          a = gMat[i-1,j]
          f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
          f2 = ifelse(a == 0, freq[i], 1-freq[i])
          gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j])
        }
      }
    }
    g[h,] = gMat[,1] + gMat[,2]
  }
  return(g)
}

GenerateHaplo <- function(freq, ld, N)
{
  g = matrix(0, nrow=N*2, ncol=length(freq))
  
  for(h in 1:N)
  {
    gMat = matrix(0, nrow=length(freq), ncol=2)
    for(i in 1:length(freq))
    {
      for(j in 1:2)
      {
        idx = ifelse(runif(1, 0, 1) < freq[i], 0, 1)
        if(i == 1)
        {
          gMat[i,j] = idx
        }
        else
        {
          d = runif(1, 0, 1)
          a = gMat[i-1,j]
          f1 = ifelse(a == 0, freq[i-1], 1-freq[i-1])
          f2 = ifelse(a == 0, freq[i], 1-freq[i])
          gMat[i,j] = ifelse(d < (f1 * f2 +ld[i-1])/f1, gMat[i-1,j], 1-gMat[i-1,j])
        }
      }
    }
    g[h*2-1,] = gMat[,1]
    g[h*2,] = gMat[,2]
  }
  return(g)
}


PlinkBedReader <- function(root, impute=F, loci=0)
{
  
  fam=read.table(paste0(root, ".fam"))
  fam$UID=paste0(fam[,1], fam[,2], ".")
  if(anyDuplicated(fam$UID) > 0) {
    stop('duplciated individual ids.')
  }
  bim=read.table(paste0(root, ".bim"))
  
  bin.connection = file(paste0(root, ".bed"), 'rb')
  test.bytes = readBin(bin.connection, what = "raw", n = 3)
  if(!identical(as.character(test.bytes), c('6c', '1b', '01'))) {
    stop('BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file')
  }
  
  if (loci <= 0)  {
    loci = nrow(bim)
  }
  len=ceiling(nrow(fam)/4)
  genotype = array(NA, dim = c(nrow(fam), loci))
  for (k in 1:loci) {
    r.bin.snp = readBin(bin.connection, what = 'raw', n = len)
    bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2, byrow = TRUE)[1:nrow(fam),]
    genotype[,k] = bin.snp[,1] + bin.snp[,2] - 10 * ((bin.snp[,1] == 1) & (bin.snp[,2] == 0))
    idx=which(genotype[,k] == -9)
    if(length(idx) > 0)
    {
      genotype[idx,k] = NA
      if(impute)
      {
        frq=mean(genotype[,k], na.rm=T)/2
        if(frq > 0)
        {
          genotype[idx,k] = rbinom(length(idx), 2, frq)          
        } else {
          if(length(idx) < nrow(genotype))
          {
            genotype[idx,k] = 0
          } else {
            genotype[idx,k] = 0
          }
        }
      }
    }
  }
  close(bin.connection)
  return(genotype)
}

lambda <- function(m, k, rep=10000)
{
  SimuLam=matrix(0, rep, 1)
  P=rbeta(rep, k, m-k+1)
  SimuLam=qchisq(P, 1, lower.tail=F)/qchisq(k/(m+1), 1, lower.tail=F)
  lambda=c(mean(SimuLam), sd(SimuLam))
  return(lambda)
}
