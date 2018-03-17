gcta='/Users/gc5k/bin/gcta_1.02/gcta_mac'
gear='java -jar /Users/gc5k/Documents/workspace/FromSVN/GEAR/gear.jar'
inode='ssh gc5k@inode.qbi.uq.edu.au'
plink='/Users/gc5k/bin/plink-1.07-mac-intel/plink'
hapgen='/Users/uqgchen5/bin/hapgen2_macosx_intel/hapgen2'
HM3='/Users/uqgchen5/bin/hapgen2_macosx_intel/HM3'

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#### input: bv, disease 0/1, sex
#### output: matrix(1:10, c(decile, or.lo, or, or.hi))
getDecileOrs = function(bv, phen, sex=NA) {
  mat = matrix(NA, 4, 10 )
  d = bv
  # find decile cutoffs
  q = quantile(d, seq(0, 1, len=11))
  dec = list()
  # group in deciles
  for(i in 1:(length(q)-1)) { dec[[i]] = which(d >= q[i] & d < q[i+1])}
  for(decile in 1:10) {
    d.top = dec[[decile]]
    d.bot = dec[[1]]
    # group high and low decile
    tb = c(d.top, d.bot)
    # control for sex, if available
    if(!all(is.na(sex)))
      df = data.frame(phen=phen[tb], group=c(rep(1, length(d.top)), rep(0, length(d.bot))), sex=sex[tb])
    else
      df = data.frame(phen=phen[tb], group=c(rep(1, length(d.top)), rep(0, length(d.bot))))
    logm = glm(phen ~ ., data=df, family=binomial(logit))
    #logm = glm(phen ~ group, data=df, family=binomial(logit))
    est = summary(logm)$coefficients[2,1]
    std.err = as.numeric(summary(logm)$coefficients[2,2])
    
    # get OR from logistic model beta
    or.adj = as.numeric(exp(est))
    # 95% confidence interval
    or.adj.lo = exp(est - 1.96 * std.err)
    or.adj.hi = exp(est + 1.96 * std.err)
    pred=1
    mat[1, decile] = decile
    mat[2, decile] = or.adj.lo
    mat[3, decile] = or.adj
    mat[4, decile] = or.adj.hi
  }
  out = data.frame(t(mat))
  names(out) = c('decile', 'OR.low', 'OR', 'OR.high')
  return(out)
}

plot.deciles.2 = function(df, name='profile score decile ORs', size=12) {
  
  require(ggplot2)
  require(grid)
  
  th2 = theme(panel.background = element_rect(fill = "white", colour='black'), text=element_text(size=size), legend.title=element_blank(), axis.title.y=element_text(vjust=2), axis.title.x=element_text(vjust=-2.4), plot.margin=unit(c(5,5,15,12),"mm"), panel.grid=element_blank(), legend.position="top")
  limits = aes(ymax=OR.high, ymin=OR.low)
  ar=.5
  ymax=10
  ro = 3
  p = ggplot(df, aes(x=decile, y=OR )) + geom_pointrange(limits, size=1) + xlab('decile compared to lowest') + ylab('odds ratio (prediction accuracy)') + scale_x_continuous(breaks=1:10) + th2 + 
    theme(aspect.ratio=ar) + theme(legend.key.size=unit(8, 'mm'), legend.position=c(.7,.3), legend.key.height=unit(3,"line")) +
    geom_hline(yintercept=1) + coord_cartesian(ylim=c(0,ymax))
  p
  
}

ccR <- function (cs1, ctrl1, cs2, ctrl2)
{
  r1=cs1/ctrl1
  r2=cs2/ctrl2
  min_ctrl=min(ctrl1, ctrl2)
  min_cs=min(cs1, cs2)
  rho=(min_ctrl*sqrt(r1*r2)+min_cs/sqrt(r1*r2))/sqrt((cs1+ctrl1)*(cs2+ctrl2))
  return(rho)
}

LamMetaPlot <- function(file, title=expression(lambda[meta]), COL=c("blue", "yellow"))
{
  lam=read.table(file, as.is=T)
  mlam=as.matrix(lam)
  mat=matrix(0, nrow=nrow(lam)^2, ncol=3)
  
  cnt=1
  for(i in 1:nrow(mlam))
  {
    for(j in 1:ncol(mlam))
    {
      mat[cnt, 1] = i
      mat[cnt, 2] = j
      mat[cnt, 3] = lam[i,j]
      cnt=cnt+1
    }
  }
  #heatmap plot
  require(ggplot2)
  # first need to reshape data to long form
  require(reshape)
  mat.frame = data.frame(mat)
  names(mat.frame)=c("Cohort1", "Cohort2", "LambdaCohort")
  ggplot(mat.frame, aes(x=Cohort1,y=Cohort2, z= LambdaCohort)) + geom_tile(aes(fill= LambdaCohort)) + scale_fill_gradient(low=COL[1], high=COL[2]) + theme_bw()
}

LamMetaPlotA <- function(file, title=expression(lambda[meta]), COL=c("blue", "yellow"))
{
  #heatmap plot
  require(ggplot2)
  # first need to reshape data to long form
  require(reshape)
  require(zoo)
  lam=read.table(file, as.is=T)
  random_matrix=as.matrix(lam)
  
  
  ## set color representation for specific values of the data distribution
  quantile_range <- quantile(random_matrix, probs = seq(0, 1, 0.02))
  
  ## use http://colorbrewer2.org/ to find optimal divergent color palette (or set own)
  color_palette <- colorRampPalette(c("#3794bf", "#FFFFFF", "#df8640"))(length(quantile_range) - 1)
  
  ## prepare label text (use two adjacent values for range text)
  label_text <- rollapply(round(quantile_range, 2), width = 2, by = 1, FUN = function(i) paste(i, collapse = " : "))
  
  ## discretize matrix; this is the most important step, where for each value we find category of predefined ranges (modify probs argument of quantile to detail the colors)
  mod_mat <- matrix(findInterval(random_matrix, quantile_range, all.inside = TRUE), nrow = nrow(random_matrix))
  
  ## remove background and axis from plot
  theme_change <- theme(
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank()
    #    axis.line = element_blank(),
    #    axis.ticks = element_blank()
    #    axis.text.x = element_blank(),
    #    axis.text.y = element_blank(),
    #    axis.title.x = element_blank(),
    #    axis.title.y = element_blank()
  )
  
  ## output the graphics
  MM=melt(mod_mat)
  idx1=which(MM[,1] <= MM[,2] &MM[,3] > MM[1,3])
  idx2=which(MM[,1] > MM[,2] & MM[,3] < MM[1,3])
  MM1=MM[c(idx1, idx2),]
  
  
  # I only used data from 'random_matrix'
  MM2 = melt(unname(random_matrix))
  # here I store the coordinates of a cell which contain 1 in the original data, so I can see what the value is after the transformation
  midref = MM2[MM2[,1] != MM2[,2] & MM2[,3]==1,][1,1:2]
  MM2[MM[,1] <= MM[,2] & MM[,3] >= MM[1,3], 3] = NA
  MM2[MM[,1] > MM[,2] & MM[,3] < MM[1,3], 3] = NA
  # quantile normalisation
  MM2$value = rank(MM2$value, na.last='keep')
  # transforming back to the original scale; this is cheating somewhat, only the min and max values correspond to the values on the legend
  MM2$value = (MM2$value-1)*(max(random_matrix)-min(random_matrix))/(max(MM2$value, na.rm=T)) + min(random_matrix)
  # look up the mid point after transformation; I don't use this here, because shifting the midpoint away from the middle makes the high and low colours look different
  midp = MM2[MM2[,1] == midref[,1] & MM2[,2] == midref[,2], 3]
  ggplot(MM2, aes(x = X2, y = X1, fill = value)) +
    geom_tile(color = "grey") +
    scale_fill_gradient2(limits=c(min(MM2[,3],na.rm=T), max(MM2[,3],na.rm=T)), low='#3794bf', mid='#FFFFFF', high='#df8640', na.value='white', midpoint=median(MM2$value, na.rm=T), name=expression(paste(lambda[meta]))) +
    # you can use this to see what if looks like if the a value of 1 in the original data is white
    #scale_fill_gradient2(limits=c(min(MM2[,3],na.rm=T), max(MM2[,3],na.rm=T)), low='#3794bf', mid='#FFFFFF', high='#df8640', na.value='white', midpoint=midp) +
    labs(x="Cohort 1", y="Cohort 2") +
    theme_change
}

metaPlot <- function(file, prune=FALSE, cohortFile=NULL, title="Meta cohort plot", sig=TRUE, xlib="Cohort")
{
  if(prune)
  {
    GIANT=matPrune(file, cohortFile)
  }
  else
  {
    GIANT=read.table(file, as.is=F, header=FALSE)
  }
  EleZ=matrix(0, length(GIANT) * (length(GIANT)-1)/2, 4)
  Ele = matrix(0, length(GIANT) * (length(GIANT)-1)/2, 4)
  z_l=qnorm(0.05/2/length(Ele))
  z_h=qnorm(1-0.05/2/length(Ele))

  cnt=0
  for(i in 2:nrow(GIANT))
  {
    for(j in 1:(i-1))
    {
      cnt=cnt+1
      EleZ[cnt,1]=i
      EleZ[cnt,2]=j
      EleZ[cnt,3]=GIANT[j,i]
      Ele[cnt,1]=i
      Ele[cnt,2]=j
      Ele[cnt,3]=GIANT[i,j]
    }
  }

  Ele[,4] = pnorm(abs(Ele[,3] - mean(Ele[,3])) / sd(Ele[,3]))

  C0="gray"
  C1="darkred"
  C2="navy"
  redfunc=colorRampPalette(c(C0, C1), bias=5, space="Lab")
  K=10
  idx=which(Ele[,3] > 0)
  ELE1=Ele[idx,]
  plot(main=title, x=NULL, y=NULL, xlab=xlib, ylab=xlib, xlim=c(-5, nrow(GIANT)+5), ylim=c(-5, nrow(GIANT)+5))
  if(nrow(ELE1) > 0)
  {
    points(ELE1[,1], ELE1[,2], col=redfunc(K)[findInterval(ELE1[,4], seq(min(ELE1[,4]), max(ELE1[,4]), length=K))], pch=20)
    if(sig)
    {
      idxZ=which(EleZ[,3] < 0)
      ELEZh=EleZ[idxZ,]
      idxZh=which(ELEZh[,3] < z_l)
      points(ELEZh[idxZh,1], ELEZh[idxZh,2], col="black", pch=0)    
    }
    legend(x=0, y=-3, legend=c(0, round(max(ELE1[,3]), digits=2)), pch=20, col=c(C0, C1), bty="n", horiz=T)    
  }

  ELE2=Ele[-idx,]
  redfunc=colorRampPalette(c(C0, C2), bias=5, space="Lab")

  if(nrow(ELE2) > 0)
  {
    points(ELE2[,2], ELE2[,1], col=redfunc(K)[findInterval(ELE2[,4], seq(min(ELE2[,4]), max(ELE2[,4]), length=K))], pch=20)
    idxZh=which(ELE2[,3] > z_h)
    points(ELE2[idxZh,1], ELE2[idxZh,2], col="black", pch=0)
    
    if(sig)
    { 
      idxZ=which(EleZ[,3] > 0)
      ELEZl=EleZ[idxZ,]
      idxZl=which(ELEZl[,3] > z_h)
      points(ELEZl[idxZl, 2], ELEZl[idxZl, 1], col="black", pch=0)
    }
    legend(x=0, y=nrow(GIANT)+3, legend=c(round(min(ELE2[,3]), digits=2), 0), pch=20, col=c(C2, C0), bty="n", horiz=T)
  }
}

matPrune <- function(file, cFile=NA)
{
  library(matrixcalc)
  CZmat=read.table(file, as.is=T)
  if(!is.na(cFile))
  {
    cohortFile = as.array(as.matrix(read.table(cFile, as.is=T))[,1])
  }
  mat=as.matrix(CZmat)
  Zmat=mat
  print(dim(Zmat))
  diag(Zmat)= 0

  for(i in 1:nrow(mat))
  {
    for(j in i:ncol(mat))
    {
      mat[i,j]=mat[j,i]
      Zmat[j,i]=Zmat[i,j]
    }
  }

  Z=qnorm(1-0.05/ (nrow(mat)* (nrow(mat)-1)/2))

  CN=0
  IDX=array(1:nrow(mat), nrow(mat))
  while( !is.positive.definite(mat) )
  {
    cnt=array(0,nrow(Zmat))
    for(i in 1:nrow(Zmat))
    {
      if(length(which(Zmat[i,] > Z)) > 0)
      {
        cnt[i] = length(which(Zmat[i,] > Z))
      }
    }
    CN = CN+1
    RevIdx=which(cnt == max(cnt))
    #  print(IDX[RevIdx[1]])
    print(IDX[RevIdx[1]])
    
    IDX = IDX[-RevIdx[1]]
    Zmat=Zmat[-RevIdx[1], -RevIdx[1]]
    mat=mat[-RevIdx[1], -RevIdx[1]]
    if(!is.na(cFile))
    {
      print(cohortFile[RevIdx[1]])
      cohortFile=cohortFile[-RevIdx[1]]      
    }
  }
#  print(IDX)
  pruneZCmat=CZmat[IDX, IDX]
  return(pruneZCmat)
}

matPruneCohort <- function(file, cFile)
{
  library(matrixcalc)
  CZmat=read.table(file, as.is=T)
  cohortFile = as.array(as.matrix(read.table(cFile, as.is=T))[,1])

  mat=as.matrix(CZmat)
  Zmat=mat
  print(dim(Zmat))
  diag(Zmat)= 0
  
  for(i in 1:nrow(mat))
  {
    for(j in i:ncol(mat))
    {
      mat[i,j]=mat[j,i]
      Zmat[j,i]=Zmat[i,j]
    }
  }
  
  Z=qnorm(1-0.05/ (nrow(mat)* (nrow(mat)-1)/2))
  
  CN=0
  IDX=array(1:nrow(mat), nrow(mat))
  while( !is.positive.definite(mat) )
  {
    cnt=array(0,nrow(Zmat))
    for(i in 1:nrow(Zmat))
    {
      if(length(which(Zmat[i,] > Z)) > 0)
      {
        cnt[i] = length(which(Zmat[i,] > Z))
      }
    }
    CN = CN+1
    RevIdx=which(cnt == max(cnt))
    
    IDX = IDX[-RevIdx[1]]
    Zmat=Zmat[-RevIdx[1], -RevIdx[1]]
    mat=mat[-RevIdx[1], -RevIdx[1]]
    cohortFile=cohortFile[-RevIdx[1]]
  }
  #  print(IDX)
  cFile = as.array(as.matrix(read.table(cFile, as.is=T))[,1])
  print(paste("Removed ", nrow(CZmat) - length(IDX), "cohorts."))

  return(cohortFile)
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

genomeReader <- function(root)
{
  gFile = paste0(root, ".genome")
  genome = read.table(gFile, as.is=T, header=T)

  idFile = read.table(paste0(root, ".fam"), as.is=T)

  mat=matrix(0, nrow(idFile), nrow(idFile))
  diag(mat)=1
  idx1=1
  L=nrow(idFile)
  for(i in 1:(nrow(idFile)-1))
  {
    idx2 = idx1 + (L-i) -1
    mat[i, (i+1):L] = genome[idx1:idx2, 10]
    mat[(i+1):L, i] = genome[idx1:idx2, 10]
    idx1 = idx1 + (L-i)
  }
  return(mat)
}

grmBinReader <- function(root)
{
  binID=paste0(root, ".grm.id")
  idF=file(binID, "r")
  grmID=read.table(binID, as.is=T)

  binGRM=paste0(root, ".grm.bin")
  binF=file(binGRM, "rb")
  grm=readBin(binF, n=nrow(grmID)*(nrow(grmID)+1)/2, what=numeric(0), size=4)
  mat=matrix(0, nrow=nrow(grmID), ncol=nrow(grmID))
  
  c1=0
  c2=0
  gmat=matrix(0, nrow(grmID), nrow(grmID))
  for(i in 1:nrow(grmID))
  {
    c1= i * (i-1)/2 +1
    c2= i * (i+1)/2
    gmat[i, 1:i] = grm[c1:c2]
    gmat[1:i, i] = grm[c1:c2]
  }
  close(idF)
  close(binF)
  return(gmat)
}

Hong23 <-function(K, case, control) {
  z=dnorm(qnorm(K))
  p=case/(case+control)
  H23=K*(1-K)*K*(1-K) / (z^2 * p * (1-p))
  return(H23)
}

Hong23h2 <-function(h2o, K, case, control) {
  z=dnorm(qnorm(K))
  p=case/(case+control)
  H23=h2o * K*(1-K)*K*(1-K) / (z^2 * p * (1-p))
  return(H23)
}

Fst <- function(geno, group) {
  #geno is a matrix for genotypes, each row represents an individuals
  #each col represents a locus
  
  GM=apply(geno, 2, mean)/2
  VT=2*GM*(1-GM)
  
  w=table(group)
  cl=names(w)
  wt=w/length(group)
  vt=matrix(0, length(cl), ncol(geno))
  
  for(i in 1:length(cl)) {
    idx = which(group==cl[i])
    g = apply(geno[idx,], 2, mean)/2
    vt[i,] = 2*g*(1-g)*wt[i]
  }
  VS=apply(vt, 2, sum)
  fst = 1-VS/VT
  return(fst)
}

OverLapIndex <- function(s1, s2) {
  all=c(s2, s1)
  flag=duplicated(all)
  flag=flag[(1+length(s2)):length(all)]
  idx=which(flag)
  return(idx)
}

manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
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

miamiPlot <- function(dataframe, P1=NULL, P2=NULL, AnoCol="red", AnoCol1="blue", AnoCol2="yellow",colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, annotate1=NULL, annotate2=NULL, ...) {
  if(is.null(P1))
  {
    P1 = "P"
  }
  pidx1=which(colnames(dataframe) == P1)
  
  if(is.null(P2))
  {
    P2 = "P2"
  }
  pidx2=which(colnames(dataframe) == P2)
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & P1 %in% names(d) & P2 %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d[,pidx1])
  d$logp2 = log10(d[,pidx2])
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  if (min(d$logp2) > -8) {
    ymin = -8
  } else {
    ymin = floor(min(d$logp2))
  }
  
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
    with(d, plot(main=title, pos, logp, ylim=c(-1*max(abs(ymin),abs(ymax)), max(abs(ymin),abs(ymax))), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  } else {
    with(d, plot(main=title, pos, logp, ylim=c(-1*max(abs(ymin),abs(ymax)), max(abs(ymin),abs(ymax))), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      with(d[d$CHR==i, ],points(pos, logp2, col=colors[icol], ...))
      icol=icol+1
    }
  }
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col=AnoCol, pch=16, ...))
    with(d.annotate, points(pos, logp2, col=AnoCol, pch=16, ...))
#    with(d.annotate, points(pos, logp, col="grey", cex=1, pch=1))
#    with(d.annotate, points(pos, logp2, col="grey", cex=1, pch=1))
    
  }

  if (!is.null(annotate1)) {
    d.annotate=d[which(d$SNP %in% annotate1), ]
    with(d.annotate, points(pos, logp, col=AnoCol1, pch=16))
#    with(d.annotate, points(pos, logp, col="grey", cex=1, pch=1))
    
  }

  if (!is.null(annotate2)) {
    d.annotate=d[which(d$SNP %in% annotate2), ]
    with(d.annotate, points(pos, logp2, col=AnoCol2, pch=16))
#    with(d.annotate, points(pos, logp2, col="grey", cex=1, pch=1))
    
  }
  
  #  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (!is.null(genomewideline)) {
    abline(h=genomewideline, col="gray")
    abline(h=-1*genomewideline, col="gray")
    
  } else {
    abline(h=-log10(0.05/nrow(d)), col="gray")
    abline(h=log10(0.05/nrow(d)), col="gray")
  }
  abline(h=0, col="white", lwd=2, lty=2)
}


manhattanTrans <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, ...) {
  
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
    with(d, plot(bty='n', main=title, y=pos, x=logp, xlim=c(0,ymax), xlab=expression(-log[10](italic(p))), ylab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }  else {
    with(d, plot(bty='n', main=title, y=pos, x=logp, xlim=c(0,ymax), xlab=expression(-log[10](italic(p))), ylab="Chromosome", yaxt="n", type="n", ...))
    axis(2, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(y=pos, x=logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(y=pos, x=logp, col="green3", ...))
  }
  #  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (!is.null(genomewideline)) {
    abline(v=genomewideline, col="gray")
  } else {
    abline(v=-log10(0.05/nrow(d)), col="gray")    
  }
}

manhattanP <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=NULL, title="", annotate=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$P))
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
    with(d, plot(main=title, pos, P, ylim=c(0,ymax), ylab=expression(italic(p)), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }  else {
    with(d, plot(main=title, pos, P, ylim=c(0,ymax), ylab=expression(italic(p)), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, P, col=colors[icol], ...))
      icol=icol+1
    }
  }
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, P, col="green3", ...))
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
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(italic("F")[italic("ST")])), xlab=paste("Chromosome",unique(d$CHR),"position"), main=title,...))
  }  else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(paste(italic("F")[italic("ST")])), xlab="Chromosome", xaxt="n", type="n", main=title,...))
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
