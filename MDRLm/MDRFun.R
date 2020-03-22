MultiCode <- function(geno, comb)
{
  gCode=array(0, dim=nrow(geno))
  for(j in 1:nrow(geno))
  {
    for(k in 1:(length(comb)) )
    {
      if(k == 1)
      {
        gCode[j] = geno[j, comb[k]]
      }  else {
        gCode[j] = paste0(gCode[j], geno[j, comb[k]])
      }
    }
  }
  return(gCode)
}

GetScheme <- function(gCode, Y)
{
  localY=scale(Y)
  code=names(table(gCode))
  scheme=array(0, dim=length(code))
  for(j in 1:length(code))
  {
    idx=which(gCode == code[j])
    m=mean(localY[idx])
    scheme[j]=ifelse(m > 0, 1, 0)
  }
  return (list(scheme, code))
}

MDR <- function(scheme, MG, gCode)
{
  xTr=array(rbinom(length(gCode), 1, 0.5), dim=length(gCode))
  N=0
  for(j in 1:length(scheme))
  {
    idxTr=which(gCode == MG[j])
    if(length(idxTr) > 0)
    {
      xTr[idxTr] = scheme[j]
      N = N+length(idxTr)      
    }
  }
  print(N)
  return(xTr)
}

LossFun <- function(gCode, gCodeT)
{
  uT=unique(gCode)
  uTr=unique(gCodeT)
  intU=intersect(uT, uTr)

  
  LF=0
  if(length(intU) > 0)
  {
    N=0
    for(i in 1:length(intU))
    {
      idx=which(gCode == intU[i])
      N = N + length(idx)
    }
    LF=N/length(gCode)
  }
  return(LF)
}

MDRScoreTrain <- function(scheme, MG, gCode, Y)
{
  TP=0
  FP=0

  FN=0
  TN=0

  localY=scale(Y)
  for(i in 1:length(scheme))
  {
    idx=which(gCode==MG[i])
    iP=which(localY[idx] > 0)
    iN=which(localY[idx] < 0)
    if(scheme[i] == 1)
    {
      if(length(iP) > 0)
      {
        TP = TP + sum(localY[idx[iP]])
      }
      if(length(iN) > 0)
      {
        FN = FN + sum(localY[idx[iN]])
      }
    }
    if(scheme[i] == 0 )
    {
      if(length(iN) > 0)
      {
        TN = TN + sum(localY[idx[iN]])
      }
      if(length(iP) > 0)
      {
        FP = FP + sum(localY[idx[iP]])
      }
    }
  }
  print(paste(TP, TN, FP, FN))
  A=(TP+abs(TN))/(TP + abs(TN) + FP + abs(FN))
  return(A)
}

MDRScoreTest <- function(scheme, MG, gCode, Y)
{
  TP=0
  FP=0

  FN=0
  TN=0

  localY=scale(Y)
  for(i in 1:length(scheme))
  {
    idx=which(gCode==MG[i])
    if(length(idx)>0)
    {
      iP=which(localY[idx] > 0)
      iN=which(localY[idx] < 0)
      if(scheme[i] == 1)
      {
        if(length(iP) > 0)
        {
          TP = TP + sum(localY[idx[iP]])
        }
        if(length(iN) > 0)
        {
          FN = FN + sum(localY[idx[iN]])
        }
      }
      if(scheme[i] == 0 )
      {
        if(length(iN) > 0)
        {
          TN = TN + sum(localY[idx[iN]])
        }
        if(length(iP) > 0)
        {
          FP = FP + sum(localY[idx[iP]])
        }
      }
    }
  }
  
  unknownMG=setdiff(unique(gCode), MG)
  ns=0
  if(length(unknownMG) > 0)
  {
    for(i in 1:length(unknownMG))
    {
      idx=which(gCodeT == unknownMG[i])
      ns = ns+sum(Y[idx])
    }
  }
  print(paste(TP, TN, FP, FN, ns))
  A=(TP+abs(TN))/(TP + abs(TN) + FP + abs(FN) + abs(ns))
  return(A)
}

CalLD <- function(freq, dprime)
{
  if(length(freq) == length(dprime))
  {
    ld=array(dim=length(dprime)-1)
  }
  else if(length(dprime) == (length(freq) -1) )
  {
    ld=array(dim=length(dprime))
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

SimuQt <- function(frq, size, locIdx=c(1), gFun=c(1), intercept=0, gbeta=1, cBeta=NULL, rsq=0.5) {
  bs=10^c((length(locIdx)-1):0)
  gMat=matrix(0, size, length(frq))
  for(i in 1:nrow(gMat)) {
    gMat[i,]=rbinom(length(frq), 2, frq)
  }

  gStack=gMat[,locIdx]%*% matrix(bs, nrow=length(bs), 1)
  bv=(gStack %in% gFun)*gbeta
  cMat=matrix(1, nrow=nrow(gMat), 1)
  cB=c(intercept)
  if(!is.null(cBeta)) {
    cMat=cbind(cMat, matrix(rnorm(length(bv)*length(cBeta)), length(bv), length(cBeta)))
    cB=c(cB, cBeta)
  }

  bv=bv+cMat%*%matrix(cB, nrow=length(cB), 1)
  y=bv+rnorm(length(bv), sd=sqrt(var(bv)/rsq * (1-rsq)))
  dat=list("gMat"=gMat, "y"=y, "X"=cMat)
}

SimuCC <- function(frq, csN, ctrlN, locIdx=c(1), gFun=c(1), intercept=-4.5, gbeta=1, cBeta=NULL) {
  bs=10^c((length(locIdx)-1):0)
  csn=0
  ctrln=0

  Y=array(0, dim=csN+ctrlN)
  gMat=matrix(0, nrow=csN+ctrlN, ncol=length(frq))
  cMat=matrix(1, nrow=csN+ctrlN, 1)
  cB=c(intercept)
  if(!is.null(cBeta)) {
    cMat=cbind(cMat, matrix(0, csN+ctrlN, length(cBeta)))
    cB=c(cB, cBeta)
  }

  cnt=0
  cnTotal=0
  while((cnt)<(csN+ctrlN)) {
    cnTotal=cnTotal+1
    g=rbinom(length(frq), 2, frq)
    gstack=sum(g[locIdx] * bs)
    
    bv=intercept
    if(!is.null(cBeta)) {
      x=rnorm(length(cBeta))
      bv=bv+x*cBeta
    }

    if(gstack %in% gFun) {
      bv=bv+gbeta
    }
    pv=exp(bv)/(1+exp(bv)) #logistic 
    ind=rbinom(1, 1, pv)

    if(ind == 1) { #if case
      if(csn < csN) {
        csn=csn+1
        cnt=cnt+1
        if(!is.null(cBeta)) {
          cMat[csn,2:ncol(cMat)]=x
          gMat[csn,]=g
        }
        Y[csn]=1
      }
    } else { #control
      if(ctrln < ctrlN) {
        ctrln=ctrln+1
        cnt=cnt+1
        if(!is.null(cBeta)) {
          cMat[ctrln+csN,2:ncol(cMat)]=x
          gMat[ctrln+csN,]=g
        }
        Y[ctrlN+csN]=0
      }
    }
  }
  dat=list("y"=Y, "X"=cMat, "gMat"=gMat, "cT"=cnTotal)
}
