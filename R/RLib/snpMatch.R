Flip<-function(a) {
  f=a
  if(a=="A") {
    f="T"
  } else if(a=="C") {
    f="G"
  } else if(a=="G") {
    f="C"
  } else {
    f= "A"
  }
  f
}

Confusion<-function(a,b) {
  snp1=a
  snp2=b
  if(snp1>snp2) {
    t=snp1
    snp1=snp2
    snp2=t
  }
  flag = F
  if(snp1=="A" && snp2=="T") {
    flag = T
  }
  if(snp1=="C" && snp2=="G") {
    flag = T
  }
  flag
}

IsBiallelic<-function(La1, La2, Lb1, Lb2) {
  a1=La1;
  a2=La2;
  if(a1>a2) {
    t=a1
    a1=a2
    a2=t
  }
  L1=paste(a1,a2, sep="")
  b1=Lb1
  b2=Lb2
  if(b1>b2) {
    t=b1
    b1=b2
    b2=t
  }
  L2=paste(b1,b2, sep="")
  #both biallelic
  if(a1!=a2 && b1 != b2) {
    flag=F

    if(L1=="AC") {
      if(L2=="AC" || L2=="GT") {
        flag=T
      }
    }
    if(L1=="AG") {
      if(L2=="AG" || L2=="CT") {
        flag=T
      }
    }
    if(L1=="AT") {
      if(L2=="AT") {
        flag=T
      }
    }
    if(L1=="CG") {
      if(L2=="CG") {
        flag=T
      }
    }
    if(L1=="CT") {
      if(L2=="CT" || L2=="AG") {
        flag=T
      }
    }
    if(L1=="GT") {
      if(L2=="GT" || L2=="AC") {
        flag=T
      }
    }
  } else if(a1!=a2 && b1==b2) { #L1 biallelic
    flag == T
    if(L1=="AT") {
      if(b1=="C" || b1=="G") {
        flag=F
      }
    }
    if(L1=="CG") {
      if(b1=="A" || b1=="T") {
        flag=F
      }
    }
  } else if(a1==a2 && b1!=b2) { #L2 biallelic
    flag = T
    if(L2=="AT") {
      if(a1=="C" || a1=="G") {
        flag=F
      }
    } else if(L2=="CG") {
      if(a1=="A" || a1=="T") {
        flag=F
      }
    }
  } else if(a1==a2 && b1==b2) { # L2 both monomorphic
    flag=T
    if(a1=="A") {
      if(b1=="C" || b1=="G") {
        flag=F
      }
    } else if(a1=="C") {
      if(b1=="A" || b1=="T") {
        flag=F
      }
    } else if(a1=="G") {
      if(b1=="A" || b1=="T") {
        flag=F
      }
    } else if(a1=="T") {
      if(b1=="C" || b1=="G") {
        flag=F
      }
    }
  }
  flag
}
