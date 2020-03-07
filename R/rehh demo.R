library(rehh)
make.example.files()
head(read.table("map.inp"))

hap<-data2haplohh(hap_file="bta12_cgu.hap",map_file="map.inp",
                  recode.allele=TRUE,chr.name=12)

#example haplohh object (280 haplotypes, 1424 SNPs) see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12)
#computing EHH statistics for the focal SNP at position 456
#which display a strong signal of selection
res.ehh<-calc_ehh(haplohh_cgu_bta12,mrk=456)

res.ehh$ehh[1:2,454:458]
res.ehh$nhaplo_eval[1:2,454:458]
res.ehh$freq_all1
res.ehh$ihh

#computing EHH statistics for the focal SNP at position 456
#which display a strong signal of selection
res.ehhs<-calc_ehhs(haplohh_cgu_bta12,mrk=456)

res.ehhs$EHHS_Sabeti_et_al_2007[453:459]
res.ehhs$EHHS_Tang_et_al_2007[453:459]
res.ehhs$nhaplo_eval[453:459]
res.ehhs$IES_Tang_et_al_2007
res.ehhs$IES_Sabeti_et_al_2007

res.scan<-scan_hh(haplohh_cgu_bta12)
dim(res.scan)

#NA
head(res.scan)
system.time(res.scan<-scan_hh(haplohh_cgu_bta12))

foo<-function(haplo){
  res.ihh=res.ies=matrix(0,haplo@nsnp,2)
  for(i in 1:length(haplo@position)){
    res.ihh[i,]=calc_ehh(haplo,mrk=i,plotehh=FALSE)$ihh
    tmp=calc_ehhs(haplo,mrk=i,plotehhs=FALSE)
    res.ies[i,1]=tmp$IES_Tang_et_al_2007
    res.ies[i,2]=tmp$IES_Sabeti_et_al_2007
  }
  list(res.ies=res.ies,res.ihh=res.ihh)
}
system.time(res.scan2<-foo(haplohh_cgu_bta12))

#NA
sum(res.scan2$res.ihh[,1]!=res.scan[,4]) + sum(res.scan2$res.ihh[,2]!=res.scan[,5]) +
  sum(res.scan2$res.ies[,1]!=res.scan[,6]) + sum(res.scan2$res.ies[,2]!=res.scan[,7])

#error
for(i in 1:29){
  hap_file=paste("hap_chr_",i,".pop1",sep="")
  data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
  res<-scan_hh(data)
  if(i==1){wg.res<-res}else{wg.res<-rbind(wg.res,res)}
}
wg.ihs<-ihh2ihs(wg.res)

data(wgscan.cgu)
ihs.cgu<-ihh2ihs(wgscan.cgu)
head(ihs.cgu$iHS)
head(ihs.cgu$frequency.class)

ihsplot(ihs.cgu,plot.pval=TRUE,ylim.scan=2,main="iHS (CGU cattle breed)")

#error
for(i in 1:29){
  hap_file=paste("hap_chr_",i,".pop1",sep="")
  data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
  res<-scan_hh(data)
  if(i==1){wg.res.pop1<-res}else{wg.res.pop1<-rbind(wg.res.pop1,res)}
  hap_file=paste("hap_chr_",i,".pop2",sep="")
  data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
  res<-scan_hh(data)
  if(i==1){wg.res.pop2<-res}else{wg.res.pop2<-rbind(wg.res.pop2,res)}
}
wg.rsb<-ies2rsb(wg.res.pop1,wg.res.pop2)

data(wgscan.cgu) ; data(wgscan.eut)
cguVSeut.rsb<-ies2rsb(wgscan.cgu,wgscan.eut,"CGU","EUT")
head(cguVSeut.rsb)

rsbplot(cguVSeut.rsb,plot.pval=TRUE)

#error
for(i in 1:29){
  hap_file=paste("hap_chr_",i,".pop1",sep="")
  data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
  res<-scan_hh(data)
  if(i==1){wg.res.pop1<-res}else{wg.res.pop1<-rbind(wg.res.pop1,res)}
  hap_file=paste("hap_chr_",i,".pop2",sep="")
  data<-data2haplohh(hap_file="hap_file","snp.info",chr.name=i)
  res<-scan_hh(data)
  if(i==1){wg.res.pop2<-res}else{wg.res.pop2<-rbind(wg.res.pop2,res)}
}
wg.xpehh<-ies2xpehh(wg.res.pop1,wg.res.pop2)

data(wgscan.cgu) ; data(wgscan.eut)
cguVSeut.xpehh<-ies2xpehh(wgscan.cgu,wgscan.eut,"CGU","EUT")
head(cguVSeut.xpehh)

xpehhplot(cguVSeut.xpehh,plot.pval=TRUE)
plot(cguVSeut.rsb[,3],cguVSeut.xpehh[,3],xlab="Rsb",ylab="XPEHH",pch=16,
     cex=0.5,cex.lab=0.75,cex.axis=0.75)
abline(a=0,b=1,lty=2)
distribplot(ihs.cgu$iHS[,3],xlab="iHS")

data(haplohh_cgu_bta12)
layout(matrix(1:2,2,1))
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=1,nmrk_l=20,nmrk_r=20,
                    main="Bifurcation diagram (RXFP2 SNP on BTA12): Ancestral Allele")
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=2,nmrk_l=20,nmrk_r=20,
                    main="Bifurcation diagram (RXFP2 SNP on BTA12): Derived Allele")