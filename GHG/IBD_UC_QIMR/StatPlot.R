
#############
BinFileName="IBD_UC_clean1.grm.bin"
BinFile=file(BinFileName, "rb")
id=read.table("IBD_UC_clean1.grm.id")
grm=readBin(BinFile, n=nrow(id)*(nrow(id)+1)/2, what=numeric(0), size=4)
i=seq(1, nrow(id))
dg=i*(i+1)/2
grm1=grm[-dg]
png("GRM.png", width=800, height=800)
hist(grm1, main="Genetic relatedness", xlab="genetic relatedness")
dev.off()

source("~/bin/MyLib/shotgun.R")
png("GWAS.png", height=600, width=800)
assoc=read.table("IBD_UC_clean1.assoc", as.is=T, header=T)
idx=which(assoc$CHR >= 1 & assoc$CHR <=22)
assoc=assoc[idx,]
manhattan(assoc, title="QIMR IBD UC GWAS")
dev.off()

idx=which(-log10(assoc$P) > 7 & assoc$CHR!=6)
gwasAb=assoc[idx,]
write.table(gwasAb, "StrangeGWASSum.txt", row.names=F, col.names=T, quote=F)
write.table(gwasAb[,2], "StrangeGWASLoci.txt", row.names=F, col.names=F, quote=F)
