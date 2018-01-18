RevSNP=c()

frq=read.table("IBD_UC_tmp.frq", as.is=T, header=T)
frqCut=0.01
frqidx=which(frq$MAF > frqCut)
png("MAF.png", width=800, height=800)
layout(matrix(1:2, 1, 2))
hist(frq$MAF, main="Distribution of MAF for all SNPs", xlab="MAF")
hist(frq$MAF[frqidx], main="Distribution of MAF (MAF > 0.01)", xlab="MAF")
dev.off()
DelSNP1=frq$SNP[-frqidx]
write.table(frq$SNP[-frqidx], paste0("MAF", frqCut, ".txt"), row.names=F, col.names=F, quote=F)

hdcut=0.000001
hardy=read.table("IBD_UC_tmp.hwe", as.is=T, header=T)
hardyU=hardy[hardy$TEST=="UNAFF",]
hdIdx=which(hardyU$P < hdcut)
write.table(hardyU$SNP[hdIdx], paste0("hwe",hdcut, ".txt"), row.names=F, col.names=F, quote=F)
DelSNP2=(hardyU$SNP[hdIdx])

RevSNP=union(DelSNP1, DelSNP2)

miss=0.01
MS=read.table("IBD_UC_tmp.lmiss", as.is=T, header=T)
lociIdx=which(MS$F_MISS > 0.05)
RevSNP=union(MS$SNP[lociIdx], RevSNP)

write.table(RevSNP, "RevSNP.txt", row.names=F, col.names=F, quote=F)

png("Inbreeding.png", width=800, height=800)
het=read.table("IBD_UC_tmp.het", as.is=T, header=T)
hist(het$F, main="Inbreeding (Plink --het)", xlab="Inbreeding")
dev.off
