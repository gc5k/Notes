#########
#D1.3-1
#LD
#########
source("shotgun.R")
layout(matrix(1:2, 1, 2))
fq=runif(100, 0.1, 0.5)
dp=runif(99, 0.9, 1)
Gdp=GenerateGenoDprime(fq, dp, 1000)
cdp=cor(Gdp)
heatmap(cdp, Colv = NA, Rowv = NA, symm = T, col=topo.colors(32, 0.3))

d=Dprime2LD(fq, dp)
Gd=GenerateGeno(fq, d, 1000)
cd=cor(Gd)
heatmap(cd, Colv = NA, Rowv = NA, symm = T, col=topo.colors(32, 0.3))


################hapmap
ceu=read.table("../HM/ceu_r.raw", as.is = T, header = T)
cceu=cor(ceu[,-c(1:6)])
heatmap(cceu, Colv = NA, Rowv = NA, labRow = NA, labCol = NA, col=topo.colors(32, 0.3))

chb=read.table("../HM/chb_r.raw", as.is = T, header = T, col=topo.colors(32, 0.3))
cchb=cor(chb[,-c(1:6)])
heatmap(cchb, Colv = NA, Rowv = NA, labRow = NA, labCol = NA, col=topo.colors(32, 0.3))

yri=read.table("../HM/yri_r.raw", as.is = T, header = T, col=topo.colors(32, 0.3))
cyri=cor(yri[,-c(1:6)])
heatmap(cyri, Colv = NA, Rowv = NA, labRow = NA, labCol = NA, col=topo.colors(32, 0.3))

###############hapmap freq chr 22
ceufq=read.table("../HM/ceu_22.frq", as.is = T, header = T)
chbfq=read.table("../HM/chb_22.frq", as.is = T, header = T)
yrifq=read.table("../HM/yri_22.frq", as.is = T, header = T)

layout(matrix(1:3, 1, 3))
hist(ceufq$MAF, xlab="CEU", main = "CEU MAF CHR 22")
hist(chbfq$MAF, xlab="CHB", main = "CHB MAF CHR 22")
hist(yrifq$MAF, xlab="YRI", main = "YRI MAF CHR 22")

##############hapmap freq genome
ceuGfq=read.table("../HM/ceu_fq.frq", as.is = T, header = T)
chbGfq=read.table("../HM/chb_fq.frq", as.is = T, header = T)
yriGfq=read.table("../HM/yri_fq.frq", as.is = T, header = T)

layout(matrix(1:3, 1, 3))
hist(ceuGfq$MAF, xlab="CEU", main = "CEU MAF", breaks = 50)
hist(chbGfq$MAF, xlab="CHB", main = "CHB MAF", breaks = 50)
hist(yriGfq$MAF, xlab="YRI", main = "YRI MAF", breaks = 50)
