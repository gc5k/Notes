FN="Arab295"
PC=5

source("EigenGWAS_Friends.R")

gwas=read.table("Arab452SetMAFclean_Naive_1.assoc.linear", as.is = T, header = T)
####GRM stats
#grmStats(FN)
grmStats(FN, "pdf") #save as pdf automatically

############PC plot
#pcPlot(FN)
pcPlot(FN, "pdf")
pcMatPlot(FN, 5, COL="green", ma=0.2)

####EigenValue plot
#EigenValuePlot(FN, PC)
EigenQQPlot(FN, 1)
EigenValuePlot(FN, PC, "pdf")
DeepEigenValuePlot(FN, 1, seq(0.1, 0.9, 0.05))
###EigenGWAS plot
#EigenGWASPlot(FN, 1)

for(i in 1:PC)
{
  EigenGWASPlot(FN, i, "pdf")
  SWEigenGWASPlot(FN, i, 10, "pdf") 
  #sliding window, @parameters, 1 file name, 2 pc, 3 kb, 4 pdf or not
}

####pheno eigen
SWPhenoEigenGWASPlot("Arab452SetMAFclean_Naive_1.assoc.linear", FN, 10)

