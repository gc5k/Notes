FN="Arab295"
PC=5

source("EigenGWAS_Friends.R")

####GRM stats
#grmStats(FN)
grmStats(FN, "pdf") #save as pdf automatically

############PC plot
#pcPlot(FN)
pcPlot(FN, "pdf")

####EigenValue plot
#EigenValuePlot(FN, PC)
EigenValuePlot(FN, PC, "pdf")

###EigenGWAS plot
#EigenGWASPlot(FN, 1)

for(i in 1:PC)
{
  EigenGWASPlot(FN, i, "pdf")
}
