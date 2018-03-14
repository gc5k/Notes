install.packages("~/git/EigenGWASFriends_0.1.0.tar.gz", repos = NULL, type = "source")
library("EigenGWASFriends")

################
# Day 1PM Arab
################
FN="arab"
RunEigenGWAS(FN, PC = 5, inbred = T)

pcMatPlot(FN, c(1:5))

EigenValuePlot(FN, 5)

EigenGWASPlot(FN, pc = 1)

SWEigenGWASPlot(FN, pc = 1, kb = 10)

################
# Day 1PM Arab
################
FN2="ceu_tsi"
RunEigenGWAS(FN2, PC = 5, inbred = T)

pcMatPlot(FN2, c(1:5))

EigenValuePlot(FN2, 5)

EigenGWASPlot(FN2, pc = 1)

SWEigenGWASPlot(FN2, pc = 1, kb = 10)
