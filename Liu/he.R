he=read.table("HE_Zero_QTL_LE_T_Standard_10000_simulation(1).txt", header=T, as.is=T)
he$P=pnorm(he$Tvalue)
layout(matrix(1:2, 1, 2))
plot(he$B, he$X.B_variance...1.2., xlab="B", ylab="SE")
plot(-log10(sort(runif(nrow(he)))), -log10(sort(he$P)), xlab="Expected P", ylab="Observed P")
abline(a=0, b=1)


he1=read.table("HE_Zero_QTL_LE_T_Standard_plus_10000_simulation(1).txt", header=T, as.is=T)
he1$P=pnorm(he1$Tvalue)
layout(matrix(1:2, 1, 2))
plot(he1$B, he1$X.B_variance...1.2., xlab="B", ylab="SE")
plot(-log10(sort(runif(nrow(he1)))), -log10(sort(he1$P)), xlab="Expected P", ylab="Observed P")
abline(a=0, b=1)

