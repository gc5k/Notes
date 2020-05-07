arg=commandArgs(TRUE)

Up=as.numeric(arg[1])
File=arg[2]
write.table(seq(1:Up), File, row.names=F, col.names=F, quote=F)
