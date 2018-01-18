args=commandArgs(TRUE)
root=args[1]
rep=args[2]

rsub='Rscript /clusterdata/gc5k/bin/rsub.R'
for(i in 1:rep)
{
  CMD=paste(rsub, "Rscript eNet.R", root, i)
  CMD=paste0(CMD, " @", root, "_", i, "eNet", " %30G")
  print(CMD)
  system(CMD)
}
