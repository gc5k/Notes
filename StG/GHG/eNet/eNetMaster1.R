args=commandArgs(TRUE)
geno=args[1]
pheno=args[2]
trt=args[3]
rep=as.numeric(args[4])
type=args[5]

rsub='Rscript /clusterdata/gc5k/bin/rsub.R'
Str=unlist(strsplit(geno, "/"))
root=Str[length(Str)]
for(i in 1:rep)
{
  if(type == "B")
  {
    pheIdx=3
  }
  else
  {
    pheIdx=3+i
  }
  trtIdx=3+i
  CMD=paste(rsub, "Rscript eNet1.R", geno, pheno, pheIdx, trt, trtIdx, i, type)
  CMD=paste0(CMD, " @", root, ".", type, ".", i, " %30G")
  print(CMD)
  system(CMD)
}
