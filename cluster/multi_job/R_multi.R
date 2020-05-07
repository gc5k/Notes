
for(i in 1:5) {
  myDir=paste0("T",i) #dir name
  dir.create(myDir) #create dir
  setwd(myDir) #shift dir


#####make script
  str = paste("#PBS -N Multi_show", i, "\n#PBS -l nodes=1:ppn=1\n#PBS -l mem=4G\n#PBS -l walltime=10:00:00\n#PBS -j oe\ncd $PBS_O_WORKDIR\n") #parameter for shell job
  str = paste0(str, "Rscript ../R_write.R ", i, " T", i, ".txt")

  write.table(str, paste0("Mjob", i, ".sh"), row.names = F, col.names = F, quote = F)
  pbs = paste0("qsub Mjob", i, ".sh")
  system(pbs)
  setwd("..")
}
