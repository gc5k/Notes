make_shell<-function(pbs_file, job_name, cmd, mem) 
{

  write(paste("#PBS -N ", job_name), pbs_file, append=FALSE)
  write("#PBS -l walltime=200:00:00", pbs_file, append=TRUE)
  write(paste("#PBS -l vmem=", mem, sep=""), pbs_file, append=TRUE)
  write("cd $PBS_O_WORKDIR", pbs_file, append=TRUE)
  write(cmd, pbs_file, append=TRUE)

}

gear='java -Xmx5G -jar /public/home/gc5k/bin/gear.jar'
gcta='/public/home/gc5k/bin/gcta64'
plink='/public/home/gc5k/bin/plink-1.07-x86_64/plink'
###################################################################

args=commandArgs(TRUE)
CMD=c()
nm=c()
mem="30G"
runit=TRUE

for (i in 1:length(args)) {
  
  flag = TRUE
  breakFlag = FALSE  
  cmdbit=args[i]
  if (substr(args[i], nchar(args[i]), nchar(args[i])) ==",") {
    breakFlag = TRUE
    cmdbit = substr(args[i], 1, nchar(args[i])-1)
  }
  
  if (substr(cmdbit, 1, 1) == "@") {
    nm=substr(cmdbit, 2, nchar(cmdbit))
    flag = FALSE
  } else if (substr(cmdbit, 1, 1) == "%") {
    mem=substr(cmdbit,2, nchar(cmdbit))
    flag = FALSE
  } else if (substr(cmdbit, 1, 2) == "TT" && nchar(cmdbit)==2) {
    runit=FALSE
    flag = FALSE
  } else if (cmdbit == 'gear' ) {
    cmd=gear
  } else if (cmdbit == 'gcta') {
    cmd=gcta
  } else if (cmdbit == 'plink') {
    cmd=plink
  } else {
    cmd=cmdbit
  }
  
  if(flag) {
    CMD=paste(CMD, cmd)
    if(breakFlag) {
      CMD=paste(CMD, "\n")
    }
  }
}

print(CMD)
print(nm)
if (length(nm)>0) {
  pbs_file=paste(nm, ".sh", sep="")
  job_name=paste(nm, sep="")
} else {
  pbs_file=paste(args[1], args[length(args)], ".sh", sep="-")
  job_name=paste(args[1], args[length(args)], sep="-")
}

make_shell(pbs_file, job_name, CMD, mem)
sys_qsub=paste("qsub", pbs_file)
if(runit) {
  system(sys_qsub)
}
