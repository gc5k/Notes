make_shell<-function(pbs_file, job_name, cmd, mem, clean) {
  
  write("#$ -cwd", pbs_file, append=FALSE)
  write(paste("#$ -l vf=", mem, sep=""), pbs_file, append=TRUE)
  write(paste("#$ -l h_vmem=", mem, sep=""), pbs_file, append=TRUE)
  write(paste("#$ -N", job_name), pbs_file, append=TRUE)
  
  write("#$ -m eas", pbs_file, append=TRUE)
  write("#$ -M chen.guobo@foxmail.com", pbs_file, append=TRUE)
  write(cmd, pbs_file, append=TRUE)
  
}

gmdr='java -Xmx5G -jar /clusterdata/gc5k/bin/gmdr.jar'
gcta='/clusterdata/gc5k/bin/gcta64'
gcta_test='/clusterdata/gc5k/bin/gcta64_test'
gctaM='/ibscratch/wrayvisscher/jyang/bin/gcta64_test'
plink='/clusterdata/gc5k/bin/plink-1.07-x86_64/plink'
polygenic='java -jar /clusterdata/gc5k/bin/polygenic.jar'
HE='/clusterdata/gc5k/bin/HE.jar'
gear='java -Xmx10G -jar /clusterdata/gc5k/bin/gear.jar'
###################################################################

args=commandArgs(TRUE)
CMD=c()
nm=c()
mem="10G"
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
  } else if (cmdbit == 'gmdr') {
    cmd=gmdr
  } else if (cmdbit == 'HE' ) {
    cmd=HE
  } else if (cmdbit == 'gear' ) {
    cmd=gear
  }else if (cmdbit == 'gctaM' ) {
    cmd=gctaM
  } else if (cmdbit == 'gcta_test' ) {
    cmd=gcta_test
  } else if (cmdbit == 'gcta') {
    cmd=gcta
  } else if (cmdbit == 'plink') {
    cmd=plink
  } else if (cmdbit == 'polygenic') {
    cmd=polygenic
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
