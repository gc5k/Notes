#PBS -N helloR
#PBS -l nodes=1:ppn=1
#PBS -l mem=65G
#PBS -l walltime=500:00:00 
#PBS -j oe
cd $PBS_O_WORKDIR 

Rscript helloWorld.R


