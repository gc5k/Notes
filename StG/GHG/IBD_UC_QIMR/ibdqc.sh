#$ -cwd
#$ -l vf=20G
#$ -l h_vmem=20G
#$ -N ibdqc
#$ -m eas
#$ -M guobo.chen@uq.edu.au

#add batch to phenotype file
 Rscript FindBatch.R
#Merge
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBDBatch1_rev --bmerge IBDBatch2_rev.bed IBDBatch2_rev.bim IBDBatch2_rev.fam --make-bed --out IBDMerge 

#Common snp
 java -Xmx15G -jar /clusterdata/gc5k/bin/gear.jar comsnp --bfiles IBDBatch1_rev IBDBatch2_rev --out IBDBatch1_2
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBDMerge --extract IBDBatch1_2.comsnp --make-bed --out IBDMerge_tmp0
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBDMerge_tmp0 --exclude IBDMerge.nof --make-bed --out IBDMerge_tmp

#Update status
 Rscript UpdateIDs.R
# 1 only keep individuals has MLIDs
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBDMerge_tmp --remove BadSample.txt --make-bed --out IBDMerge_tmp2
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBDMerge_tmp2 --remove UnknownInd.txt --make-bed --out IBDMergeComSNPKeep

#remove outliers
 java -Xmx15G -jar /clusterdata/gc5k/bin/gear.jar comsnp --bfiles IBDMergeComSNPKeep /clusterdata/gc5k/ibscratch/gc5k/bin/HM3_founders_noATGC_autosome_naive_imputed --out ibd_hapmap
 java -Xmx15G -jar /clusterdata/gc5k/bin/gear.jar propc --bfile IBDMergeComSNPKeep --extract-score ibd_hapmap.comsnp --score /clusterdata/gc5k/ibscratch/gc5k/bin/HM3_SNP.blup20 --out IBDMergeComSNPKeep
 java -Xmx15G -jar /clusterdata/gc5k/bin/gear.jar propc --bfile /clusterdata/gc5k/ibscratch/gc5k/bin/HM3_founders_noATGC_autosome_naive_imputed --extract-score ibd_hapmap.comsnp --score /clusterdata/gc5k/ibscratch/gc5k/bin/HM3_SNP.blup20 --out HapMap

 Rscript PCA.R
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBDMergeComSNPKeep --remove Admixed.txt --make-bed --out IBD_UC_tmp

# standard QC
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp --missing --out IBD_UC_tmp
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp --freq --out IBD_UC_tmp
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp --hardy --out IBD_UC_tmp
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp --het --out IBD_UC_tmp

 Rscript statQC.R
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp --exclude RevSNP.txt --make-bed --out IBD_UC_tmp2
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp2 --exclude gwas_strange_loci.txt --make-bed --out IBD_UC_tmp3

#pseudo gwas for two batches
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp3 --pheno IBD_UC.phe --mpheno 10 --assoc --out Batch_GWAS
 Rscript BatchGWAS.R
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_tmp3 --exclude Batch_GWAS_loci.txt --make-bed --out IBD_UC_clean

# realcheck
 java -Xmx15G -jar /clusterdata/gc5k/bin/gear.jar --realcheck --bfile IBD_UC_clean --realcheck-marker-number 10000 --out IBD_UC_clean
 /ibscratch/wrayvisscher/jyang/bin/gcta64_test --bfile IBD_UC_clean --make-grm --out IBD_UC_clean --autosome
 /ibscratch/wrayvisscher/jyang/bin/gcta64_test --grm IBD_UC_clean --grm-cutoff 0.05 --make-grm --out IBD_UC_clean1 --autosome
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_clean --keep IBD_UC_clean1.grm.id --make-bed --out IBD_UC_clean1
 /ibscratch/wrayvisscher/jyang/bin/gcta64_test --bfile IBD_UC_clean1 --make-grm --out IBD_UC_clean1 --autosome

# GWAS
# /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_clean1 --fisher --out IBD_UC_clean1
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --noweb --bfile IBD_UC_clean1 --assoc --out IBD_UC_clean1
 Rscript StatPlot.R
