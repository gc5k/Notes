#$ -cwd
#$ -l vf=30G
#$ -l h_vmem=30G
#$ -N buildsrc
#$ -m eas
#$ -M guobo.chen@uq.edu.au

####################remove atgc for iChip
 java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile IBD_UQ_001 --order-ind IBD_UQ_001_0.2.grm.id --remove-atgc --make-bed --out IBD_UQ_001_noatgc
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile IBD_UQ_001_noatgc --update-map ichip.hg19.update.id.txt --update-name --make-bed --out IBD_UQ_001_tmp --noweb
 awk '{print $2}' IBD_UQ_001_tmp.bim | sort | uniq -d > iChip_dup_snp.txt
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile IBD_UQ_001_tmp --exclude iChip_dup_snp.txt --make-bed --out IBD_UQ_001_pre --noweb
 rm IBD_UQ_001_tmp.*
 rm IBD_UQ_001_noatgc.*

####################remove atgc for gwas
 grep -v NA$ cd0.dat3 > cd1.keep
 java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile cd0 --order-ind cd1.keep --remove-atgc --make-bed --out cd1

 grep -v NA$ uc0.dat > uc1.keep
 java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile uc0 --order-ind uc1.keep --remove-atgc --make-bed --out uc1

####################find common snps between gwas
 awk '{print $2}' cd1.bim > gSNP.txt
 awk '{print $2}' uc1.bim >> gSNP.txt
 sort gSNP.txt | uniq -d > gwasCommonSNP.txt
 rm gSNP.txt

####################extract common snps
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile cd1 --extract gwasCommonSNP.txt --make-bed --out cd1ComSNP --noweb
 rm cd1.*
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile uc1 --extract gwasCommonSNP.txt --make-bed --out uc1ComSNP --noweb
 rm uc1.*

####################extract IGbackboneSNP.txt
 cat gwasCommonSNP.txt > igComSNP.txt
 awk '{print $2}' IBD_UQ_001_pre.bim >> igComSNP.txt
 sort igComSNP.txt | uniq -d > IGbackboneSNP.txt
 rm igComSNP.txt

#####################realcheck
 head -1000 IGbackboneSNP.txt > RealSNP1K.txt
 java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --realcheck --bfile IBD_UQ_001_pre --bfile2 cd1ComSNP --realcheck-snps RealSNP1K.txt --realcheck-threshold-lower 0.3 --out ichip2cd
 java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --realcheck --bfile IBD_UQ_001_pre --bfile2 uc1ComSNP --realcheck-snps RealSNP1K.txt --realcheck-threshold-lower 0.3 --out ichip2uc


#1 GWAS vs iChip
##generate CD iChip
awk 'NR>1' plink.pheno.txt | awk '{print $1, $2, $3}' > plink.cd.txt
Rscript pheMatch.R 0.9 ichip2cd.real plink.cd.txt cd0.dat3 cd_ichip_keep.phe cd_gwas_keep.phe
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile IBD_UQ_001_pre --keep cd_ichip_keep.phe --make-bed --out tIBD_cd --noweb
java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile tIBD_cd --order-ind cd_ichip_keep.phe --make-bed --out Ocd1iChip
rm tIBD_cd.*
cp cd_ichip_keep.phe Ocd1iChip.phe

##generate CD gChip
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile cd1ComSNP --keep cd_gwas_keep.phe --make-bed --out tgwas_cd --noweb
java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile tgwas_cd --order-ind cd_gwas_keep.phe --make-bed --out Ocd1gChip
rm tgwas_cd.*
cp cd_gwas_keep.phe Ocd1gChip.phe

##generate UC iChip
awk 'NR>1' plink.pheno.txt | awk '{print $1, $2, $4}' > plink.uc.txt
Rscript pheMatch.R 0.9 ichip2uc.real plink.uc.txt uc0.dat uc_ichip_keep.phe uc_gwas_keep.phe
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile IBD_UQ_001_pre --keep uc_ichip_keep.phe --make-bed --out tIBD_uc --noweb
java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile tIBD_uc --order-ind uc_ichip_keep.phe --make-bed --out Ouc1iChip
rm tIBD_uc.*
cp uc_ichip_keep.phe Ouc1iChip.phe

##genarate UC gchip
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile uc1ComSNP --keep uc_gwas_keep.phe --make-bed --out tgwas_uc --noweb
java -Xmx25G -jar /clusterdata/gc5k/bin/gear.jar --bfile tgwas_uc --order-ind uc_gwas_keep.phe --make-bed --out Ouc1gChip
rm tgwas_uc.*
cp uc_gwas_keep.phe Ouc1gChip.phe


#2 make gwas+ichip
##merge cd gwas and ichip
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile Ocd1gChip --exclude IGbackboneSNP.txt --make-bed --out tcd1ComSNPRevBackbone --noweb
cp Ocd1iChip.bed tOcd1iChip.bed
cp Ocd1iChip.bim tOcd1iChip.bim
cp tcd1ComSNPRevBackbone.fam tOcd1iChip.fam
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile tcd1ComSNPRevBackbone --bmerge tOcd1iChip.bed tOcd1iChip.bim tOcd1iChip.fam --make-bed --out Ocd1IGSNP --noweb
rm tOcd1iChip.*
rm tcd1ComSNPRevBackbone.*
cp cd_gwas_keep.phe Ocd1IGSNP.phe

##merge uc gwas and ichip
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile Ouc1gChip --exclude IGbackboneSNP.txt --make-bed --out tuc1ComSNPRevBackbone --noweb
cp Ouc1iChip.bed tOuc1iChip.bed
cp Ouc1iChip.bim tOuc1iChip.bim
cp tuc1ComSNPRevBackbone.fam tOuc1iChip.fam
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile tuc1ComSNPRevBackbone --bmerge tOuc1iChip.bed tOuc1iChip.bim tOuc1iChip.fam --make-bed --out Ouc1IGSNP --noweb
rm tOuc1iChip.*
rm tuc1ComSNPRevBackbone.*
cp uc_gwas_keep.phe Ouc1IGSNP.phe


#3 make backbone
##make cd backbone
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile Ocd1gChip --extract IGbackboneSNP.txt --make-bed --out OBKcd1gChip --noweb
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile Ocd1iChip --extract IGbackboneSNP.txt --make-bed --out OBKcd1iChip --noweb
cp cd_gwas_keep.phe OBKcd1gChip.phe
cp cd_ichip_keep.phe OBKcd1iChip.phe

##make uc backbone
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile Ouc1gChip --extract IGbackboneSNP.txt --make-bed --out OBKuc1gChip --noweb
/clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile Ouc1iChip --extract IGbackboneSNP.txt --make-bed --out OBKuc1iChip --noweb
cp uc_gwas_keep.phe OBKuc1gChip.phe
cp uc_ichip_keep.phe OBKuc1iChip.phe


#4 make Large training set
###################extract iChip CD Large
awk 'NR>1' plink.pheno.txt | awk '{print $1, $2, $3}' | grep -v 0$ > plink.cd.txt
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile IBD_UQ_001_pre --keep plink.cd.txt --make-bed --out iChipCD_tmp0 --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile iChipCD_tmp0 --remove cd_ichip_keep.phe.dis --make-bed --out iChipCD_tmp --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile iChipCD_tmp --remove ANZAC.txt --make-bed --out LcdiChip --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile LcdiChip --extract IGbackboneSNP.txt --make-bed --out LBKcdiChip --noweb

#make iChip ANZAC
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile iChipCD_tmp --keep ANZAC.txt --make-bed --out AZcdiChip --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile AZcdiChip --extract IGbackboneSNP.txt --make-bed --out AZBKcdiChip --noweb

###################extract iChip UC Large
awk 'NR>1' plink.pheno.txt | awk '{print $1, $2, $4}' | grep -v 0$ > plink.uc.txt
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile IBD_UQ_001_pre --keep plink.uc.txt --make-bed --out iChipUC_tmp0 --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile iChipUC_tmp0 --remove uc_ichip_keep.phe.dis --make-bed --out iChipUC_tmp --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile iChipUC_tmp --remove ANZAC.txt --make-bed --out LuciChip --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile LuciChip --extract IGbackboneSNP.txt --make-bed --out LBKuciChip --noweb

#make iChip ANZAC
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile iChipUC_tmp --keep ANZAC.txt --make-bed --out AZuciChip --noweb
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile AZuciChip --extract IGbackboneSNP.txt --make-bed --out AZBKuciChip --noweb

###################extract gChip CD Large
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile cd1ComSNP --remove cd_gwas_keep.phe.dis --make-bed --out LcdgChip --noweb

###################extract backbone CD Large 
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile LcdgChip --extract IGbackboneSNP.txt --make-bed --out LBKcdgChip --noweb

###################extract gChip UC Large
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile uc1ComSNP --remove uc_gwas_keep.phe.dis --make-bed --out LucgChip --noweb

###################extract backbone UC Large
 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile LucgChip --extract IGbackboneSNP.txt --make-bed --out LBKucgChip --noweb


#5 generate projected pc from popref
 awk 'NR>1' /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup | awk '{print $1}' > score.snp
 awk '{print $2}' Ocd1gChip.bim >> score.snp
 sort score.snp | uniq -d > score.keep 

 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile Ocd1gChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out Ocd1gChip
 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile Ouc1gChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out Ouc1gChip

 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile LcdgChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out LcdgChip
 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile LucgChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out LucgChip

 /clusterdata/gc5k/bin/plink-1.07-x86_64/plink --bfile /clusterdata/gc5k/ibscratch/gc5k/bin/popres_maf0.01_info0.8 --extract score.keep --make-bed --out pop_600k --noweb
 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile pop_600k --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out popres

 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile LBKcdiChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out LBKcdiChip
 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile LcdiChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out LcdiChip

 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile LBKuciChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out LBKuciChip
 java -Xmx35G -jar /clusterdata/gc5k/bin/gear.jar profile --bfile LuciChip --extract-score score.keep --score /clusterdata/gc5k/ibscratch/gc5k/bin/popres_HM3.blup --no-weight --out LuciChip

#6 CVsplit
Rscript CVSplit.R OBKcd1iChip.fam OBKcd1gChip.fam Ocd1iChip.fam Ocd1gChip.fam
Rscript CVSplit.R OBKuc1iChip.fam OBKuc1gChip.fam Ouc1iChip.fam Ouc1gChip.fam
Rscript 5FoldCVSplit.R LBKcdiChip.fam LBKcdgChip.fam LcdiChip.fam LcdgChip.fam
Rscript 5FoldCVSplit.R LBKuciChip.fam LBKucgChip.fam LuciChip.fam LucgChip.fam

#7 make adjustment
Rscript CovAdj.R OBKcd1iChip.trt OBKcd1gChip.trt Ocd1iChip.trt Ocd1gChip.trt Ocd1gChip.profile
Rscript CovAdj.R OBKuc1iChip.trt OBKuc1gChip.trt Ouc1iChip.trt Ouc1gChip.trt Ouc1gChip.profile

Rscript LCovAdj.R LBKcdgChip.trt LcdgChip.profile
Rscript LCovAdj.R LcdgChip.trt LcdgChip.profile
Rscript LCovAdj.R LBKucgChip.trt LucgChip.profile
Rscript LCovAdj.R LucgChip.trt LucgChip.profile

Rscript LCovAdj.R LBKcdiChip.trt LBKcdiChip.profile
Rscript LCovAdj.R LcdiChip.trt LcdiChip.profile
Rscript LCovAdj.R LBKuciChip.trt LBKuciChip.profile
Rscript LCovAdj.R LuciChip.trt LuciChip.profile

