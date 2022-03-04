#!/bin/bash
#PBS -N abcd_ancestry_preprocessing
#PBS -e ~/shell_logs/qc/abcd_ancestry_preprocessing
#PBS -o ~/shell_logs/qc/abcd_ancestry_preprocessing
#PBS -l walltime=30:00:00,nodes=1:ppn=4
#PBS -S /bin/bash


Roxanna Korologou-Linden 10/06/2020
Adapted code from PLINK QC by Hannah Meyer

###################################################################
##set dirs and filenames
###################################################################

refname="all_phase3"
high_ld="~/ABCD_cleaning/reference/ext_data/high-LD-regions-hg19-GRCh37.txt"
qcdir="${HOME}/ABCD_cleaning/qcdir"
refdir="${HOME}/ABCD_cleaning/reference"
name="ABCD_data_with_sex.clean"
MYEX="${HOME}/ABCD_cleaning/reference/plink2"


####################################################################
##Download and decompress reference dataset (1000 genomes)
####################################################################


cd $refdir
pgen=https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/op9osq6luy3pjg8/all_phase3.pvar.zst?dl=1
sample=https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1
wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst

$MYEX --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

wget $pvar
$MYEX 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst

wget $sample
$MYEX 'phase3_corrected.psam?dl=1' all_phase3.psam




#######################################################################################################################################################
### Filter reference and study data for non A-T or G-C SNPs
#######################################################################################################################################################

##Identify duplicate variants in reference dataset
awk 'BEGIN {OFS="\t"}  a[$2]++; (a[$2] == 2) {print $2}' \
    $refdir/$refname.bim  > \
    $refdir/$refname.duplicate_variants                                                                    ##resulting file has 2,057,345 duplicated variants                                                            

##Identify duplicate variants in ABCD dataset (resulting file 3459 duplicate variants)
awk 'BEGIN {OFS="\t"}  a[$2]++; (a[$2] == 2) {print $2}' \
    $qcdir/$name.bim  > \
    $qcdir/$name.duplicate_variants                                                                        ##resulting file has 2,715 duplicated variants

##Create genotype files with unique variants in reference dataset                                          
$MYEX --bfile  $refdir/$refname \
      --exclude $refdir/$refname.duplicate_variants \
      --make-bed \
      --out $refdir/$refname.no_duplicate_variants                                                         ##following removal of duplicated variants, there are 82,301,086 SNPs in reference            

$MYEX --bfile  $qcdir/$name \
      --exclude $qcdir/$name.duplicate_variants \
      --make-bed \
      --out $qcdir/$name.no_duplicate_variants                                                             ##following removal of duplicated variants, 375,441 SNPs remain in genotype file



##Create file with AC AND GT SNPS  (resulting file has 35218 AC GT snps)
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.no_duplicate_variants.bim  > \
    $qcdir/$name.ac_gt_snps                                                                                ##There are 24,466 AC and GT SNPs in the study dataset

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.no_duplicate_variants.bim  > \
    $refdir/$refname.ac_gt_snps                                                                            ##There are 12,361,963 AC and GT SNPs in the reference dataset

#####################################################################################################################################################################################
##Create genotype files with unique variants in reference dataset
###########################################################

$MYEX --bfile  $refdir/$refname.no_duplicate_variants \
      --exclude $refdir/$refname.ac_gt_snps \
      --make-bed \
      --out $qcdir/$refname.no_duplicate_variants_ac_gt_snps
mv  $qcdir/$refname.no_duplicate_variants_ac_gt_snps.log $qcdir/plink_log/$refname.no_duplicate_variants_ac_gt_snps.log   ##After removal of duplicates and AC and GT SNPs, there are 69,939,123 SNPs

#####################################################################################################################################################################################
##Create genotype files with unique variants in study dataset 
###########################################################

$MYEX --bfile $qcdir/$name.no_duplicate_variants --exclude $qcdir/$name.ac_gt_snps  --make-bed --out $qcdir/$name.no_duplicate_variants.ac_gt_snps #350,975 variants in file
mv  $qcdir/$name.no_duplicate_variants.ac_gt_snps.log $qcdir/plink_log/$name.no_duplicate_variants.ac_gt_snps.log         

#####################################################################################################################################################################################
### Prune study data
##We will do principal component analysis on genetic variants that are pruned for variants in linkage disequlibrium (LD) with an $r^2 >0.2$ in a 50kb window. The LD-pruned dataset is generated below, using plink --indep-pairwise to compute the LD-variants; additionally exclude range is used to remove genomic
ranges of known high-LD structure. This file was originally provided by Anderson (look in plinkQC documentation for further info)
#####################################################################################################################################################################################


plink --bfile  $qcdir/$name.no_duplicate_variants.ac_gt_snps \
      --exclude range  $high_ld \
      --indep-pairwise 50 5 0.2 \
      --out $qcdir/$name.no_duplicate_variants.ac_gt_snps
mv  $qcdir/$name.no_duplicate_variants.ac_gt_snps.log $qcdir/plink_log/$name.no_duplicate_variants.ac_gt_snps.log  ##152,643 pruned SNPs

#####################################################################################################################################################################################
##Extract full data on pruned variants for ABCD
#####################################################################################################################################################################################

plink --bfile  $qcdir/$name.no_duplicate_variants.ac_gt_snps \
      --extract $qcdir/$name.no_duplicate_variants.ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned
mv  $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned.log $qcdir/plink_log/$name.no_duplicate_variants.ac_gt_snps.pruned.log  


######################################################################################################################################################################################
### Filter reference data for the same SNP set as in study 
###We will use the list of pruned variants from the study sample to reduce the reference dataset to the size of the study samples:
#######################################################################################################################################################################################

$MYEX --bfile  $qcdir/$refname.no_duplicate_variants_ac_gt_snps  \
      --extract $qcdir/$name.no_duplicate_variants.ac_gt_snps.prune.in  \
      --make-bed \
      --out $qcdir/$refname.no_duplicate_variants.ac_gt_snps.pruned
mv  $qcdir/$refname.no_duplicate_variants.ac_gt_snps.pruned.log $qcdir/plink_log/$refname.no_duplicate_variants.ac_gt_snps.pruned.log   ##resulting file has 150,423 SNPs

########################################################################################################################################################################################
### Check and correct chromosome mismatch
##check that the variant IDs of the reference data have the same chromosome ID as the study data. For computing the genetic PC, the annotation is not important, however, merging the files via PLINK will only work for variants with perfectly
matching attributes. For simplicity, we update the pruned reference dataset. We'll ignore XY-encoded sex chromosomes (via `sed -n '/^[XY]/!p'`).
########################################################################################################################################################################################


awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned.bim $qcdir/$refname.no_duplicate_variants.ac_gt_snps.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.no_duplicate_variants.no_ac_gt_snps.toUpdateChr                                                  ##4,126 SNPs needing chromosome number updated

plink --bfile $qcdir/$refname.no_duplicate_variants.ac_gt_snps.pruned \
      --update-chr $qcdir/all_phase3.no_duplicate_variants.no_ac_gt_snps.toUpdateChr 1 2 \
      --make-bed \
      --out $qcdir/$refname.no_duplicate_variants.ac_gt_snps.updateChr
mv $qcdir/$refname.no_duplicate_variants.ac_gt_snps.updateChr.log $qcdir/plink_log/$refname.no_duplicate_variants.ac_gt_snps.updateChr.log   ##this returns the reference file with chromosome number for SNPs needing updating, updated.

##############################################################################################################################################################################################
### Position mismatch
##############################################################################################################################################################################################

#Similar to the chromosome matching, we use an awk-script to find variants with mis-matching chromosomal positions.
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned.bim $qcdir/$refname.no_duplicate_variants.ac_gt_snps.pruned.bim > \
    $qcdir/$refname.no_duplicate_variants.ac_gt_snps.toUpdatePos                                                                   ##364 SNPs have have a position mismatch between study and reference dataset.
 
###############################################################################################################################################################################################
### Possible allele flips
###############################################################################################################################################################################################
#Unlike chromosomal and base-pair annotation, mismatching allele-annotations will not only prevent the plink --merge, but also mean that it is likely that actually a different genotype was measured. Initially, we can use the following awk-script to check if non-matching allele codes are a simple case of allele
#flips. 
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned.bim $qcdir/$refname.no_duplicate_variants.ac_gt_snps.pruned.bim >  $qcdir/$refname.no_duplicate_variants.ac_gt_snps.toFlip    ##178 SNPs need flipping

###############################################################################################################################################################################################
### Update positions and flip alleles
###############################################################################################################################################################################################


plink --bfile $qcdir/$refname.no_duplicate_variants.ac_gt_snps.updateChr \
      --update-map $qcdir/$refname.no_duplicate_variants.ac_gt_snps.toUpdatePos 1 2 \
      --flip $qcdir/$refname.no_duplicate_variants.ac_gt_snps.toFlip \
      --make-bed \
      --out $qcdir/$refname.no_duplicate_variants.ac_gt_snps.flipped                                                               ##150,423 variants in file
mv $qcdir/$refname.no_duplicate_variants.ac_gt_snps.flipped.log $qcdir/plink_log/$refname.no_duplicate_variants.ac_gt_snps.flipped.log

#############################################################################################################################################################################################
### Remove mismatches
##############################################################################################################################################################################################

Any alleles that do not match after allele flipping, are identified and removed
from the reference dataset.
```{bash mismatch, eval=FALSE}
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned.bim $qcdir/$refname.no_duplicate_variants.ac_gt_snps.flipped.bim > \
    $qcdir/$refname.no_duplicate_variants.ac_gt_snps.mismatch                                                                     ##549 mismatches

plink --bfile $qcdir/$refname.no_duplicate_variants.ac_gt_snps.flipped \
      --exclude $qcdir/$refname.no_duplicate_variants.ac_gt_snps.mismatch \
      --make-bed \
      --out $qcdir/$refname.no_duplicate_variants.ac_gt_snps.clean                                                                ##149,874 SNPs in the final reference dataset to be merge following mismatch removal
mv $qcdir/$refname.no_duplicate_variants.ac_gt_snps.clean.log $qcdir/plink_log/$refname.no_duplicate_variants.ac_gt_snps.clean.log


##############################################################################################################################################################################################
################MY ADDED STEP where i excluded mismatched from study dataset
###############################################################################################################################################################################################

plink --bfile $qcdir/$name.no_duplicate_variants.ac_gt_snps.pruned \
      --exclude $qcdir/$refname.no_duplicate_variants.ac_gt_snps.mismatch \
      --make-bed \
      --out $qcdir/$name.no_duplicate_variants.ac_gt_snps.clean                                                                   ##152,094 SNPs after we remove the mismatched SNPs


plink --bfile $qcdir/$name.no_duplicate_variants.ac_gt_snps.clean  \
      --bmerge $qcdir/all_phase3.no_duplicate_variants.ac_gt_snps.clean.bed $qcdir/all_phase3.no_duplicate_variants.ac_gt_snps.clean.bim \
         $qcdir/all_phase3.no_duplicate_variants.ac_gt_snps.clean.fam  \
      --make-bed \
      --out $qcdir/$name.no_duplicate_variants.ac_gt_snps.merge.$refname                                                          ##152,094 SNPs
mv $qcdir/$name.no_duplicate_variants.ac_gt_snps.merge.$refname.log $qcdir/plink_log



#############################################################################################################################################################################################
## PCA on the merged data
##############################################################################################################################################################################################

#We can now run principal component analysis on the combined dataset using
#plink --pca which returns a .eigenvec file with the family and individual ID
#in columns 1 and 2, followed by the first 20 principal components. 
plink --bfile $qcdir/$name.no_duplicate_variants.ac_gt_snps.merge.$refname \
      --pca \
      --out $qcdir/$name.no_duplicate_variants.ac_gt_snps.merge.$refname
mv $qcdir/$name.no_duplicate_variants.ac_gt_snps.merge.$refname $qcdir/plink_log

