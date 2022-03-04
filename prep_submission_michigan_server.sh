#!/bin/bash
#PBS -N abcd_qc
#PBS -e ~/shell_logs/qc/imputation_files
#PBS -o ~/shell_logs/qc/imputation_files
#PBS -l ~/walltime=02:00:00,nodes=1:ppn=4
#PBS -S /bin/bash

cd ~
qcdir="./ABCD_cleaning/qcdir"
plink2="./plink2"

module add  apps/vcftools-0.1.17.0
module add apps/tabix-0.2.6

for chr in {1..23}; do $plink2 --bfile $qcdir/ABCD_data_with_sex.clean --chr $chr  --recode vcf -snps-only just-acgt --make-bed --out $qcdir/ABCD_data_chr${chr};done 
for chr in {1..23}; do vcf-sort $qcdir/ABCD_data_chr$chr.vcf | bgzip -c > $qcdir/ABCD_data_chr$chr.vcf.gz;done 