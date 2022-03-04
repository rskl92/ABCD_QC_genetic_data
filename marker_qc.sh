#!/bin/bash
#PBS -N abcd_qc
#PBS -e ~/shell_logs/qc/marker_qc_abcd
#PBS -o ~/shell_logs/qc/marker_qc_abcd
#PBS -l walltime=01:00:00,nodes=1:ppn=4
#PBS -S /bin/bash


######################################################################################################################################################################################
#MARKER-LEVEL QC
######################################################################################################################################################################################
##FOR MARKER-LEVEL QC, USE PLINK2 AS --hardy/--hwe skips chrY and chrM, but uses a specialized algorithm for chrX which takes sex into account.
#####################################################################################################################################################################################
cd ~

qcdir="./ABCD_cleaning/qcdir"
plink2="./plink2"
name="ABCD_data_with_sex"


$plink2 --bfile $qcdir/ABCD_data_with_sex --remove $qcdir/ABCD_data_with_sex.remove.IDs --not-chr 25-26 --maf 0.01 --hwe 5e-07  --geno 0.05 --make-bed --out $qcdir/ABCD_data_with_sex.clean



