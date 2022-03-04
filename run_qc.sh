#!/bin/bash
#PBS -N abcd_qc
#PBS -e ~/shell_logs/qc/abcd_qc
#PBS -o ~/shell_logs/qc/abcd_qc
#PBS -l walltime=220:00:00,nodes=1:ppn=4
#PBS -S /bin/bash

cd "~/ABCD_cleaning/qcdir/"


##update sex by extracting it from the phenotypic dataset and remove plate 461
plink --bfile ABCD_release_2.0.1_r1 --remove plate_461.txt --update-sex pedsex_abcd_cohort.txt --make-bed --out ABCD_data_with_sex


workDir="~/ABCD_cleaning/qcdir"
cd $workDir
module load  languages/R-3.6.2-gcc9.1.0
Rscript $workDir/perindividual_qc.R
