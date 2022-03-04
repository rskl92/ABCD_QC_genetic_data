#!/bin/bash
#PBS -N abcd_qc
#PBS -e ~/shell_logs/qc/ancestry
#PBS -o ~/shell_logs/qc/ancestry
#PBS -l walltime=03:00:00,nodes=1:ppn=4
#PBS -S /bin/bash


workDir="~/ABCD_cleaning/qcdir"
cd $workDir
module load  languages/R-3.6.2-gcc9.1.0
Rscript $workDir/ancestry_check.R
