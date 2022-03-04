# Quality Control of ABCD Genetic Data Release v. 2019.05
Date of release: Tues May 16 12:00:00 PST 2019
Web-page: abcdstudy.org
doi: http://dx.doi.org/10.15154/1503209

## Details on genotype platform

1. genotyping platform
    Affymetrix NIDA SmokeScreen Array
    The sample preparation and genotyping are performed by Rutgers RUCDR. This includes:
      - extraction kit: Chemagen bead based/Chemagic STAR DNA Saliva4k Kit (CMG-1755-A)
      - processing: DNA fragmentation, labeling, ligation & hybridization
      - equipment: Affimetrix GeneTitan Instrument.
    Additional information can be found in the NIMH experiment description #1194.
    The smokescreen array contains 733,293 SNPs (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4769529/)

The initial dataset provided by the ABCD (ABCD_release_2.0.1_r1) included 517,724 genetic variants chr1-23,25-26. 

## Introduction

This document describes the quality control procedure applied to the Adolescent Brain Cognitive Development study (ABCD) genetic data. A cohort description can be found here 
(https://www.sciencedirect.com/science/article/pii/S1878929317301822). The sample quality control of the genetic data was undertaken usin the PLINK QC package by Hannah Meyer 
(https://cran.r-project.org/web/packages/plinkQC/vignettes/plinkQC.pdf).





### 2.3 Ancestry restrictions 

The 1000 genomes reference panels was downloaded from (https://www.cog-genomics.org/plink/2.0/resources) and processed using the pipeline in this vignette
(https://meyer-lab-cshl.github.io/plinkQC/articles/Genomes1000.html). To estimate joint principal components of the reference and study population, we’ll need to combine the 
two datasets. The plink –merge function allows this merge, but requires the variants in the datasets to be matching by chromosome, position and alleles. 
We merge the cleaned dataset (after removal of individuals and variants failing QC checks) with the reference dataset.

```bash
cd ~/ABCD_cleaning/scripts
qsub ancestry_pre_processing.sh
```

We identified individuals of European ancestry by combining the genotypes of ABCD with genotypes of 1000 genomes phase 3, consisting of individuals from known ethnicities. 
Principal compoent analysis using this genotype panel can be used to identify population structure down to the level of 1000 genomes (i.e. large scale continental ancestry).
To identify these individuals, we used check_ancestry implemented in PLINK QC. It uses principal components 1 and 2 to find the centre of the European reference samples.
We performed PCA analysis on the pruned ABCD dataset (No of variants=152,094)
All study samples whose euclidean distance from the centre falls outside a specified radius are considered non-European. 

We run a bash script to submit a R script (ancestry_check.R) implementing plink QC, to estimate ancestries.

```bash
cd ~/ABCD_cleaning/scripts
qsub run_check_ancestry.sh
```


• ABCD.european.txt – 5,300 individuals to include                                                                                                                                      
• ABCD.non_european.txt - 4,607 individuals to exclude




## Individual/sample level QC

The plate 461 was found to be especially problematic by the ABCD study team.They recommended researchers to avoid using genotype data from plate 461, or at least including 
plate number as the covariates. Hence, we removed individuals from this plate entirely (N=56). Furthermore, the genetic data file did not include the sex of the subjects.
The sex was extracted from phenotypic datasets and the .fam file was updated, accordingly.                                                                                         
The per-individual quality control with perIndividualQC wraps around these functions in this instance: (i) check_sex: for the identification of individuals with discordant sex 
information, (ii) check_heterozygosity_and_missingness: for the identification of individuals with outlying missing genotype and/or heterozygosity rates, 
and (iii)check_relatedness: for the identification of related individuals. 
The removal of related individuals depends on future use of the genetic data. The check_relatedness is used and plink --genome calculates identity by state (IBS) for each pair of individuals based on the average proportion of alleles 
shared at genotyped SNPs. The degree of recent shared ancestry (IBD) is estimated from genome-wide IBS. The proportion of IBD between two individuals is returned as PI_HAT.  
The function finds pairs of samples whose proportion of IBD is greater than the specified high IBDTh (0.1875). For pairs of individuals that do not have additional relatives 
in the dataset, the individual with the highest genotype missingness rate is selected and returned as the individual failing the relatedness check. 
We use plink 1.90b as --het has not been implemented in plink2 yet.

Filters: Females (X chrom homozygosity<0.2) Males (X chrom homozygosity>0.8); Missing genotype rate: 0.03; Heterozygosity rates: +/- 3 SD from mean heterozygosity rate



We run this bash script to submit perindividual_qc.R to the queue,to perform individual QC.

```bash
cd ~/ABCD_cleaning/scripts
qsub run_qc.sh
```


The resulting files:

• ABCD_data_with_sex.fail-sexcheck.txt – 210 individuals to exclude – this list of individuals has been derived by comparing genetic sex (inferred data using the homozygosity
rates across all X-chromosomal genetic markers with the expected rates (typically <0.2 for females and >0.8 for males). Hence, individals in this file have an assigned sex 
(PEDSEX in the .fam file) that is different to the sex inferred from the homozygosity rates (SNPSEX).                                                        
•ABCD_data_with_sex.fail-het.txt – 246 individuals to exclude - individuals that are outliers in heterozygosity and missing rates.                                               
•ABCD_data_with_sex.fail-imiss.txt - 283 individuals to exclude - individuals with outlying missing genotype rates.                                                              
•ABCD_data_with_sex.remove.IDs - 664 non-overlapping individuals to exclude - individuals in the above three files with duplicates removed. 
•ABCD_data_with_sex.fail-IBD.IDs (N=1,756)
•ABCD_data_with_sex.fail-ancestry (N=4,607)




## Marker-level QC 

Here, we will exclude the individuals who failed the per-individual quality control in the previous step and filter the variants. This step will not use PLINKQC
(which currently uses PLINK 1.90b),as PLINK2 implements a specialised algorithm (Grafelman and Weir's  extended chrX Hardy-Weinberg exact test), which takes male frequencies
into account (https://pubmed.ncbi.nlm.nih.gov/27071844/). We run the script marker_qc.sh containing the code below. In this step, we also remove the pseudo-autosomal
chromosal (Chr 25) and well as the mitochondrial chromosome (Chr 26). This produces a cleaned dataset (we do not exclude individuals who failed ancestry or IBD in this step).
The resulting dataset contains 9,907 individuals and 377,164 SNPs.

Filters: HWE: 5e-07 ; Call rate: 0.05 ; MAF: 0.01

```bash 
cd ~

qcdir="./ABCD_cleaning/qcdir"
plink2="./plink2"
name="ABCD_data_with_sex"


$plink2 --bfile $qcdir/ABCD_data_with_sex --remove $qcdir/ABCD_data_with_sex.remove.IDs --not-chr 25-26 --maf 0.01 --hwe 5e-07 \
--geno 0.05 --make-bed --out $qcdir/ABCD_data_with_sex.clean

```


We convert the binary PLINK files to a suitable format for imputation (.vcf) and remove non-ACGT SNPs using the prep_submission_michigan_server.sh - with the line of code below.

```bash

for chr in {1..23}; do $plink2 --bfile $qcdir/ABCD_data_with_sex.clean --chr $chr  --recode vcf -snps-only just-acgt --make-bed --out $qcdir/ABCD_data_chr${chr};done 
for chr in {1..23}; do vcf-sort $qcdir/ABCD_data_chr$chr.vcf | bgzip -c > $qcdir/ABCD_data_chr$chr.vcf.gz;done 
```







