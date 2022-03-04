library(plinkQC)
library(data.table)
##############################################################################################################################################################################################
##set dirs and filenames
##############################################################################################################################################################################################

refname="all_phase3"
qcdir="~/ABCD_cleaning/qcdir"
refdir="~/ABCD_cleaning/reference"
name="ABCD_data_with_sex"
MYEX="~/ABCD_cleaning/reference/plink2"
prefixMergedDataset="ABCD_data_with_sex.clean.no_duplicate_variants.ac_gt_snps.merge.all_phase3"
cleaneddataset<- fread("~/ABCD_cleaning/qcdir/ABCD_data_with_sex.clean.fam")

exclude_ancestry <- evaluate_check_ancestry(indir=qcdir, qcdir=qcdir, name=name,
refSamplesFile=paste(refdir, "/1000g_ID2POP.txt",
sep=""),
refColorsFile=paste(refdir, "/1000g_PopColors.txt",
sep=""),
prefixMergedDataset=prefixMergedDataset,
interactive=FALSE)


##SAVE PLOT OF PCA ANALYSIS
ggplot2::ggsave("p_ancestry.png")

##WRITE LIST OF IDS of EUROPEAN ANCESTRY 
europeans <- cleaneddataset[!cleaneddataset$V2 %in% exclude_ancestry$fail_ancestry,]
europeans <- europeans[,1:2]
colnames(europeans) <- c("FID","IID")

write.table(exclude_ancestry$fail_ancestry, "ABCD.non_european.txt",quote=FALSE,row.names=FALSE)
write.table(europeans, "ABCD.european.txt",quote=FALSE,row.names=FALSE)


