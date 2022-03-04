#SCRIPT TO CLEAN ABCD GENETIC DATA using PLINK QC package
#set required directories and dataset names

library(plinkQC)
cd ~
qcdir="./ABCD_cleaning/qcdir"
plink2="./ABCD_cleaning/reference/plink2"
path2plink="./ABCD_cleaning/qcdir"
name="ABCD_data_withsex"
refdir="./ABCD_cleaning/reference"




######################################################################################################################################################################################
#INDIVIDUAL-LEVEL QC
######################################################################################################################################################################################
##FOR INDIVIDUAL-LEVEL QC, USE PLINK 1.90b as plink2 has not implemented --het yet
#####################################################################################################################################################################################

fail_individuals <- perIndividualQC(indir=qcdir, qcdir=qcdir, 
                                    name,
				    do.run_check_relatedness=TRUE, imissTh=0.03,
                                    path2plink=path2plink, verbose=TRUE, interactive=FALSE)

######################################################################################################################################################################################
##SAVE DESC STATS OF NO OF PARTICIPANTS FAILING QC, OVERALL
#####################################################################################################################################################################################

ggplot2::ggsave("fail_individuals_ABCD.png")

overview_individuals <- overviewPerIndividualQC(fail_individuals,
                                                interactive=FALSE)
ggplot2::ggsave("overview_individuals_ABCD.png")

