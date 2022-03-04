library(data.table)

setwd("~")
setwd("./ABCD_cleaning/qcdir")
eaf_abcd <-fread("allphase3_clean_freq.frq")
eaf_allphase3 <- fread("abcd_clean_freq.frq")

colnames(eaf_abcd)[5] <- "freq_abcd"
eaf_abcd$identifier <- "abcd"
eaf_allphase3$identifier <- "allphase3"
colnames(eaf_allphase3)[5] <- "freq_allphase3"
data_merge <- merge(eaf_abcd,eaf_allphase3,by="SNP")


tiff("./ABCD_cleaning/qcdir/scatterplot.tiff")
plot <- plot(x=data_merge$freq_allphase3, y=data_merge$freq_abcd,xlab="EAF_1000g",ylab="EAF_ABCD")
dev.off()