rm(list=ls())
gc()
setwd("/Users/tq20202/Desktop/DICE")

f<-list.files("raw data")
setwd("/Users/tq20202/Desktop/DICE/raw data")
data<-c()
for (i in 1:length(f)){
  temp<-read.table(file=f[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp$cell<-strsplit(f[i],split = "[.]")[[1]][1]
  data<-rbind(data,temp)
}
data$exposure<-paste(data$gene,data$cell,sep="_")
data$se <- abs(as.numeric(data$beta))/qnorm(1-as.numeric(data$p)/2)
data$se<-data$se+0.00001
write.table(data, file="/Users/tq20202/Desktop/DICE/rawdata.txt", row.names=F, col.names=T, sep="\t", quote=F)

library(TwoSampleMR)
data <- read.table(file="/Users/tq20202/Desktop/DICE/rawdata.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
a <- data[data$p<=5e-8,]
dat <- data.frame(SNP = a$SNP,
                  chr = a$chr,
                  pos = a$position,
                  beta = as.numeric(a$beta),
                  se = as.numeric(a$se),
                  effect_allele = a$ALT,
                  other_allele = a$REF,
                  pval = as.numeric(a$p),
                  Phenotype = a$exposure,
                  samplesize = 1544)
dat <- format_data(dat, type="exposure", phenotype_col = "Phenotype",chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize")
# LD clumping 
dat <-clump_data(data, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, clump_p2 = 1,pop = "EUR")
write.table(dat, file="/Users/tq20202/Desktop/DICE/exposure.txt", row.names=F, col.names=T, sep="\t", quote=F)

