#########beta compare##################################
data<-read.table(file="/Users/tq20202/Desktop/DICE/raw data/B_CELL_NAIVE.txt", header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
data<-data[data$p<1e-5,]
naiveB<-read.csv("/Users/tq20202/Desktop/oneK1K/data/naiveB.csv")
naiveB<-naiveB[naiveB$P_VALUE<1e-5,]
data$match<-paste(data$SNP,data$gene,sep="_")
naiveB$match<-paste(naiveB$RSID,naiveB$GENE,sep="_")
i<-intersect(data$match,naiveB$match)
dice<-data[data$match %in% i,]
dice<-dice[!duplicated(dice$match),]
onek1k<-naiveB[naiveB$match %in% i,]
onek1k$A1_FREQ<-1-onek1k$A2_FREQ_ONEK1K
onek1k$beta<-as.numeric(onek1k$SPEARMANS_RHO)/sqrt(2*as.numeric(onek1k$A1_FREQ)*(1-as.numeric(onek1k$A1_FREQ)))
onek1k<-onek1k[!duplicated(onek1k$match),]
dice<-dice[order(dice$match),]
onek1k<-onek1k[order(onek1k$match),]

library(ggplot2)
for (i in 1:nrow(dice)){
  if (dice[i,"ALT"]!=onek1k[i,"A1"]){onek1k[i,"beta"]<-(-1)*onek1k[i,"beta"]}
}
dd<-data.frame(beta_dice=dice$beta,beta_oneK1K=onek1k$beta,alt_dice=dice$ALT,alt_onek1k=onek1k$A1,ref_dice=dice$REF,ref_onek1k=onek1k$A2)
head(dd)
ggplot(data = dd, mapping = aes(x = beta_dice, y = beta_oneK1K)) + geom_point(size = 3, alpha=0.3)




