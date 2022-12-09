rm(list=ls())
gc()
setwd("/Users/tq20202/Desktop/oneK1K")

d1 <- read.delim("~/Desktop/oneK1K/data_5000000.tsv",header = T)
naiveB<-d1[d1$CELL_TYPE=="Naïve/Immature B Cell",]
naiveB<-naiveB[,1:20]
write.csv(naiveB,"naiveB.csv")
memoryB1<-d1[d1$CELL_TYPE=="Memory B Cell",]
memoryB1<-memoryB1[,1:20]
d2 <- read.delim("~/Desktop/oneK1K/data_10000000.tsv",header = F)
memoryB2<-d2[d2$V2=="Memory B Cell",]
colnames(memoryB2)<-colnames(memoryB1)
memoryB<-rbind(memoryB1,memoryB2)
write.csv(memoryB,"memoryB.csv")
CD4effect1<-d2[d2$V2=="CD4 Effector memory/TEMRA",]
d3 <- read.delim("~/Desktop/oneK1K/data_15000000.tsv",header = F)
CD4effect2<-d3[d3$V2=="CD4 Effector memory/TEMRA",]
CD4effect<-rbind(CD4effect1,CD4effect2)
colnames(CD4effect)<-colnames(naiveB)
write.csv(CD4effect,"CD4effect.csv")
CD4naive1<-d3[d3$V2=="CD4 Naive/Central memory T cell",]
d4 <- read.delim("~/Desktop/oneK1K/data_20000000.tsv",header = F)
CD4naive2<-d4[d4$V2=="CD4 Naive/Central memory T cell",]
CD4naive<-rbind(CD4naive1,CD4naive2)
colnames(CD4naive)<-colnames(naiveB)
write.csv(CD4naive,"CD4naive.csv")
CD4sox41<-d4[d4$V2=="CD4 SOX4 T cell",]
d5 <- read.delim("~/Desktop/oneK1K/data_25000000.tsv",header = F)
CD4sox42<-d5[d5$V2=="CD4 SOX4 T cell",]
CD4sox4<-rbind(CD4sox41,CD4sox42)
colnames(CD4sox4)<-colnames(naiveB)
write.csv(CD4sox4,"CD4sox4.csv")
CD8effect1<-d5[d5$V2=="CD8 Effector memory",]
d6 <- read.delim("~/Desktop/oneK1K/data_30000000.tsv",header = F)
CD8effect2<-d6[d6$V2=="CD8 Effector memory",]
CD8effect<-rbind(CD8effect1,CD8effect2)
colnames(CD8effect)<-colnames(naiveB)
write.csv(CD8effect,"CD8effect.csv")
CD8naive1<-d6[d6$V2=="CD8 Naive/Central memory T cell",]
d7 <- read.delim("~/Desktop/oneK1K/data_35000000.tsv",header = F)
CD8naive2<-d7[d7$V2=="CD8 Naive/Central memory T cell",]
CD8naive<-rbind(CD8naive1,CD8naive2)
colnames(CD8naive)<-colnames(naiveB)
write.csv(CD8naive,"CD8naive.csv")
CD8s100b1<-d7[d7$V2=="CD8 S100B T cell",]
d8 <- read.delim("~/Desktop/oneK1K/data_40000000.tsv",header = F)
CD8s100b2<-d8[d8$V2=="CD8 S100B T cell",]
CD8s100b<-rbind(CD8s100b1,CD8s100b2)
colnames(CD8s100b)<-colnames(naiveB)
write.csv(CD8s100b,"CD8s100b.csv")
dend1<-d8[d8$V2=="Dendritic Cell",]
d9 <- read.delim("~/Desktop/oneK1K/data_45000000.tsv",header = F)
dend2<-d9[d9$V2=="Dendritic Cell",]
dend<-rbind(dend1,dend2)
colnames(dend)<-colnames(naiveB)
write.csv(dend,"dend.csv")
mono1<-d9[d9$V2=="Classic Monocyte",]
d10 <- read.delim("~/Desktop/oneK1K/data_50000000.tsv",header = F)
mono2<-d10[d10$V2=="Classic Monocyte",]
mono<-rbind(mono1,mono2)
colnames(mono)<-colnames(naiveB)
write.csv(mono,"mono.csv")
nonclamono1<-d10[d10$V2=="Non-classic Monocyte",]
d11 <- read.delim("~/Desktop/oneK1K/data_55000000.tsv",header = F)
nonclamono2<-d11[d11$V2=="Non-classic Monocyte",]
nonclamono<-rbind(nonclamono1,nonclamono2)
colnames(nonclamono)<-colnames(naiveB)
write.csv(nonclamono,"nonclamono.csv")
killer1<-d11[d11$V2=="Natural Killer Cell",]
d12 <- read.delim("~/Desktop/oneK1K/data_60000000.tsv",header = F)
killer2<-d12[d12$V2=="Natural Killer Cell",]
killer<-rbind(killer1,killer2)
colnames(killer)<-colnames(naiveB)
write.csv(killer,"killer.csv")
killerrecru1<-d12[d12$V2=="Natural Killer Recruiting Cell",]
d13 <- read.delim("~/Desktop/oneK1K/data_64374850.tsv",header = F)
killerrecru2<-d13[d13$V2=="Natural Killer Recruiting Cell",]
killerrecru<-rbind(killerrecru1,killerrecru2)
colnames(killerrecru)<-colnames(naiveB)
write.csv(killerrecru,"killerrecru.csv")
plasma<-d13[d13$V2=="Plasma Cell",]
colnames(plasma)<-colnames(naiveB)
write.csv(plasma,"plasma.csv")











