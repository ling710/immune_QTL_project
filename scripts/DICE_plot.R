setwd("/Users/tq20202/Desktop/DICE")
res<-read.csv("topfinding.csv")
list<-unique(res[res$PPH4>=0.8,"gene"])
fname<-list.files("mrres")
plotres<-c()
setwd("/Users/tq20202/Desktop/DICE/mrres")
for (i in 1:length(fname)){
  temp<-read.table(file=fname[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp$gene<-sapply(strsplit(temp$id,split="_"),"[",1)
  temp<-temp[temp$gene %in% list,]
  if(nrow(temp)>0){temp$cancer<-strsplit(fname[i],split="_")[[1]][2]}
  plotres<-rbind(plotres,temp)
}
write.csv(plotres,"/Users/tq20202/Desktop/DICE/plot/data.csv")

library(ggplot2)
setwd("/Users/tq20202/Desktop/DICE/plot")
plotdata<-read.csv("data.csv")
g<-unique(plotdata$gene)
for (i in 1:length(g)){
  dataset<-plotdata[plotdata$gene==g[i],]
  p <- ggplot(dataset, aes(or, cancer,col=cell))+ggtitle(paste(g[i],"MR forest plot in DICE cohert",sep=" "))
  p<-p + geom_point(size=3.6) +geom_errorbarh(aes(xmax =or_uci95, xmin = or_lci95), height = 0.4) +scale_x_continuous(limits= c(0.5, 1.5), breaks= seq(0.5, 1.5, 0.5)) +geom_vline(aes(xintercept = 1)) +xlab('Odds Ratio ') + ylab(' Cancer type')
  png(filename = paste(g[i],"MR forest plot in DICE cohert.png",sep=" "),width = 2000,height = 1200,res = 300)
  print(p)
  dev.off()
}








