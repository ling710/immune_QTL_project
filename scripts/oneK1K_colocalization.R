##################coloc###################################
library(coloc)
##coloc function##
coloc.analysis <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,s){
  dataset1 <- list(beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1)
  #type cc (case control) for binary study
  dataset2 <- list(beta=beta2, varbeta=se2^2,MAF=MAF2, type="cc",s=s,N=N2)
  #setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12)
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  return(df)
}

setwd("/Users/tq20202/Desktop/oneK1K")
fdrres<-read.csv("fdrres.csv")
setwd("/Users/tq20202/Desktop/oneK1K/data")
bres<-c()
for (i in 20:nrow(fdrres)){
  snp<-fdrres[i,"SNP"]
  chr<-fdrres[i,"chr"]
  pos<-fdrres[i,"pos"]
  celltype=fdrres[i,"cell"]
  cancertype=fdrres[i,"outcome"]
  file<-read.csv(fdrres[i,"cell_file"])
  data<-file[file$CHR==chr & file$POS>pos-500000 & file$POS<pos+500000,]
  data$beta<-data$SPEARMANS_RHO/sqrt(2*data$A2_FREQ_ONEK1K*(1-data$A2_FREQ_ONEK1K))
  data$z<-qnorm(data$P_VALUE,lower.tail = F)
  data$se <- abs(data$beta/data$z)
  data$se<-data$se+0.00001
  data$n=982
  outid<-fdrres[i,"outcome_id"]
  for(k in 1:nrow(data)){
    if (data[k,"A2_FREQ_ONEK1K"]>0.5){oa<-data[k,"A2"];ea<-data[k,"A1"];data[k,"A1"]<-oa;data[k,"A2"]<-ea;data[k,"A2_FREQ_ONEK1K"]=1-data[k,"A2_FREQ_ONEK1K"]}
  }
  exp=data.frame(SNP=data$RSID,effect_allele=data$A2,eaf=data$A2_FREQ_ONEK1K,other_allele=data$A1,beta=data$beta,se=data$se,p=data$P_VALUE,n=data$n)
  out <- extract_outcome_data(snps=exp$SNP, outcomes= outid)
  out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele.outcome,eaf=out$eaf.outcome,other_allele=out$other_allele.outcome,beta=out$beta.outcome,se=out$se.outcome,p=out$pval.outcome,n=outcome_samplesize[which(outcome_id==outid)],case=ncase[which(outcome_id==outid)])
  inter=intersect(exp$SNP,out$SNP)
  exp<-exp[exp$SNP %in% inter,]
  exp<-exp[order(exp$SNP),]
  exp <- exp[!duplicated(exp$SNP),]
  out<-out[out$SNP %in% inter,]
  out<-out[order(out$SNP),]
  out <- out[!duplicated(out$SNP),]
  pname1<-paste("/Users/tq20202/Desktop/oneK1K/colocalization/pwcoco/",snp,"_",celltype,"_",cancertype,"_exp.txt",sep="")
  pname2<-paste("/Users/tq20202/Desktop/oneK1K/colocalization/pwcoco/",snp,"_",celltype,"_",cancertype,"_out.txt",sep="")
  write.table(exp, file=pname1, row.names=F, col.names=T, sep="\t", quote=F)
  write.table(out, file=pname2, row.names=F, col.names=T, sep="\t", quote=F)
  df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$eaf),MAF2=as.numeric(out$eaf),N1=exp$n,N2=out$n,s=out$case/out$n)
  df <- na.omit(df)
  #coloc analysis
  result <- coloc.analysis(df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, df$N1, df$N2, df$s) 
  result <- data.frame(result,celltype=celltype,cancertype=cancertype,snp=snp)
  bres <- rbind(bres,result)
}
write.table(bres, file="/Users/tq20202/Desktop/oneK1K/colocalization/coloc_res.txt", row.names=F, col.names=T, sep="\t", quote=F)

#####LDCHECK#######################################
library(TwoSampleMR)
library(MendelianRandomization)
library(ieugwasr)
setwd("/Users/tq20202/Desktop/oneK1K")
fdrres<-read.csv("fdrres.csv")
setwd("/Users/tq20202/Desktop/oneK1K/data")
ldres<-c()
for (i in 1:nrow(fdrres)){
  rsid=fdrres[i,"SNP"]
  chr=fdrres[i,"chr"]
  pos=fdrres[i,"pos"]
  filename=fdrres[i,"cell_file"]
  outid=fdrres[i,"outcome_id"]
  file<-read.csv(filename)
  file<-file[file$CHR==chr & file$POS>pos-500000 & file$POS<pos+500000,]
  assoc<-extract_outcome_data(snps=file$RSID, outcomes= outid)
  if (nrow(assoc)!=0){
    assoc <- assoc[order(assoc$pval.outcome),]
    data <- assoc[assoc$pval.outcome<1E-3,]
    if(nrow(data)>=500){data <- data[1:499,] }
    print(nrow(data))
    
    if (nrow(data)!=0){
      data<-na.omit(data[,1:9])
      snp <- append(data$SNP, rsid)
      a <- NULL
      attempts <- 0
      while(attempts<=10){    
        a <- ld_matrix(variants=snp,pop = "EUR")
        if(is.null(a)){attempts<-attempts+1}else{break}
      }
      if(is.null(nrow(a))==TRUE){c<-cbind(snp=rsid, ld_snp=rsid, ld_r2=1)} 
      else {col.index <- which(grepl(rsid,colnames(a)))
      if (length(col.index)>0){
        b <- (a[,col.index])^2
        b <- b[order(b)]
        b <- b[(length(b)-1)] 
        c <- cbind(snp=rsid, ld_snp=names(b), ld_r2=as.numeric(b))
      } else {c <- cbind(snp=rsid, ld_snp="NA", ld_r2="NA")}
      }
    } else {
      c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  } else {c <- c(snp=rsid, ld_snp="NA",ld_r2="NA")}
  
  ldres <- rbind(ldres,c)
}
write.csv(ldres,"/Users/tq20202/Desktop/oneK1K/colocalization/ldres.csv")
  
  
  
  
  
  









