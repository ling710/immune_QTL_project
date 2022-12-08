rm(list=ls())
gc()
setwd("/Users/tq20202/Desktop/main")
f<-list.files("exposure")
setwd("/Users/tq20202/Desktop/main/exposure")
data<-c()
for (i in 1:length(f)){
  temp<-read.table(file=f[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  temp$cell<-strsplit(f[i],split="_ld")[[1]][1]
  data<-rbind(data,temp)
}
data <- data.frame(SNP = data$SNP,
                  chr = data$chr.exposure,
                  pos = data$pos.exposure,
                  beta = as.numeric(data$beta.exposure),
                  se = as.numeric(data$se.exposure),
                  effect_allele = data$effect_allele.exposure,
                  other_allele = data$other_allele.exposure,
                  pval = as.numeric(data$pval.exposure),
                  Phenotype = paste(data$exposure,data$cell,sep="_"),
                  samplesize = 1544)
data <- format_data(data, type="exposure", phenotype_col = "Phenotype",chr_col = "chr",pos_col = "pos",samplesize_col = "samplesize")
write.csv(data,"/Users/tq20202/Desktop/main/exposure.csv")

###########MR analysis###################################
rm(list=ls())
gc()
setwd("/Users/tq20202/Desktop/main")

library(TwoSampleMR)
library(MendelianRandomization)

outcome_id <- c("ukb-b-2160","ieu-b-4953","ieu-a-965","ieu-a-967","ieu-b-4963","bbj-a-140",
                "ukb-b-1251","ukb-d-C43","finn-b-C3_STOMACH","ukb-b-16890","ukb-b-8777",
                "ukb-b-20145","bbj-a-117","finn-b-C3_GBM","ieu-b-4912","ukb-b-8193")
outcome_name <- c("prostate_cancer","liver_cancer","lung_cancer","Squamous cell lung_cancer",
                  "Ovarian_cancer","Pancreatic_cancer","rectum_cancer","skin_cancer",
                  "stomach_cancer","breast_cancer","cervical_cancer","colon_cancer",
                  "Esophageal_cancer","Brain_cancer","head and neck_cancer","bladder_cancer")
outcome_samplesize <- c(463010,372184,18336,18313,199741,196187,463010,361194,218792,462933,
                        462933,462933,197045,218701,373122,462933)
ncase <- c(3436,168,3442,3275,1218,442,1470,1672,633,10303,1889,1494,1300,91,1106,1101)
ncontrol <- c(459574,372016,14894,15038,198523,195745,461540,359522,218159,452630,461044,
              461439,195745,218792,372016,461832)

exposure<-read.csv("exposure.csv")
  for (j in 1:length(outcome_id)){
    exposure_dat<-exposure[exposure$pval.exposure<=5e-8,]
    ###F-statistics for each SNPs
    exposure_dat <- cbind(exposure_dat,fstatistics=1)
    for (s in 1:nrow(exposure_dat)){
      z <- exposure_dat[s,"beta.exposure"]/exposure_dat[s,"se.exposure"]
      pve <- z^2/(1544+z^2)
      exposure_dat[s,"fstatistics"] <- (1544-2)*pve/(1-pve)
    }
    print(min(exposure_dat$fstatistics))
    print(max(exposure_dat$fstatistics))
    exposure_dat <- exposure_dat[exposure_dat$fstatistics>10,]
    
    ##############2SMR#######################################################
    #########################################################################
    outcome <- extract_outcome_data(snps=exposure_dat$SNP, outcomes= outcome_id[j])
    outcome<-outcome[!duplicated(outcome$SNP),]
    inter<-intersect(exposure_dat$SNP,outcome$SNP)
    exposure_dat <- exposure_dat[exposure_dat$SNP %in% inter,]
    outcome_dat <- outcome[outcome$SNP %in% inter,]
    
    harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
    #harmdat$mr_keep<-TRUE  #consider ambiguous
    
    test<-harmdat
    test$new<- paste(test$SNP,test$exposure,sep="_")
    test<-test[,-39]
    test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                        beta_col = "beta.exposure", se_col = "se.exposure",
                        effect_allele_col = "effect_allele.exposure",
                        other_allele_col = "other_allele.exposure",
                        phenotype_col = "new",samplesize_col = "samplesize.exposure")
    testout  <- outcome_dat
    inter<-intersect(test$SNP,testout$SNP)
    test <- test[test$SNP %in% inter,]
    testout <- testout[testout$SNP %in% inter,]
    testharm <- harmonise_data(test, testout, action=2)
    testharm$mr_keep<-TRUE
    testres<-mr(testharm)
    ###steiger filtering
    testres <- generate_odds_ratios(testres)
    testres <- cbind(SNP=1,testres,beta.exp=1,beta.out=1,se.exp=1,se.out=1,eaf.out=1,p.exp=1,p.out=1,samplesize.exposure=1)
    for (r in 1: nrow(testharm)){
      testres[r,1] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],"SNP"]
      testres[r,16:23] <- testharm[testharm$id.exposure==testres[r,"id.exposure"],c("beta.exposure","beta.outcome","se.exposure",
                                                                                    "se.outcome","eaf.outcome",
                                                                                    "pval.exposure","pval.outcome","samplesize.exposure")]
    }
    testres <- na.omit(testres)
    testres <- cbind(testres,rsq.exposure=1,rsq.outcome=1,samplesize.outcome=outcome_samplesize[j])
    testres$rsq.exposure <- (get_r_from_pn(p=testres$p.exp, n=testres$samplesize.exposure))^2
    testres$rsq.outcome <- (get_r_from_lor(lor=log(testres$or),af=testres$eaf.out,ncase=ncase[j],ncontrol=ncontrol[j],prevalence=ncase[j]/outcome_samplesize[j]))^2
    st <- psych::r.test(
      n = testres$samplesize.exposure, 
      n2 = testres$samplesize.outcome, 
      r12 = sqrt(testres$rsq.exposure), 
      r34 = sqrt(testres$rsq.outcome)
    )
    testres$steiger_dir <- testres$rsq.exposure > testres$rsq.outcome
    testres$steiger_pval <- pnorm(-abs(st$z)) * 2
    
    harmdat$match<-paste(harmdat$SNP,harmdat$exposure,sep="_")
    harmdat$steiger_dir <- 0
    harmdat$steiger_pval <- 0
    for (k in 1:nrow(harmdat)){
      harmdat[k,"steiger_dir"] <- testres[testres$exposure==harmdat[k,"match"],"steiger_dir"][1]
      harmdat[k,"steiger_pval"] <- testres[testres$exposure==harmdat[k,"match"],"steiger_pval"][1]
    }
    harmdat<-harmdat[harmdat$steiger_dir==1,]
    length(unique(harmdat$exposure))
    harmname<-paste("harmdata",outcome_name[j],outcome_id[j],sep="_")
    harmname<-paste(harmname,".txt",sep="")
    write.table(harmdat, file=harmname, row.names=F, col.names=T, sep="\t", quote=F)
    
    #Generalized IVW  and MR-Egger for the "MendelianRandomization" package
    exp<-unique(harmdat$exposure)
    exp1=exposure_dat
    exp2=NULL
    mrres <- c()
    for (t in 1:length(exp)){
      dat <- harmdat[harmdat$exposure==exp[t],]
      dat<-dat[!is.na(dat$SNP),]
      if (nrow(dat)==1){
        result <- mr_wald_ratio(b_exp= dat$beta.exposure, b_out=dat$beta.outcome, se_exp= dat$se.exposure, se_out= dat$se.outcome)
        result <- data.frame(id=exp[t],method="Wald_ratio",nsnp=result$nsnp,b=result$b,se=result$se,CIlower=NA,CIupper=NA,
                             pval=result$pval,intercept=NA,intercept_se=NA,inter_CIlower=NA,inter_CIupper=NA,
                             intercept_pval=NA,hetero_Q=NA,hetero_pvale=NA)}
      else if (nrow(dat)==2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
      result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                        bxse = dat$se.exposure,
                                                        by = dat$beta.outcome,
                                                        byse = dat$se.outcome,
                                                        cor = rho))
      result <- data.frame(id=exp[t],
                           method="IVW",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError,
                           CIlower=result$CILower,
                           CIupper=result$CIUpper,
                           pval=result$Pvalue,
                           intercept=NA,
                           intercept_se=NA,
                           inter_CIlower=NA,
                           inter_CIupper=NA,
                           intercept_pval=NA,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])}
      else if (nrow(dat)==2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                                                                            bxse = dat$se.exposure,
                                                                                                            by = dat$beta.outcome,
                                                                                                            byse = dat$se.outcome))
      result <- data.frame(id=exp[t],
                           method="IVW",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError,
                           CIlower=result$CILower,
                           CIupper=result$CIUpper,
                           pval=result$Pvalue,
                           intercept=NA,
                           intercept_se=NA,
                           inter_CIlower=NA,
                           inter_CIupper=NA,
                           intercept_pval=NA,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])}
      else if (nrow(dat)>2 & all(dat$SNP %in% exp1$SNP)){rho <- ld_matrix(dat$SNP, pop = "EUR")
      result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                          bxse = dat$se.exposure,
                                                          by = dat$beta.outcome,
                                                          byse = dat$se.outcome,
                                                          cor = rho))
      result <- data.frame(id=exp[t],
                           method="MR-Egger",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError.Est,
                           CIlower=result$CILower.Est,
                           CIupper=result$CIUpper.Est,
                           pval=result$Pvalue.Est,
                           intercept=result$Intercept,
                           intercept_se=result$StdError.Int,
                           inter_CIlower=result$CILower.Int,
                           inter_CIupper=result$CIUpper.Int,
                           intercept_pval=result$Pvalue.Int,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])
      result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                         bxse = dat$se.exposure,
                                                         by = dat$beta.outcome,
                                                         byse = dat$se.outcome,
                                                         cor = rho))
      result2 <- data.frame(id=exp[t],
                            method="IVW",
                            nsnp=result2$SNPs,
                            b=result2$Estimate,
                            se=result2$StdError,
                            CIlower=result2$CILower,
                            CIupper=result2$CIUpper,
                            pval=result2$Pvalue,
                            intercept=NA,
                            intercept_se=NA,
                            inter_CIlower=NA,
                            inter_CIupper=NA,
                            intercept_pval=NA,
                            hetero_Q=result2$Heter.Stat[1],
                            hetero_pvale=result2$Heter.Stat[2])
      result <- rbind(result,result2)}
      else if (nrow(dat)>2 & any(dat$SNP %in% exp2$SNP)){result <- MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure,
                                                                                                             bxse = dat$se.exposure,
                                                                                                             by = dat$beta.outcome,
                                                                                                             byse = dat$se.outcome))
      result <- data.frame(id=exp[t],
                           method="MR-Egger",
                           nsnp=result$SNPs,
                           b=result$Estimate,
                           se=result$StdError.Est,
                           CIlower=result$CILower.Est,
                           CIupper=result$CIUpper.Est,
                           pval=result$Pvalue.Est,
                           intercept=result$Intercept,
                           intercept_se=result$StdError.Int,
                           inter_CIlower=result$CILower.Int,
                           inter_CIupper=result$CIUpper.Int,
                           intercept_pval=result$Pvalue.Int,
                           hetero_Q=result$Heter.Stat[1],
                           hetero_pvale=result$Heter.Stat[2])
      result2 <- MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure,
                                                         bxse = dat$se.exposure,
                                                         by = dat$beta.outcome,
                                                         byse = dat$se.outcome))
      result2 <- data.frame(id=exp[t],
                            method="IVW",
                            nsnp=result2$SNPs,
                            b=result2$Estimate,
                            se=result2$StdError,
                            CIlower=result2$CILower,
                            CIupper=result2$CIUpper,
                            pval=result2$Pvalue,
                            intercept=NA,
                            intercept_se=NA,
                            inter_CIlower=NA,
                            inter_CIupper=NA,
                            intercept_pval=NA,
                            hetero_Q=result2$Heter.Stat[1],
                            hetero_pvale=result2$Heter.Stat[2])
      result <- rbind(result,result2)}
      mrres <- rbind(mrres,result)
    }
    #####FDR Adjustment for Pvalue######
    mrres$fdr <- p.adjust(mrres$pval, method = "fdr", n = length(mrres$pval))
    resname<-paste("mrres",outcome_name[j],outcome_id[j],sep="_")
    resname<-paste(resname,".txt",sep="")
    write.table(mrres, file=resname, row.names=F, col.names=T, sep="\t", quote=F)
  }






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

setwd("/Users/tq20202/Desktop/DICE")
res<-read.csv("fdrres.csv")
setwd("/Users/tq20202/Desktop/DICE/raw data")
bres<-c()
for (i in 1: nrow(res)){
  snp<-res[i,"SNP"]
  chr<-res[i,"chr"]
  pos<-res[i,"pos"]
  cell<-res[i,"cell"]
  cellfile<-res[i,"cellfile"]
  outid<-res[i,"outcome_id"]
  f<-read.table(file=cellfile, header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F)
  data<-f[f$chr==chr & f$position>pos-500000 & f$position<pos+500000,]
  data$se<-abs(data$beta)/qnorm(1-data$p/2)
  data$se<-data$se+0.00001
  data$n=1544
  exp=data.frame(SNP=data$SNP,effect_allele=data$ALT,other_allele=data$REF,effect_allele_freq=0,beta=data$beta,se=data$se,p=data$p,n=data$n)
  out <- extract_outcome_data(snps=exp$SNP, outcomes= outid)
  temp=intersect(exp$SNP,out$SNP)
  exp<-exp[exp$SNP %in% temp,]
  exp <- exp[order(exp$SNP),]
  exp <- exp[!duplicated(exp$SNP),]
  out<-out[out$SNP %in% temp,]
  out <- out[order(out$SNP),]
  out <- out[!duplicated(out$SNP),]
  for (j in 1:nrow(exp)){
    e1<-exp[j,"effect_allele"]
    e2<-out[j,"effect_allele.outcome"]
    if (e1==e2){exp[j,"effect_allele_freq"]=out[j,"eaf.outcome"]}else{exp[j,"effect_allele_freq"]=1-out[j,"eaf.outcome"]}
  }
  out=data.frame(SNP=out$SNP,effect_allele=out$effect_allele.outcome,other_allele=out$other_allele.outcome,effect_allele_freq=out$eaf.outcome,beta=out$beta.outcome,se=out$se.outcome,p=out$pval.outcome,n=outcome_samplesize[which(outcome_id==outid)],case=ncase[which(outcome_id==outid)])
  expfile<-paste("/Users/tq20202/Desktop/DICE/pwcoco/data/",i,"_exp.txt",sep="")
  outfile<-paste("/Users/tq20202/Desktop/DICE/pwcoco/data/",i,"_out.txt",sep="")
  write.table(exp, file=expfile, row.names=F, col.names=T, sep="\t", quote=F)
  write.table(out, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
  df <- data.frame(SNP=exp$SNP,beta1=as.numeric(exp$beta),beta2=as.numeric(out$beta),se1=as.numeric(exp$se),se2=as.numeric(out$se),MAF1=as.numeric(exp$effect_allele_freq),MAF2=as.numeric(out$effect_allele_freq),N1=exp$n,N2=out$n,s=out$case/out$n)
  df <- na.omit(df)
  #coloc analysis
  result <- coloc.analysis(df$beta1, df$beta2, df$se1, df$se2, df$MAF1, df$MAF2, N1=df$N1, N2=df$N2, s=df$s) 
  result <- data.frame(result,num=i)
  bres <- rbind(bres,result)
}
write.table(bres, file="/Users/tq20202/Desktop/DICE/coloc_res.txt", row.names=F, col.names=T, sep="\t", quote=F)












