rm(list=ls())
gc()
library(TwoSampleMR)
library(MendelianRandomization)
setwd("/Users/tq20202/Desktop/oneK1K")
fname<-list.files("data")
setwd("/Users/tq20202/Desktop/oneK1K/data")
data<-c()
for (i in 1:length(fname)){
  exposure_dat<-read.csv(fname[i])
  exposure_dat$beta<-exposure_dat$SPEARMANS_RHO/sqrt(2*exposure_dat$A2_FREQ_ONEK1K*(1-exposure_dat$A2_FREQ_ONEK1K))
  exposure_dat$z<-qnorm(exposure_dat$P_VALUE,lower.tail = F)
  exposure_dat$se <- abs(exposure_dat$beta/exposure_dat$z)
  exposure_dat$se=exposure_dat$se+0.00001
  exposure_dat<-exposure_dat[exposure_dat$P_VALUE<=5e-8,]
  exposure_dat$pve<-exposure_dat$z^2/(982+exposure_dat$z^2)
  exposure_dat$fstatistics<-(982-2)*exposure_dat$pve/(1-exposure_dat$pve)
  exposure_dat<-exposure_dat[exposure_dat$fstatistics>10,]
  tmp <- data.frame(SNP = exposure_dat$RSID,beta = exposure_dat$beta,se = exposure_dat$se,
                    effect_allele = exposure_dat$A2,eaf = exposure_dat$A2_FREQ_ONEK1K,
                    other_allele = exposure_dat$A1,pval = exposure_dat$P_VALUE,
                    chr = exposure_dat$CHR,pos = exposure_dat$POS,samplesize = 982,
                    phenotype_id = paste(exposure_dat$GENE,exposure_dat$CELL_TYPE,sep="_"))
  exposure_dat <- format_data(tmp, type="exposure",effect_allele_col = "effect_allele",other_allele_col = "other_allele",phenotype_col = "phenotype_id",samplesize_col = "samplesize",chr_col = "chr",pos_col = "pos")
  exposure_dat <-clump_data(exposure_dat,clump_kb=10000,clump_r2=0.001,clump_p1=1,clump_p2=1,pop="EUR")
  data<-rbind(data,exposure_dat)
} 
for(i in 1:nrow(data)){
  if (data[i,"eaf.exposure"]>0.5){oa<-data[i,"effect_allele.exposure"];ea<-data[i,"other_allele.exposure"];data[i,"other_allele.exposure"]<-oa;data[i,"effect_allele.exposure"]<-ea;data[i,"eaf.exposure"]=1-data[i,"eaf.exposure"]}
}
write.csv(data,"/Users/tq20202/Desktop/oneK1K/exposure.csv")

exposure_dat<-read.csv("/Users/tq20202/Desktop/oneK1K/exposure.csv")
setwd("/Users/tq20202/Desktop/oneK1K/mrres")
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
for (j in 1:length(outcome_id)){
  outcome <- extract_outcome_data(snps=exposure_dat$SNP, outcomes= outcome_id[j])
  outcome<-outcome[!duplicated(outcome$SNP),]
  inter<-intersect(exposure_dat$SNP,outcome$SNP)
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% inter,]
  outcome_dat <- outcome[outcome$SNP %in% inter,]
  harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)
  harmdat$mr_keep<-TRUE  #consider ambiguous
  
  test<-harmdat
  test$new<- paste(test$SNP,test$exposure,sep="_")
  test<-test[,-40]
  test <- format_data(test, type="exposure", snp_col = "SNP", pval_col = "pval.exposure",
                      beta_col = "beta.exposure", se_col = "se.exposure",
                      effect_allele_col = "effect_allele.exposure",
                      eaf_col = "eaf.exposure",
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
  mrres<-generate_odds_ratios(mrres)
  mrres$fdr <- p.adjust(mrres$pval, method = "fdr", n = length(mrres$pval))
  resname<-paste("mrres",outcome_name[j],outcome_id[j],sep="_")
  resname<-paste(resname,".txt",sep="")
  resname2<-paste("harmdat",outcome_name[j],outcome_id[j],sep="_")
  resname2<-paste(resname2,".txt",sep="")
  write.table(mrres, file=resname, row.names=F, col.names=T, sep="\t", quote=F)
  write.table(harmdat, file=resname2, row.names=F, col.names=T, sep="\t", quote=F)
}

#####FDR result########
setwd("/Users/tq20202/Desktop/oneK1K")
fname<-list.files("mrres")
setwd("/Users/tq20202/Desktop/oneK1K/mrres")
fdrres<-c()
for (i in 1:length(fname)){
  t<-read.table(file=fname[i], header=T, sep="\t", quote="", stringsAsFactors=F, check.names=F) 
  t<-t[t$fdr<0.05,]
  if(nrow(t)>0){t$cancer<-fname[i]}
  fdrres<-rbind(fdrres,t)
}
write.csv(fdrres,"/Users/tq20202/Desktop/oneK1K/fdrres.csv")









                           