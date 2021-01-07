

# Analysis Date: September 4, 2019
# Author: Paul Tran
# Title: MCG survival and update combined pheno table

rm(list=ls(all=TRUE))
library(survival)
library(survminer)
############Section 1: read in TCGA#############
setwd("/Users/paultran/Box/TCGA Cancer Classification Project/Brain ns processing 20190814")
#read in data
morepheno<-read.csv("05Brain_MCGandTCGA_PhenoData_20190823v4.csv",row.names = 1)
ncounter_pheno<-read.csv("02bBrain_PhenoData_20190814_matched.csv")

#rename group 2b to 4
morepheno$Unsup.DBU<-as.vector(morepheno$Unsup.DBU)
morepheno$Unsup.DBU[morepheno$Unsup.DBU=="Grp2b"]<-"Grp4"
morepheno$Unsup.DBU[morepheno$Unsup.DBU=="Grp2a"]<-"Grp2"
morepheno$Unsup.DBU<-as.factor(morepheno$Unsup.DBU)
table(morepheno$Unsup.DBU)
table(morepheno$Sup168.lSVC.call)
morepheno$Sup168.lSVC.call<-gsub("Grp2a","Grp2",morepheno$Sup168.lSVC.call)
morepheno$Sup168.lSVC.call<-gsub("Grp2b","Grp4",morepheno$Sup168.lSVC.call)
table(morepheno$Sup168.lSVC.call)

#simplify sup calls
morepheno$Sup168.lSVC.call.simp<-as.vector(morepheno$Sup168.lSVC.call)
morepheno$Sup168.lSVC.call.simp[grep(".", morepheno$Sup168.lSVC.call.simp,fixed = T)]<-"Ambi"
table(morepheno$Sup168.lSVC.call.simp)

#calc surv info
ncounter_pheno$Date.of.Initial.Diagnosis[ncounter_pheno$Date.of.Initial.Diagnosis==""]<-NA
ncounter_pheno$Date.of.Last.Contact[ncounter_pheno$Date.of.Last.Contact==""]<-NA
ncounter_pheno$Vital.Status[ncounter_pheno$Vital.Status==""]<-NA
ncounter_pheno$Date.of.Last.Contact<-strptime(ncounter_pheno$Date.of.Last.Contact,format = "%m/%d/%Y")
ncounter_pheno$Date.of.Initial.Diagnosis<-strptime(ncounter_pheno$Date.of.Initial.Diagnosis,format = "%m/%d/%Y")

ncounter_pheno$Survival.months<-(ncounter_pheno$Date.of.Last.Contact-ncounter_pheno$Date.of.Initial.Diagnosis)/30

ncounter_pheno$Vital.Status<-as.vector(ncounter_pheno$Vital.Status)
ncounter_pheno$Vital.Status[grep("Dead",ncounter_pheno$Vital.Status)]<-1
ncounter_pheno$Vital.Status[grep("Alive",ncounter_pheno$Vital.Status)]<-0
ncounter_pheno$Vital.Status<-as.numeric(ncounter_pheno$Vital.Status)

#update pheno tables with MCG surv info
morepheno[morepheno$Study=="MCG",]->MCG
cbind(rownames(MCG),as.vector(ncounter_pheno$Study.ID))
MCG<-cbind.data.frame(MCG,as.numeric(ncounter_pheno$Survival.months),ncounter_pheno$Vital.Status)
colnames(MCG)[c(27,28)]<-c("Survival.months","Vital.Status")

cbind(rownames(MCG),rownames(morepheno)[morepheno$Study=="MCG"])
morepheno$Survival.months[morepheno$Study=="MCG"]<-MCG[,27]
morepheno$Vital.status[morepheno$Study=="MCG"]<-MCG[,28]


#add histo to combined pheno table
pheno<-read.csv("08Brain_MCGandTCGA_PhenoDatav5.csv",row.names = 1)
head(cbind(rownames(pheno),rownames(morepheno)),100)
morepheno$Histology<-c(pheno$Histology,ncounter_pheno$Histo.Behavior.ICD.O.3)

#write new pheno table
write.csv(morepheno,"08Brain_MCGandTCGA_PhenoDatav5.csv")


###MCG survival plots #############
pdf("08_MCGsurvivalplots.pdf")
msurv<-Surv(MCG[,27],MCG[,28])
fit<-survfit(msurv~Sup168.lSVC.call.simp,data=MCG)

ggsurvplot(
  fit,
  pval = T,
  linetype = "strata",
  risk.table = T,
  palette = c("grey",
              "#984EA3",
              "#FF7F00",
              "#F781BF",
              "#000000"),
  xlab="Time (months)"
) 

#rm ambi, combine 2-4
table(MCG$Sup168.lSVC.call.simp)

MCGsub<-MCG[MCG$Sup168.lSVC.call.simp!="Ambi",]
MCGsub$Sup_grps_comb<-MCGsub$Sup168.lSVC.call.simp
MCGsub$Sup_grps_comb[MCGsub$Sup_grps_comb!="Grp1"]<-"notGrp1"
msurv<-Surv(MCGsub[,27],MCGsub[,28])
fit<-survfit(msurv~Sup_grps_comb,data=MCGsub)
fit
summary(coxph(msurv~Sup_grps_comb,data=MCGsub))
ggsurvplot(
  fit,
  pval = F,
  linetype = "strata",
  risk.table = F,
  palette = c("#984EA3",
              "red"),
  newpage=T,
  xlab="Time (months)"
) 
ggsave("../Brain 20190913/04_MCGsupmodel_survival.png",
       width = 3,
       height = 2.5,
       dpi = 300,
       units = "in")
#PGS
MCGsub$PGS[MCGsub$Sup_grps_comb=="notGrp1"]
MCGsub$Supgrp_PGS<-MCGsub$Sup_grps_comb
MCGsub$Supgrp_PGS[MCGsub$Sup_grps_comb=="notGrp1" & MCGsub$PGS == 1]<-"notGrp1_PGShigh"
MCGsub$Supgrp_PGS[MCGsub$Sup_grps_comb=="notGrp1" & MCGsub$PGS == 0]<-"notGrp1_PGSLow"
fit<-survfit(msurv~Supgrp_PGS,data=MCGsub)

ggsurvplot(
  fit,
  pval = T,
  linetype = "strata",
  risk.table = T,
  palette = c("#984EA3",
              "orange",
              "red"),
  newpage=T,
  xlab="Time (months)"
) 

dev.off()
