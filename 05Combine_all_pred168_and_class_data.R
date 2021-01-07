
# Analysis Date: August 20, 2019
# Author: Paul Tran
# Title: Aggregate all prediction and classification data for TCGA and MCG samples

rm(list=ls(all=TRUE))

setwd("C:/Users/ptran/Box Sync/TCGA Cancer Classification Project/Brain ns processing 20190814")

#read all excel sheets function
library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#read in python outputs
tcga_unamb<-read_excel_allsheets("05Combined_validated_3grps_supervised_pred_168_genes_2019_08_22_15_40_47_726538.xlsx")
tcga_amb<-read_excel_allsheets("05Combined_validated_ambiguous_supervised_pred_168_genes_2019_08_22_15_40_47_726538.xlsx")
MCG<-read_excel_allsheets("05Combined_validated_MCG_supervised_pred_168_genes_2019_08_22_15_40_47_726538.xlsx")

#read in pheno data for order
all_pheno<-read.csv("03bBrain_MCGandTCGA_PhenoData_20190821.csv",row.names = 1)

#make combined full pred table
all_pred_Data<-rbind.data.frame(tcga_unamb$`split 1`[1:1003],
                 tcga_unamb$`split 2`[1:1003],
                 tcga_unamb$`split 3`[1:1003],
                 tcga_unamb$`split 4`[1:1003],
                 tcga_amb$`Ambiguous Predictions`[1:1003],
                 MCG$`MCG Predictions`[1:1003])

#reorder to match other tables
all_pred_Data<-all_pred_Data[match(rownames(all_pheno),all_pred_Data$sample_id),]
head(cbind(all_pred_Data$sample_id,rownames(all_pheno)),100)
tail(cbind(all_pred_Data$sample_id,rownames(all_pheno)),100)
rownames(all_pred_Data)<-all_pred_Data$sample_id
all_pred_Data<-all_pred_Data[,-1]

all_counts<-c()
for(i in 1:dim(all_pred_Data)[1]){
  counts<-c(length(which(as.character(all_pred_Data[i,3:1002])=="1")),
            length(which(as.character(all_pred_Data[i,3:1002])=="2a")),
            length(which(as.character(all_pred_Data[i,3:1002])=="2b")),
            length(which(as.character(all_pred_Data[i,3:1002])=="3")))
  all_counts<-rbind.data.frame(all_counts,counts)
  print(paste0("Counting Row #",i))
}

all_counts<-all_counts/1000
colnames(all_counts)<-c("Grp1_freq","Grp2a_freq","Grp2b_freq","Grp3_freq")
cbind.data.frame(all_counts,all_pred_Data)->all_pred_Data
all_pred_Data[1:10,1:10]

#combine with pheno
colnames(all_pheno)
colnames(all_pred_Data[,1:6])
all_pheno<-cbind.data.frame(all_pheno,all_pred_Data[,1:6])
colnames(all_pheno)<-c("Study","Unsup.DBU",         
                       "IDH.codel.subtype","Sup.UMAP", 
                       "Grp1_freq","Grp2a_freq",       
                       "Grp2b_freq","Grp3_freq",        
                       "Sup.Ens.lSVC","Sup.Ens.lSVC.Percent.Predicted")
#write out data,
write.csv(all_pred_Data,"05_lSVC1000_MCGandTCGA_predictions168.csv")
write.csv(all_pheno,"05Brain_MCGandTCGA168_PhenoData_20190822.csv")



table(all_pheno$Unsup.DBU,all_pheno$Sup.UMAP)
table(all_pheno$Unsup.DBU,all_pheno$Sup.Ens.lSVC)
table(all_pheno$Sup.UMAP,all_pheno$Sup.Ens.lSVC)
