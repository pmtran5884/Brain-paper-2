
# Analysis Date: August 6, 2019
# Author: Paul Tran
# Title: Brain Ca Nanostring Replicability

"pseudocode
read in data
identify all samp names with duplicates
code each duplicate as day to day instrumental replicate vs different core replicate based on samp name
plot all replicates
1)correlation heatmap with clustering
2)pairs plots
3)average correlation for replicates, vs average correlation non-reps, dot/boxplot"

rm(list=ls(all=TRUE))
require(tidyverse)
require(corrplot)

setwd("C:/Users/ptran/Box Sync/TCGA Cancer Classification Project/Brain ns processing 20190802")
mydata<-read.csv("01Brain_Nanostring_Data_PostQC_20190805_withHK.csv",row.names = 1)
glimpse(mydata)

colnames(mydata)
oldsamps_ind<-c(1:34,65,69,70,71)
colnames(mydata)[oldsamps_ind]
colnames(mydata)[-oldsamps_ind]


subj_ids<-c(unlist(lapply(strsplit(colnames(mydata),".",fixed = T), function(l) l[[2]])))
pheno<-cbind.data.frame(colnames(mydata),subj_ids)
pheno$`colnames(mydata)`<-as.vector(pheno$`colnames(mydata)`)
pheno$subj_ids<-as.vector(pheno$subj_ids)


samps_withreps<-unique(pheno$subj_ids[duplicated(pheno$subj_ids)])
pheno_reps<-pheno[pheno$subj_ids%in%samps_withreps,]
pheno_reps<-pheno_reps[order(pheno_reps$subj_ids),]
# write.csv(pheno_reps,"02BrainNs_withreps.csv",row.names = F)

pheno_reps_anno<-read.csv("02BrainNs_withreps_annotated.csv")

####day-to-day/same core/nanostring replicability#####
samp_names_rep_dd<-as.vector(pheno_reps_anno$colnames.mydata.[pheno_reps_anno$replicate_type..dd.=="y"])
corrplot_samp_names_rep_dd<-as.vector(pheno_reps_anno$corrplot_samp_ids[pheno_reps_anno$replicate_type..dd.=="y"])
samp_rep_dd<-mydata[,match(samp_names_rep_dd,colnames(mydata))]
cbind(colnames(samp_rep_dd),corrplot_samp_names_rep_dd)
colnames(samp_rep_dd)<-corrplot_samp_names_rep_dd
M<-cor(samp_rep_dd)
# cor.mtest <- function(mat, ...) {
#   mat <- as.matrix(mat)
#   n <- ncol(mat)
#   p.mat<- matrix(NA, n, n)
#   diag(p.mat) <- 0
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       tmp <- cor.test(mat[, i], mat[, j], ...)
#       p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#     }
#   }
#   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#   p.mat
# }
# matrix of the p-value of the correlation
p.mat <- cor.mtest(samp_rep_dd)
# corrplot(M, method="circle", addCoef.col = "black",p.mat = p.mat, sig.level = 0.01,insig = "blank")
# corrplot.mixed(M,lower.col = "black",number.cex=0.7,tl.pos = "d",tl.cex=0.5,
#                tl.col="black",tl.srt=45)
corrplot(M,number.cex=0.7,tl.pos = "lt",tl.cex=1,
               tl.col="black",tl.srt=45,order="hclust",addrect=6)

# pairs(samp_rep_dd)
# Correlation panel
# panel.cor <- function(x, y){
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- round(cor(x, y), digits=2)
#   txt <- paste0("R = ", r)
#   cex.cor <- 2#/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor )#* r)
# }
# 
# # Create the plots
# pairs(samp_rep_dd, 
#       lower.panel = panel.cor)

cbind.data.frame(colnames(M),colnames(M),M[lower.tri(M,diag = T)])
data.frame(row=rownames(M)[row(M)], col=colnames(M)[col(M)], corr=c(M))

#pretty figures
pairs(samp_rep_dd)

colnames(samp_rep_dd)


samp_rep_dd %>% gather(Samp1, Rep1, c(BRAIN17B_1,BRAIN19_1,BRAIN205_1,BRAIN207_1,BRAIN33_1,BRAIN37_1))->combined_data1
combined_data1 %>% gather(Samp2, Rep2, c(BRAIN17B_2,BRAIN19_2,BRAIN205_2,BRAIN207_2,BRAIN33_2,BRAIN37_2))->combined_data2

Samp1<-unlist(lapply(strsplit(combined_data2$Samp1,"_"), `[[`, 1))
Samp2<-unlist(lapply(strsplit(combined_data2$Samp2,"_"), `[[`, 1))
length(which(Samp1==Samp2))

combined_data3<-combined_data2[which(Samp1==Samp2),]
ggplot(data=combined_data3,aes(x=Rep1, y=Rep2))+
  geom_point()+
  facet_wrap(~Samp1)+theme_classic2()
ggsave("SuppFig3_NSreplicability_scatterplots.png",
       width= 6,
       height = 4,
       units = "in")
#######different core/sample replicability########
samp_names_rep_core<-as.vector(pheno_reps_anno$colnames.mydata.[pheno_reps_anno$replicate_type..core.=="y"])
corrplot_samp_names_rep_core<-as.vector(pheno_reps_anno$corrplot_samp_ids[pheno_reps_anno$replicate_type..core.=="y"])
samp_rep_core<-mydata[,match(samp_names_rep_core,colnames(mydata))]
cbind(colnames(samp_rep_core),corrplot_samp_names_rep_core)
colnames(samp_rep_core)<-corrplot_samp_names_rep_core
N<-cor(samp_rep_core)
# corrplot(N, method="circle")
corrplot(N,number.cex=0.7,tl.pos = "lt",tl.cex=1,
         tl.col="black",tl.srt=45,order="hclust",addrect=6)
# pairs(samp_rep_core, 
#       lower.panel = panel.cor)
