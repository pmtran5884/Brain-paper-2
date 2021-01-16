
rm(list=ls())
library(readxl)
library(RColorBrewer)
library(pheatmap)
setwd("/Users/paultran/Box/Brain paper 2/code")
mypheno<-read_excel("Brain2_Supp_tables.xlsx", sheet = "Supplementary Table 3", skip = 1)
mygeno<-read_excel("Brain2_Supp_tables.xlsx", sheet = "Supplementary Table 4", skip = 1)

#table WHO2016
table(mypheno$WHO2016)
table(mypheno$WHO2016,mypheno$Sup168)

#match rownames
mygeno<-mygeno[match(mypheno$Nanostring.or.submitted.ID,mygeno$...1),-1]
rownames(mygeno)<-mypheno$Nanostring.or.submitted.ID

#define heatamp colors
int <- 0.01
pairs.breaks <- c(seq(-1, -int, by=0.0005), seq(-int,int, by=0.1), seq(int, 1, by=0.0005))
colfunc<-colorRampPalette(c("blue", "white", "red"))
mycol<-colfunc(length(pairs.breaks)-1)


#make gene cols for KI67, GFAP, TP53, Codel, IDH, Preservation, WHO2016, and Sup168
my_samp_col <- mypheno[,c(4,5,7:12)]
my_samp_col <- data.frame(apply(my_samp_col,2,as.factor))
my_samp_col <- droplevels.data.frame(my_samp_col)
colnames(my_samp_col)[1]<-"TP"
rownames(my_samp_col) <- mypheno$Nanostring.or.submitted.ID

#choose colors for everything else
mut_cols<-c(`1`="#000000",
            `0`="grey",
            `NA`="#FFFFFF")
annotation_colors<- list(
  # IDH_codel_subtype = c(`IDHwt`="#E41A1C",
  #                       `IDHmut-codel`="#377EB8",
  #                       `IDHmut-non-codel`="#4DAF4A",
  #                       `Unknown`="#FFFFFF"),
  #WHO2016 = c(brewer.pal(9,"Set1"),"grey"),
  TP = c(`TP1`="#984EA3",
                         `TP2`="#FF7F00",
                         `TP3`="#F781BF",
                         `TP4`="#000000",
                         `Ambi`="grey"),
  Preservation = c(FFPE = "pink",
                   Frozen = "cyan"),
  KI67 = mut_cols,
  GFAP = mut_cols,
  TP53 = mut_cols,
  Codel = mut_cols,
  IDH = mut_cols

)

#make heatmap
pheatmap(as.matrix(t(mygeno)),
                       annotation_col = my_samp_col,
                       annotation_colors = annotation_colors,
                       scale = "none",color = mycol,
                       cluster_rows = T,cluster_cols = F,
                       fontsize_row = 2, show_colnames = F,
                       clustering_method = "ward.D2",
                       clustering_distance_rows = "correlation",
                       na_col = "#808080",annotation_legend = F,
                       width = 8,height = 10,
                       border_color = NA,
                       filename = "09_Brain_MCG_heatmap.png",
                       res=300)

pheatmap(as.matrix(t(mygeno)),
         annotation_col = my_samp_col,
         annotation_colors = annotation_colors,
         scale = "none",color = mycol,
         cluster_rows = T,cluster_cols = F,
         fontsize_row = 2, show_colnames = F,
         clustering_method = "ward.D2",
         clustering_distance_rows = "correlation",
         na_col = "#808080",annotation_legend = T,
         width = 8,height = 10,
         filename = "09_Brain_MCG_heatmap_withanno.png",
         res=300)
