#R_code (v4.2.0) for Figure 2b
#load all package
library(pheatmap)
setwd("Fig. 2/")
#AAI value (genome_shared_gene >= 10%)
data <- read.delim('aai.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
matrix_data<-reshape(data, idvar = "qname", timevar = "tname", direction = "wide")
colnames(matrix_data) <- sub("aai.", "", colnames(matrix_data))
matrix_data[is.na(matrix_data)] <- 0
write.csv(matrix_data, file = "matrix.csv", row.names = FALSE)
data<-read.csv('matrix.csv',row.names = 1)
annotation_row<- read.delim("Taxonomy.txt",header=T,sep="\t",row.names=1)                                                                                         
annotation_row <- as.data.frame(annotation_row)
ann_colors = list(Category=c(Chitovirales="#FFD5D5",Asfuvirales="#FFF3F3",Pimascovirales="#FEF3E6",
                             Pandoravirales="#F2EFF8",Algavirales="#EDF6FF",Imitervirales="#D9EAD3"))
color_palette <- colorRampPalette(c("#F9FCFD", "#A2D6F2", "#004c8c"))(100)
pheatmap(data,scale = "none",color = color_palette ,breaks = seq(0, 100),
         cluster_cols=F,cluster_row=F,annotation_colors = ann_colors,
         annotation_row=annotation_row,
         show_colnames=F,show_rownames=F,fontsize= 11)

