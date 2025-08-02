#R_code (v4.2.0) for Figure 2
#load all package
library(pheatmap)
#work_dir
setwd("Fig. 2/")
#AAI value (genome_shared_gene >= 10%)
data <- read.delim('aai.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#matrix
matrix_data<-reshape(data, idvar = "qname", timevar = "tname", direction = "wide")
colnames(matrix_data) <- sub("aai.", "", colnames(matrix_data))
matrix_data[is.na(matrix_data)] <- 0
write.csv(matrix_data, file = "matrix.csv", row.names = FALSE)
#plot
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

#R_code (v4.2.0) for Figure S1
#load all package
library(ggplot2)
library(readxl)
setwd("Fig. S1/")
data<-read_excel('Fig. S1.xlsx')
data$Genome <- factor(data$Genome, levels = unique(data$Genome))
ggplot(data)+
  geom_bar(stat="identity",width=0.75, position="dodge", aes(x=Genome, y=markers),alpha=0.8,fill="#E0E0E0",size=0.2)+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(color="black"))+
  theme_test()+theme(axis.title.x = element_text(size=8),
                     axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.y = element_blank(),
                     axis.ticks = element_line(size = 0.25),
                     axis.ticks.length = unit(0.04, "cm"),
                     panel.border = element_rect(size=0.4))+
  theme(legend.position='none')+
  coord_flip()+ylim(0,11)
ggsave("bar.pdf", width= 1,height= 8) 

#R_code (v4.2.0) for Figure S5
#load all package
library(pheatmap)
setwd("Fig. S5/")
#plot
data <- read.delim('AAI_imi.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
matrix_data<-reshape(data, idvar = "qname", timevar = "tname", direction = "wide")
colnames(matrix_data) <- sub("aai.", "", colnames(matrix_data))
matrix_data[is.na(matrix_data)] <- 0
write.csv(matrix_data, file = "matrix_imi.csv", row.names = FALSE)

data<-read.csv('matrix_imi.csv',row.names = 1)
color_palette <- colorRampPalette(c("#FCFCFC", "#BEBEBE", "#3C3C3C"))(100)
pheatmap(data,scale = "none",color = color_palette ,breaks = seq(0, 100),
         cluster_cols=F,cluster_row=F,
         show_colnames=F,show_rownames=F,fontsize= 11,border_color = NA)

