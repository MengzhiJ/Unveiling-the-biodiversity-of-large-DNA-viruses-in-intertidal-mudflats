#R_code (v4.2.0) for Figure 3
#load all package
library(ggplot2)
library(reshape2)
library(dplyr)
library(readxl)
library(tidyverse)
library(scales)
library(ggforce)
library(ggsci)
library(ape)
library(picante)

setwd("Fig. 3/")
#number of host
data<-read_excel("Fig. 3.xlsx")
data$Virus_tax <- factor(data$Virus_tax,levels=c("Biggiephage","Jabbarphage","Judaphage","Whopperphage","Kabirphage","Mahaphage","Unclassified"))
data$Host_tax <- factor(data$Host_tax,levels=c('Planctomycetota', 'Nitrospirota','Chloroflexota','Thermotogota','Patescibacteria', 'Actinobacteriota','Cyanobacteria','Firmicutes','Bacteroidota','Proteobacteria'))
ggplot(data, aes(x=Host_tax, y=Number))+
  geom_bar(stat="identity", position="stack", aes(fill=Virus_tax),alpha=0.8,width=0.8)+
  scale_fill_manual(breaks=c("Biggiephage","Jabbarphage","Judaphage","Whopperphage","Kabirphage","Mahaphage"),
    values = c("#78be59","#f5c533","#CA8EFF", "#f6b26b","#9cceff","#2894ff"))+
  labs(x="",y="# of viral genomes")+
  theme(legend.background=element_rect(color="white"),)+
  theme(legend.key = element_blank())+
  theme(axis.line = element_line(color=))+
  theme_test()+coord_flip()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))+
  ylim(0,30)+theme(legend.position = "none")

#plot phage number
sales <- c(79.5,20.5)
names<-c("a","b")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),
               stat="pie",alpha=0.9)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("a","b"),values=c("#7AC5CD","#f6b26b"))

#calculate phage PD
tree<-read.tree(file="Terl_trim.faa.contree")
edges <- tree$edge
edge_lengths <- tree$edge.length
tip_labels <- tree$tip.label

get_node_label <- function(node) {
  if (node <= length(tip_labels)) {
    return(tip_labels[node])
  } else {
    return(paste("Node", node, sep = "_"))
  }
}

edge_info <- data.frame(
  Start = sapply(edges[, 1], get_node_label),
  End = sapply(edges[, 2], get_node_label),
  Length = edge_lengths
)

#plot phage PD
sales <- c(38.3,36.6,25.1)
names<-c("a","b","c")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0,0),fill=names),stat="pie",alpha=0.9)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("a","b","c"),values=c("#7AC5CD","#E0E0E0","#f6b26b"))

#plot pd_host_tax
data<-read_excel("Fig. 3.xlsx")
data$Host_tax <- factor(data$Host_tax,levels=c('Planctomycetota', 'Nitrospirota','Chloroflexota',
                                               'Thermotogota','Patescibacteria',
                                               'Actinobacteriota','Cyanobacteria','Firmicutes',
                                               'Bacteroidota','Proteobacteria'))
ggplot(data, aes(x=Host_tax, y=PD))+
  geom_bar(stat="identity", position="stack", fill="#d0d0d0",alpha=0.8,width=0.8)+
  labs(x="",y="Phylogenetic diversity (PD)")+
  theme_test()+coord_flip()+
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        legend.text = element_text(size=12),
         legend.title = element_text(size=12))+ylim(0,15)+
  theme(legend.position = "none")+
  theme(axis.ticks.y.left = element_blank())

#R_code (v4.2.0) for Figure S2
#load all package
library(ggplot2)
library(pheatmap)
library(readxl)

setwd("Fig. S2/")

data<-read.csv('tRNA.csv',row.names = 1)
pheatmap(data,scale="none",
         cluster_cols=F,cluster_row=F,
         colorpanel(70,low="white",mid="#97CBFF",high="#0072E3"),
         legend=T,show_colnames=T,show_rownames=T,fontsize= 8,
         border_color = "#E0E0E0", cellwidth = 10, cellheight = 10,filename = "heatmap.pdf")

data<-read_excel('Number.xlsx')
ggplot(data)+
  geom_bar(stat="identity",width=0.75, position="dodge", aes(x=Genome, y=Number),alpha=0.8,fill="#E0E0E0",size=0.2)+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(color="black"))+
  theme_test()+theme(panel.border = element_rect(size=0.6))+
  theme(legend.position='none')+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank() 
  ) +coord_flip()+ylim(0,70)
ggsave("bar.pdf", width= 1,height= 7) 

#R_code (v4.2.0) for Figure S3
#load all package
library(pheatmap)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(readxl)
library(splines)

setwd("Fig. S3/")
#Fig. S3
#AAI
data <- read.delim('aai.tsv', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
matrix_data<-reshape(data, idvar = "qname", timevar = "tname", direction = "wide")
colnames(matrix_data) <- sub("aai.", "", colnames(matrix_data))
matrix_data[is.na(matrix_data)] <- 0
write.csv(matrix_data, file = "matrix.csv", row.names = FALSE)
data<-read.csv('matrix.csv',row.names = 1)
pheatmap(data,scale = "none", colorpanel(70,low="#FCFCFC",mid="#BEBEBE",high="#3C3C3C"),
         cluster_cols=F,cluster_row=F,
         border_color = "#3C3C3C", show_colnames=T,show_rownames=T,
         fontsize= 13,cellwidth = 20, cellheight = 20,filename = "heatmap.pdf")
#SNV_time
data<-read_excel("Fig. S3.xlsx")
data$Time <- factor(data$Time,levels=c("June_2020","Aug_2020","Oct_2020","Dec_2020","Feb_2021","April_2021","June_2021","July_2023"))
ggplot(data, mapping = aes(x=Time,y=SNV_count))+
  geom_boxplot(aes(group=Time),color="black",width=0.4,outlier.color="white",
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Time,y=SNV_count),fill="gray",width=0.15,shape=21,size=1.3,alpha=0.7)+
  labs(x="Sampling time",y="SNV count")+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#7AC5CD',fill="#7AC5CD",
              alpha=0.2,aes(group=1),formula = y ~ bs(x, df = 4))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=11.5),axis.title.y = element_text(size=11.5),
        axis.text.x = element_text(hjust =0.5,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        panel.border = element_rect(size=1.2))

#Coverage_time
data<-read_excel("Fig. S3.xlsx")
data$Time <- factor(data$Time,levels=c("June_2020","Aug_2020","Oct_2020","Dec_2020","Feb_2021","April_2021","June_2021","July_2023"))
ggplot(data, mapping = aes(x=Time,y=Coverage))+
  geom_boxplot(aes(group=Time),color="black",width=0.4,outlier.color="white",
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Time,y=Coverage),fill="gray",width=0.15,shape=21,size=1.3,alpha=0.7)+
  labs(x="Sampling time",y="Normalized coverage")+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color='#7AC5CD',fill="#7AC5CD",
              alpha=0.2,aes(group=1),formula = y ~ bs(x, df = 4))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=11.5),axis.title.y = element_text(size=11.5),
        axis.text.x = element_text(hjust =0.5,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        panel.border = element_rect(size=1.2))

#SNV_depth
data<-read_excel("Fig. S3.xlsx")
data$Depth <- factor(data$Depth,levels=c("0-20","20-40","40-60","60-80","80-100"))
ggplot(data, mapping = aes(x=Depth,y=SNV_count))+
  geom_boxplot(aes(group=Depth),color="black",width=0.4,outlier.color="white",
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Depth,y=SNV_count),fill="gray",width=0.15,shape=21,size=1.3,alpha=0.7)+
  labs(x="Sampling depth",y="SNV count")+
  geom_smooth(method = 'gam',se=T,size=1,fullrange=T,color='#7AC5CD',fill="#7AC5CD",
              alpha=0.2,aes(group=1),formula = y ~ s(x, k = 3))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=11.5),axis.title.y = element_text(size=11.5),
        axis.text.x = element_text(hjust =0.5,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        panel.border = element_rect(size=1.2))+ylim(20,300)

#Coverage_depth
data<-read_excel("Fig. S3.xlsx")
data$Depth <- factor(data$Depth,levels=c("0-20","20-40","40-60","60-80","80-100"))
ggplot(data, mapping = aes(x=Depth,y=Coverage))+
  geom_boxplot(aes(group=Depth),color="black",width=0.4,outlier.color="white",
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Depth,y=Coverage),fill="gray",width=0.15,shape=21,size=1.3,alpha=0.7)+
  labs(x="Sampling depth",y="Normalized Coverage")+
  geom_smooth(method = 'gam',se=T,size=1,fullrange=T,color='#7AC5CD',fill="#7AC5CD",
              alpha=0.2,aes(group=1),formula = y ~ s(x, k = 3))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=11.5),axis.title.y = element_text(size=11.5),
        axis.text.x = element_text(hjust =0.5,size=10,colour = 'black'),
        axis.text.y=element_text(size=10,colour = 'black'),
        panel.border = element_rect(size=1.2))
