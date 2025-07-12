#R_code (v4.2.0) for Figure 4
#load all package
library(Hmisc)
library(psych)
library(data.table)
library(dplyr)
library(vegan)
library(reshape2)
library(tidyverse)
library(scales)
library(ggforce)
library(ggsci)
library(ggpmisc)
library(ggraph)
library(tidygraph)
library(igraph)
library(ggplot2)
library(readxl)
library(ggpubr)

setwd('Fig. 4/')

#correlation between euk_NCLDV lineage
votu<-t(read.delim('huge_mantel.txt',row.names = 1))
votu_bray<-vegdist(votu,method = "bray")
geo<-read.csv('distance.csv',row.names = 1)
geo<-vegdist(geo,method = 'euclidean')
euk<-read.delim('euk_mantel.txt',row.names = 1)

extract_columns <- function(df, pattern) {
  cols_to_extract <- grep(pattern, names(df), value = TRUE)
  if (length(cols_to_extract) == 0) {
    return(data.frame())
  }
  df_extracted <- df %>%
    select(all_of(cols_to_extract))
  return(df_extracted)
}

euk_lineage<- extract_columns(euk, "Amoebozoa")
Amoebozoa<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Annelida")
Annelida<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Arthropoda")
Arthropoda<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Ascomycota")
Ascomycota<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Bacillariophyta")
Bacillariophyta<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Basidiomycota")
Basidiomycota<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Bigyra")
Bigyra<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Bolidophyceae")
Bolidophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Cercozoa")
Cercozoa<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Chlorophyta")
Chlorophyta<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Choanoflagellata")
Choanoflagellata<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Chytridiomycota")
Chytridiomycota<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Ciliophora")
Ciliophora<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Cryptophyceae")
Cryptophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Dictyochophyceae")
Dictyochophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Dinophyceae")
Dinophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Discoba")
Discoba<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Eustigmatophyceae")
Eustigmatophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Fornicata")
Fornicata<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Haptophyta")
Haptophyta<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Microsporidia")
Microsporidia<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Mucoromycota")
Mucoromycota<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Nematoda")
Nematoda<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Pelagophyceae")
Pelagophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Phaeophyceae")
Phaeophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Pinguiophyceae")
Pinguiophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Raphidophyceae")
Raphidophyceae<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Rhodophyta")
Rhodophyta<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Rotifera")
Rotifera<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Streptophyta")
Streptophyta<-vegdist(euk_lineage,method = "bray")
euk_lineage<- extract_columns(euk, "Zoopagomycota")
Zoopagomycota<-vegdist(euk_lineage,method = "bray")

euk_bray<-list(Amoebozoa,Annelida,Arthropoda,Ascomycota,Bacillariophyta,Basidiomycota,Bolidophyceae,Cercozoa,Chlorophyta,Choanoflagellata,
               Chytridiomycota,Ciliophora,Cryptophyceae,Dictyochophyceae,Dinophyceae,Discoba,Eustigmatophyceae,Fornicata,
               Haptophyta,Microsporidia,Mucoromycota,Nematoda,Pelagophyceae,Phaeophyceae,
               Raphidophyceae,Rhodophyta,Rotifera,Streptophyta,Zoopagomycota,Bigyra,Pinguiophyceae)

name<-cbind(Amoebozoa,Annelida,Arthropoda,Ascomycota,Bacillariophyta,Basidiomycota,Bolidophyceae,Cercozoa,Chlorophyta,Choanoflagellata,
            Chytridiomycota,Ciliophora,Cryptophyceae,Dictyochophyceae,Dinophyceae,Discoba,Eustigmatophyceae,Fornicata,
            Haptophyta,Microsporidia,Mucoromycota,Nematoda,Pelagophyceae,Phaeophyceae,
            Raphidophyceae,Rhodophyta,Rotifera,Streptophyta,Zoopagomycota,Bigyra,Pinguiophyceae)
names(euk_bray)<-colnames(name)[1:31]


r<-c()
p<-c()
for (i in 1:31){
  r[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 9999,method="pearson",na.rm=T)$statistic
  p[i]<-mantel.partial(votu_bray,euk_bray[[i]],geo,permutations = 9999,method="pearson",na.rm=T)$signif
}


p<-as.data.frame(p)
r<-as.data.frame(r)
rownames(p)<-colnames(name)
rownames(r)<-colnames(name)

#correlation between NCLDV-euk pairs
combine<-read.delim('combine.txt',row.names = 1)
combine<-as.data.frame(t(combine))
cor<-rcorr(as.matrix(combine),type='spearman')
r<-cor$r
p<-cor$P
r<-r[1:15,16:299]
p<-p[1:15,16:299]
p_adjust<-p.adjust(p,method="bonferroni")
data1=melt(r,value.name = "r")
data2=data.frame(p_adjust)
data<-cbind(data1,data2)
filtered_rows <- data[abs(data$r)>0.8&data$p< 0.05, ]
clean_data<-na.omit(filtered_rows)

#Fig.4a
#bar
data<-read.delim('Fig. 4a.txt',row.names = 1)
ggplot(data, aes(x=correlation_r,y=..density..))+  
  geom_histogram(binwidth = 0.02,alpha=0.3,colour="black",size=0.2,fill="gray")+
  geom_density(alpha=.45,fill="gray",color="gray",lwd=0.22)+
  theme_test()+
  theme(axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
        axis.text.y=element_text(size=7,colour = 'black'),
        axis.ticks = element_line(size = 0.25),
        axis.ticks.length = unit(0.05, "cm"),
        panel.border = element_rect(size=0.45))+
  labs(x="Mantel correlation R",y="Density")+ylim(0,6)+ 
  scale_x_continuous(breaks = seq(-1, 1, by = 0.2))+
  geom_vline(xintercept=c(0.53),linetype=2,color="#66B3FF",size=0.4)

#pie
sales <- c(45.2,54.8)
names<-c("Known","Unknown")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),stat="pie",alpha=0.6)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Known","Unknown"),values=c("#66B3FF","lightgray"))

#Fig.4b
data<-read.delim('Fig. 4b.txt',row.names = 1)
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(fill="lightgray",size = 1.2, alpha = 0.85,shape=21,stroke=0.22)+
  labs(x="Abundance of eukaryotic hosts",y="Abundance of NCLDVs")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#66B3FF",method = 'lm',se=T,size=0.4,fullrange=T,fill="lightgray",alpha=0.25)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+
  ylim(0,6)+xlim(0,8)+
  stat_poly_eq(aes(label = paste(..rr.label.., ..p.value.label.., sep = "~`,`~")),size=2.3,
               label.x = )   

#Fig.4c
node <- read.csv("node.csv")
edge <- read.csv("edge.csv")
network <- tbl_graph(nodes = node, edges = edge, directed = FALSE)
node_colors <- c(
  "Virus" = "#66B3FF",   
  "Algae" = "#4CAF50",   
  "Animal" = "#FF9224",  
  "Fungi" = "#d3a4ff",   
  "Protozoa" = "#FFD306")

edge_colors <- c(
  "Known" = "#84C1FF",   
  "Unknown" = "lightgray")

ggraph(network, layout = "circle") + 
  geom_edge_arc(aes(edge_width = number,edge_color=lineage),curvature = 0.1,
                alpha = 0.4) +
  scale_edge_width(range = c(0.8, 2)) +
  geom_node_point(aes(fill=type,shape=shape), size = 7) + 
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = node_colors) + 
  scale_edge_color_manual(values = edge_colors) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
  theme_void()

#Fig.s3
#algae
data<-read_excel('Fig. S3.xlsx')
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(color="#4CAF50",size = 1, alpha = 0.45,shape=16,stroke=0.22)+
  labs(x="Abundance of Algae",y="Abundance of NCLDVs")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#4CAF50",method = 'lm',se=F,size=0.4,fullrange=T,fill="#4CAF50",alpha=0.1)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+
  ylim(0,8)+xlim(0,9)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),size=2.3,
               label.x = )  
ggsave("Fig.S3a.pdf",width=2,height=2.2,path="D:/")

#fungi
data<-read_excel('Fig. S3.xlsx')
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(color="#d3a4ff",size = 1, alpha = 0.45,shape=16,stroke=0.22)+
  labs(x="Abundance of Fungi",y="Abundance of NCLDVs")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#d3a4ff",method = 'lm',se=F,size=0.4,fullrange=T,fill="#d3a4ff",alpha=0.1)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),size=2.3,
           label.x = )  
ggsave("Fig.S3b.pdf",width=2,height=2.2,path="D:/")

#Animal
data<-read_excel('Fig. S3.xlsx')
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(color="#FF9224",size = 1, alpha = 0.45,shape=16,stroke=0.22)+
  labs(x="Abundance of Animal",y="Abundance of NCLDVs")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#FF9224",method = 'lm',se=F,size=0.4,fullrange=T,fill="#FF9224",alpha=0.1)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+ylim(0,7)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),size=2.3,
           label.x = )  
ggsave("Fig.S3C.pdf",width=2,height=2.2,path="D:/")

#Protozoa
data<-read_excel('Fig. S3.xlsx')
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(color="#FFD306",size = 1, alpha = 0.45,shape=16,stroke=0.22)+
  labs(x="Abundance of Protozoa",y="Abundance of NCLDVs")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="#FFD306",method = 'lm',se=F,size=0.4,fullrange=T,fill="#FFD306",alpha=0.1)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+ylim(0,6)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),size=2.3,
           label.x = )  

#Fig.s4
setwd('Fig. S4')
data<-read_excel('Fig. S8.xlsx')
ggplot(data, mapping = aes(x=Host,y=Virus))+
  geom_point(fill="lightgray",size = 1, alpha = 0.45,shape=21,stroke=0.22)+
  labs(x="Abundance of hosts",y="Abundance of jumbo phages")+
  theme(axis.line = element_line(color="black"))+
  geom_smooth(color="black",method = 'lm',se=T,size=0.4,fullrange=T,fill="black",alpha=0.1)+
  theme_test()+ theme(legend.position ="none")+
  theme( axis.title.x = element_text(size=8),
         axis.title.y = element_text(size=8),
         axis.text.x = element_text(hjust =0.5,size=7,colour = 'black'),
         axis.text.y=element_text(size=7,colour = 'black'),
         axis.ticks = element_line(size = 0.25),
         axis.ticks.length = unit(0.05, "cm"),
         panel.border = element_rect(size=0.45))+ylim(0,120)+xlim(0,500)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),size=2.3,
           label.x = )  
