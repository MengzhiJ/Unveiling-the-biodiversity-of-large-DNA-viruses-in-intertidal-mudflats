#load all need packages
library(ggplot2)
library(dplyr)
library(readxl)
library(ggpmisc)
library(reshape2)
library(ggforce)
library(pheatmap)
library(gplots)
library(gggenes)

setwd('Fig. 5/')
#Fig.5a
data<-read.csv('KEGG.csv',row.names = 1)
pheatmap(data,scale="none",
         cluster_cols=F,cluster_row=F,
         colorpanel(70,low="white",high="#A6D0F7"),
         gaps_row = c(70),legend=F,show_colnames=F,show_rownames=T,fontsize= 11.5,
         border_color = "#E0E0E0", cellwidth = 13, cellheight = 13,filename = "heatmap.pdf")

df <- as.data.frame(data)
df$Genome <- rownames(df)
df_long <- melt(
  df, 
  id.vars = "Genome", 
  variable.name = "Pathway", 
  value.name = "Presence"
)

data<-read.csv('KEGG_long.csv')
data$Genome <- factor(data$Genome, levels = unique(data$Genome))
ggplot(data)+
  geom_point(aes(size=Count,fill=Class,x=Class,y=Genome),shape=21,alpha=1,stroke = 0.2,color="#3C3C3C")+
  scale_size(range = c(0.5, 4))+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(color="black"))+
  theme_test()+theme(panel.border = element_rect(size=0.4))+
  theme(legend.position='none')+
  theme(axis.text.x = element_blank(),  
    axis.text.y = element_blank(),  
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank())+
  scale_fill_manual(breaks=c("Amino acid metabolism","Biosynthesis of other secondary metabolites","Carbohydrate metabolism",
                             "Energy metabolism","Glycan biosynthesis and metabolism","Lipid metabolism",
                             "Metabolism of cofactors and vitamins","Metabolism of other amino acids","Metabolism of terpenoids and polyketides",
                             "Nucleotide metabolism","Xenobiotics biodegradation and metabolism"),
                    values=c("#B9B9EF","#FFA6A6","#A6D0F7","#abefac","#fcf8b1","#ffbc9f",
                             "#8dd3c7","#f4d7ff","#ffe9b8","#add8e6","#f2f2f2"))

#Fig.5b
data<-read_excel('Fig. 5.xlsx')
data$genome <- factor(data$genome, levels = unique(data$genome))
ggplot(data, aes(xmin = start, xmax = end, y = genome,fill=Type,forward = (strand == 1)))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1.5, "mm"))+
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  scale_fill_manual(values=c("#7AC5CD","lightgray"))+
  theme_genes()+
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    legend.position = "none") 

#Fig.5c
data<-read_excel('Fig. 5.xlsx')
data$Gene <- factor(data$Gene, levels = unique(data$Gene))
ggplot(data)+
  geom_point(aes(size=Number,fill=Type,x=Gene,y=Type),shape=21,alpha=1,stroke = 0.2,color="#3C3C3C")+
  scale_size(range = c(1, 4))+
  theme_test()+
  theme(legend.position='none')+
  theme(axis.text.x = element_text(angle = 90,vjust =0,size=5,colour = 'black'),
        axis.text.y=element_text(size=5,colour = 'black'),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0.04, "cm"),
        panel.border = element_rect(size=0.4))+
  scale_fill_manual(values = c("#6BB5FF","#7AC5CD"))+labs(x="",y="")


#Fig.5e
#boxplot
data<-read_excel("Fig. 5.xlsx")
ggplot(data,aes(x=Type, y=Pi))+
  geom_boxplot(aes(fill=Type),width=0.28,size=0.2,alpha=0.9,outlier.shape = NA,outlier.size = NA)+
  geom_jitter(aes(x=Type, y=Pi,fill=Type),shape=21,stroke=0.15,size=0.6,
              position=position_jitterdodge(jitter.width = 0.3),alpha=0.9)+
  scale_fill_manual(values = c("#6BB5FF","#7AC5CD"))+
  scale_color_manual(values = c("#6BB5FF","#7AC5CD"))+
  theme_test()+ 
  theme(axis.text.x = element_text(size=6,colour = 'black'),
        axis.text.y=element_text(size=6,colour = 'black'),
        axis.ticks.length = unit(0.06, "cm"),
        axis.ticks = element_line(size = 0.25),
        panel.border = element_rect(size=0.4))+
   theme(legend.position="none")+labs(x="",y="")

#pie
sales <- c(1.05,98.95)
names<-c("Positive","Purify")
share<-sales/sum(sales)*100 
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),stat="pie",alpha=0.9,color = NA)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Positive","Purify"),values=c("#6BB5FF","#ECF5FF"))

data <- data.frame(variable = c(3,31,13,19,20), 
                   group = paste0("a", 1:5))
ggplot(data, aes(x = "", y = variable, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white", size =1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("a1", "a2", "a3", "a4", "a5"),
                    values = c("#B9B9EF","#8dd3c7","#ffe9b8", "#A6D0F7", "#F0F0F0"))

data <- data.frame(variable = c(19, 6, 3,34, 10, 43, 3, 5), 
                   group = paste0("a", 1:8))
ggplot(data, aes(x = "", y = variable, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = "white", size =1) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("a1", "a2", "a3", "a4", "a5", "a6", "a7","a8"),
                    values = c("#8dd3c7", "#ffe9b8","#ffbc9f", 
                               "#A6D0F7","#B9B9EF", 
                               "#F0F0F0", "#FFA6A6","#f4d7ff"))

setwd('Fig. S10/')
#Fig.S10
#pie
sales <- c(74,26)
names<-c("Known","Unknown")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),stat="pie",alpha=0.8)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Known","Unknown"),values=c("lightgray","black"))

sales <- c(41.9,58.1)
names<-c("Known","Unknown")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),stat="pie",alpha=0.6)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Known","Unknown"),values=c("#7AC5CD","white"))

#bar
data<-read_excel('Fig. S10.xlsx')
data_m<-melt(data) 
ggplot(data_m)+
  geom_bar(stat="identity",width=0.75, position="dodge", aes(x=COG, y=value,fill=class),alpha=0.8)+
  scale_fill_manual(values = c("#6BB5FF","#7AC5CD"))+
  theme(panel.border = element_blank())+
  labs(x="",y="Relative proportion (%)")+
  theme(axis.line = element_line(color="black"))+coord_flip()+
  theme_test()+theme( axis.title.x = element_text(size=5.5),
                      axis.title.y = element_text(size=5.5),
                      axis.text.x = element_text(hjust =0.5,size=5,colour = 'black'),
                      axis.text.y=element_text(size=5,colour = 'black'),
                      axis.ticks = element_line(size = 0.25),
                      axis.ticks.length = unit(0.04, "cm"),
                      strip.text=element_text(size=5.5, margin = margin(2, 2, 2, 2)),
                      strip.background = element_rect(color = "black", size = 0.35),
                      panel.border = element_rect(size=0.4))+
  theme(legend.position='none')+ylim(0,25)+facet_wrap(~class)    

setwd('Fig. S11/')
#Fig.S11
sales <- c(72.7,27.3)
names<-c("Known","Unknown")
share<-sales/sum(sales)*100
data <- data.frame(
  sales,share,names)
ggplot()+
  geom_arc_bar(data=data,aes(x0 = 0, y0 = 0, r0 = 0, r = 1,amount=sales,explode=c(0,0),fill=names),stat="pie",alpha=0.8)+
  coord_fixed()+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("Known","Unknown"),values=c("lightgray","black"))


data<-data.frame(variable=c(561,1021,618,174,239,146,120,72,3), group = paste0("a", 1:9))
ggplot(data, aes(x = 3, y = variable, fill = group))+ geom_col() +
  coord_polar(theta = "y") +xlim(c(0.5, 4.5))+theme_void()+theme(legend.position = "none")+
  scale_fill_manual(breaks=c("a1","a2","a3","a4","a5","a6","a7","a8","a9"),
                    values=c("#7AC5CD","#00A0FF","#A0E632","#D8B0FF","#A67A13","#8E9CFF",
                             "#ff9d6f","#E6AF2D","#EDE76D"))

data<-read_excel('Fig. S11.xlsx')
data$Class <- factor(data$Class, levels = unique(data$Class))
ggplot(data)+
  geom_bar(stat="identity",width=0.75, position="dodge", aes(x=Class, y=Proportion,fill=Type),alpha=0.8,color="black",size=0.3)+
  scale_fill_manual(values = c("#6BB5FF","#7AC5CD"))+
  theme(panel.border = element_blank())+
  theme(axis.line = element_line(color="black"))+
  theme_test()+theme(panel.border = element_rect(size=0.7))+
  theme(legend.position='none')+
  theme(
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()  
  ) +ylim(0,100)
