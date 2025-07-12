library(gplots)  
library(pheatmap) 
library(readxl) 
library(vegan)
library(ggpubr)
library(ggpmisc)
library(ggplot2)
library(SoDA)
library(splines)
library(ape)
library(ggridges)
#Fig. 6a
#national_cluster
setwd("Fig. 6/")
data <- read.delim('OTU_national.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
data<-as.matrix(data)
drows<-vegdist(data,method="bray")
dcols<-vegdist(t(data),method="bray")
annotation_col<- read.delim("group_national.txt",header=T,sep="\t",row.names=1)                                                                                         
annotation_col <- as.data.frame(annotation_col)
ann_colors = list(Sample=c(DD="#6A5ACD",DY="#6495ED",QD="#1E90FF",LYG="#7CCD7C",YC="#B4EEB4",
                           NB="#EE9A00",WZ="#FFAF60",XM="#DDA0DD",
                           ST="#EEC900",ZH="#FFD700",BH="#EE7942",SY="#20B2AA"))
pheatmap(log2(data+1),scale = "none",clustering_method = "average",
         cluster_distance_cols = dcols,annotation_colors = ann_colors,annotation_col=annotation_col,
         clustering_distance_rows = drows,cluster_cols=T,cluster_row=F,treeheight_col =25,treeheight_row =25,
         show_colnames=F,show_rownames=F,border_color=NA,cellwidth = 2, 
         cellheight = 5,gaps_row = 47,file="national_cluster1.pdf")

#time_cluster
data <- read.delim('OTU_time.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
data<-as.matrix(data)
drows<-vegdist(data,method="bray")
dcols<-vegdist(t(data),method="bray")
annotation_col<- read.delim("group_time.txt",header=T,sep="\t",row.names=1)                                                                                         
annotation_col <- as.data.frame(annotation_col)
ann_colors = list(Sample=c(Aug="#D4EDFF",Oct="#8CC0FF",Dec="#5AA9FF",Feb="#2894FF",April="#0068C6",June="#003D79"))
pheatmap(log2(data+1),scale = "none",clustering_method = "average",
         cluster_distance_cols = dcols,annotation_colors = ann_colors,annotation_col=annotation_col,
         clustering_distance_rows = drows,cluster_cols=T,cluster_row=F,treeheight_col =25,treeheight_row =25,
         show_colnames=F,show_rownames=F,border_color=NA,cellwidth = 2, 
         cellheight = 5,gaps_row = 47,file="time_cluster.pdf")

#depth_cluster
data <- read.delim('OTU_depth.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
data<-as.matrix(data)
drows<-vegdist(data,method="bray")
dcols<-vegdist(t(data),method="bray")
annotation_col<- read.delim("group_depth.txt",header=T,sep="\t",row.names=1)                                                                                         
annotation_col <- as.data.frame(annotation_col)
ann_colors = list(Sample=c(Twenty="#E8F7F9",Forty="#C1EBF0",Sixty="#92D3D6",Eighty="#55B5B7",Hundred="#3F8D93"))
pheatmap(log2(data+1),scale = "none",clustering_method = "average",
         cluster_distance_cols = dcols,annotation_colors = ann_colors,annotation_col=annotation_col,
         clustering_distance_rows = drows,cluster_cols=T,cluster_row=F,treeheight_col =25,treeheight_row =25,
         show_colnames=F,show_rownames=F,border_color=NA,cellwidth = 2, 
         cellheight = 5,gaps_row = 47,file="depth_cluster.pdf")

#OTU_heatmap
data <- read.delim('OTU_all.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
data<-as.matrix(data)
drows<-vegdist(data,method="bray")
dcols<-vegdist(t(data),method="bray")
annotation_col<- read.delim("group_all.txt",header=T,sep="\t",row.names=1)                                                                                         
annotation_col <- as.data.frame(annotation_col)
ann_colors = list(Sample=c(DD="#6A5ACD",DY="#6495ED",QD="#1E90FF",LYG="#7CCD7C",YC="#B4EEB4",
                           NB="#EE9A00",WZ="#FFAF60",XM="#DDA0DD",
                           ST="#EEC900",ZH="#FFD700",BH="#EE7942",SY="#20B2AA",
                           Aug="#D4EDFF",Oct="#8CC0FF",Dec="#5AA9FF",Feb="#2894FF",April="#0068C6",June="#003D79",
                           Twenty="#E8F7F9",Forty="#C1EBF0",Sixty="#9ADFE7",Eighty="#73D3DE",Hundred="#7AC5CD"))
pheatmap(log2(data+1),scale = "none",clustering_method = "average",
         cluster_distance_cols = dcols,annotation_colors = ann_colors,annotation_col=annotation_col,
         clustering_distance_rows = drows,cluster_cols=F,cluster_row=F,treeheight_col =25,treeheight_row =25,
         show_colnames=F,show_rownames=F,border_color=NA,cellwidth = 2, 
         cellheight = 5,gaps_row = 71,gaps_col = c(93,129),file="Fig. S12a_pheatmap.pdf")

#alpha diversity
otu <- read.delim('phage_OTU_national.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu<-t(otu)
richness <- rowSums(otu > 0)
shannon_index <- diversity(otu, index = 'shannon', base = exp(1))

#alpha_time
data<-read_excel("Fig. 6.xlsx")
data$Time <- factor(data$Time,levels=c("Aug","Oct","Dec","Feb","April","June"))
ggplot(data, mapping = aes(x=Time,y=Richness))+
  geom_boxplot(aes(group=Time),color="black",width=0.4,outlier.color="white",size=0.5,
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Time,y=Richness),fill="gray",width=0.15,shape=21,size=1.4,alpha=0.7)+
  labs(x="Sampling time (month)",y="Alpha diversity (richness)")+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color="#6BB5FF",fill="#ACD6FF",
              alpha=0.2,aes(group=1),formula = y ~ bs(x, df = 3))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+ylim(5,40)+
  facet_wrap(~Tax)

#alpha_depth
data<-read_excel("Fig. 6.xlsx")
data$Depth <- factor(data$Depth,levels=c("20","40","60","80","100"))
ggplot(data, mapping = aes(x=Depth,y=Richness))+
  geom_boxplot(aes(group=Depth),color="black",width=0.4,outlier.color="white",size=0.5,
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Depth,y=Richness),fill="gray",width=0.15,shape=21,size=1.4,alpha=0.7)+
  labs(x="Sampling depth (cm)",y="Alpha diversity (richness)")+
  geom_smooth(method = 'lm',se=T,size=1,fullrange=T,color="#7AC5CD",fill="#9ADFE7",
              alpha=0.2,aes(group=1),formula = y ~ bs(x, df = 3))+
  theme_test()+ 
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)+ylim(5,40)

#alpha_national
data<-read_excel("Fig. 6.xlsx")
data$Site <- factor(data$Site,levels=c('DD','DY','QD','LYG','YC','NB','WZ','XM','ST','ZH','BH','SY'))
ggplot(data, mapping = aes(x=Site,y=Richness,fill=Site))+
  geom_boxplot(aes(group=Site),color="black",width=0.4,outlier.color="white",size=0.4,
               position = position_dodge(0.8))+
  geom_jitter(aes(x=Site,y=Richness),width=0.15,shape=21,size=1.4,alpha=0.7)+
  scale_fill_manual(breaks=c('DD','DY','QD','LYG','YC','NB','WZ','XM','ST','ZH','BH','SY'),
                    values=c("#6A5ACD","#6495ED","#1E90FF","#7CCD7C","#B4EEB4","#EE9A00","#FFAF60",
                             "#DDA0DD", "#EAC100","#FFD700","#EE7942","#20B2AA"))+
labs(x="",y="Alpha diversity (richness)")+
  theme_test()+ 
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)+coord_flip()

#Fig. 6b
#PCOA_time
otu <- read.delim('Phage_OTU_time.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- data.frame(t(otu))
otu2 <- log(otu1+1)
group <- read.delim('group_time.txt', sep = '\t', stringsAsFactors = FALSE)
dist <- vegdist(otu2, method = "bray")
pcoa_result <- pcoa(dist)
sample_site <- data.frame(PCoA1 = pcoa_result$vectors[, 1], 
                          PCoA2 = pcoa_result$vectors[, 2])

sample_site$names <- rownames(sample_site)
sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
variance_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)
PERMANOVA<-adonis2(dist ~ Sample, data = group, permutations = 999)

otu<-read.csv('PCOA.csv')
ggplot(data = otu, mapping = aes(PCoA1, PCoA2)) + 
  geom_point(aes(fill = Time,shape=Tax), size = 2, alpha = 0.9) + 
  scale_shape_manual(breaks = c('Nucleocytoviricota','Uroviricota'),values=c(21,24))+
  scale_fill_manual(breaks=c("Aug","Oct","Dec","Feb","April","June"),
                    values=c("#D4EDFF", "#8CC0FF", "#5AA9FF", "#2894FF", "#0068C6", "#003D79"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)+xlim(-0.5,0.5)+ylim(-0.4,0.4)


#PCOA_depth
otu <- read.delim('phage_OTU_depth.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- data.frame(t(otu))
otu2 <- log(otu1+1)
group <- read.delim('group_depth.txt', sep = '\t', stringsAsFactors = FALSE)
dist <- vegdist(otu2, method = "bray")
pcoa_result <- pcoa(dist)
sample_site <- data.frame(PCoA1 = pcoa_result$vectors[, 1], 
                          PCoA2 = pcoa_result$vectors[, 2])

sample_site$names <- rownames(sample_site)
sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
variance_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)
PERMANOVA<-adonis2(dist ~ Sample, data = group, permutations = 999)
write.table (sample_site,file ="PCOA_phage.csv", row.names = T, col.names =FALSE, quote =FALSE) 

otu<-read.csv('PCOA.csv')
ggplot(data = otu, mapping = aes(PCoA1, PCoA2)) + 
  geom_point(aes(fill = Depth,shape=Tax), size = 2, alpha = 0.9) + 
  scale_shape_manual(breaks = c('Nucleocytoviricota','Uroviricota'),values=c(21,24))+
  scale_fill_manual(breaks=c("Twenty","Forty","Sixty","Eighty","Hundred"),
                    values=c("#E8F7F9","#C1EBF0","#92D3D6","#55B5B7","#3F8D93"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)+ylim(-0.25,0.3)

#PCOA_national
otu <- read.delim('NCLDV_OTU_national_bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu1 <- data.frame(t(otu))
otu2 <- log(otu1+1)
group <- read.delim('group_national.txt', sep = '\t', stringsAsFactors = FALSE)
dist <- vegdist(otu2, method = "bray")
pcoa_result <- pcoa(dist)
sample_site <- data.frame(PCoA1 = pcoa_result$vectors[, 1], 
                          PCoA2 = pcoa_result$vectors[, 2])

sample_site$names <- rownames(sample_site)
sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
variance_explained <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)
PERMANOVA<-adonis2(dist ~ Sample, data = group, permutations = 999)

otu<-read.csv('PCOA.csv')
ggplot(data = otu, mapping = aes(PCoA1, PCoA2)) + 
  geom_point(aes(fill = Site,shape=Tax), size = 2, alpha = 0.9) + 
  scale_shape_manual(breaks = c('Nucleocytoviricota','Uroviricota'),values=c(21,24))+
  scale_fill_manual(breaks=c('DD','DY','QD','LYG','YC','NB','WZ','XM','ST','ZH','BH','SY'),
                    values=c("#6A5ACD","#6495ED","#1E90FF","#7CCD7C","#B4EEB4","#EE9A00","#FFAF60",
                             "#DDA0DD", "#EAC100","#FFD700","#EE7942","#20B2AA"))+  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
        axis.text.x = element_text(hjust =0.5,size=12,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+facet_wrap(~Tax)

#Fig. 6c
otu<-read.delim('NCLDV_OTU_national_bray.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
loca<-read_excel('DDR.xlsx')
samp.ck=match.name(rn.list=list(t(otu),loca))

distv<-vegdist(t(otu),method = 'bray')
distv_num<-1-as.numeric(distv)
distv_num<-log10(distv_num)
distv_num<-as.numeric(distv_num)
loca_geo<-geoXY(loca$Latitude,loca$Longitude)
rownames(loca_geo)<-rownames(loca)
dist_loca<-vegdist(loca_geo,method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)
dist<-data.frame(distv_num,dist_loca_num)
write.csv(dist,file="NCLDV_dist.csv",row.names = F)
#dist_df <- melt(as.matrix(distv))
#dist_df <- melt(as.matrix(dist_loca))

#distance_calculate_local
otu<-read.delim('Phage_local_OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
distv<-vegdist(t(otu),method = 'bray')
distv_num<-1-as.numeric(distv)
distv_num<-log10(distv_num)
distv_num<-as.numeric(distv_num)
loca_geo<-read.delim('DDR.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
dist_loca<-vegdist(loca_geo,method = 'euclidean')
dist_loca_log<-log10(dist_loca)
dist_loca_num<-as.numeric(dist_loca_log)
dist<-data.frame(distv_num,dist_loca_num)

#DDR_national
data<-read.csv('bray_national.csv')
data$Tax=factor(data$Tax, levels=c("Nucleocytoviricota","Uroviricota"))
data<-data[order(data$Tax, decreasing = TRUE), ]
ggplot(data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(color=Tax),size = 1, alpha = 0.5,shape=16)+
  theme(axis.line = element_line(color="black"))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+ 
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#7AC5CD"))+
  geom_smooth(aes(color=Tax,fill=Tax),method = 'lm',se=T,size=0.8,fullrange=F,alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+ylim(-2.3,0.3)+
  stat_poly_eq(aes(color=Tax,label = paste(..eq.label..,..p.value.., sep = "~`,`~")),size=5,
               formula = y ~ x,parse = TRUE)

#DDR_local
data<-read.csv('bray_local.csv')
data$Tax=factor(data$Tax, levels=c("Nucleocytoviricota","Uroviricota"))
data<-data[order(data$Tax, decreasing = TRUE), ]
ggplot(data,aes(x=dist_loca_num,y=distv_num))+
  geom_point(aes(color=Tax),size = 1, alpha = 0.5,shape=16)+
  theme(axis.line = element_line(color="black"))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+ 
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#7AC5CD"))+
  geom_smooth(aes(color=Tax,fill=Tax),method = 'lm',se=T,size=0.8,fullrange=F,alpha=0.2)+
  theme_test()+ theme(legend.position ="none")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+
  ylim(-0.5,0.07)+
  stat_poly_eq(aes(color=Tax,label = paste(..eq.label..,..p.value.., sep = "~`,`~")),size=5,
               formula = y ~ x,parse = TRUE)
#Fig. 6d
#Pi value_time
data<-read_excel("Fig. 6.xlsx")
data$Time <- factor(data$Time,levels=c("Aug","Oct","Dec","Feb","April","June"))
ggplot(data, aes(x = Pi, y = Time,fill=Tax)) +
  geom_density_ridges(alpha=0.8,scale=0.6,linewidth = 0.3)+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+xlim(-0.0004,0.0055)

data<-read_excel("Fig. 6.xlsx")
data$Time <- factor(data$Time,levels=c("Aug","Oct","Dec","Feb","April","June"))
ggplot(data,aes(x=Time,y=Pi,fill=Tax,group=Tax))+
  stat_summary(alpha=0.8,fun = mean,geom="bar",color="black",width=0.7,size=0.3,position = position_dodge(0.7))+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.2,size=0.3,position = position_dodge(0.7))+
  geom_jitter(aes(fill=Tax),shape=21,size=1,alpha=0.9,stroke = 0.2,position=position_jitterdodge(jitter.width = 0.25,dodge.width = 0.8))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_blank(),
        panel.border = element_rect(size=0.9),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+coord_flip()+
  scale_y_continuous(limits=c(0,0.0012),breaks = c(0,0.0006,0.0012))


#Pi value_depth
data<-read_excel("Fig. 6.xlsx")

data$Depth <- factor(data$Depth,levels=c("20","40","60","80","100"))
ggplot(data, aes(x = Pi, y = Depth,fill=Tax)) +
  geom_density_ridges(alpha=0.8,scale=0.5,linewidth = 0.3)+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")

data<-read_excel("Fig. 6.xlsx")
data$Depth <- factor(data$Depth,levels=c("20","40","60","80","100"))
ggplot(data,aes(x=Depth,y=Pi,fill=Tax,group=Tax))+
  stat_summary(alpha=0.8,fun = mean,geom="bar",color="black",width=0.6,size=0.3,position = position_dodge(0.6))+
  stat_summary(fun.data = mean_se,geom="errorbar",width=.2,size=0.3,position = position_dodge(0.7))+
  geom_jitter(aes(fill=Tax),shape=21,size=1,alpha=0.9,stroke = 0.2,position=position_jitterdodge(jitter.width = 0.25,dodge.width = 0.7))+
  labs(x="",y="")+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+
  scale_color_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                     values=c("#6BB5FF","#7AC5CD"))+
  theme(axis.line = element_line(color="black"))+
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_blank(),
        panel.border = element_rect(size=0.9),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=14),
        legend.position = "none")+coord_flip()+
  scale_y_continuous(limits=c(0,0.0014),breaks = c(0,0.0007,0.0014))

#Pi value_national
data<-read_excel("Fig. 6.xlsx")
data$Site <- factor(data$Site,levels=c('DD','DY','QD','LYG','YC','NB','WZ','XM','ST','ZH','BH','SY'))
ggplot(data, aes(x = Pi, y = Site,fill=Tax)) +
  geom_density_ridges(alpha=0.8,scale=0.9,linewidth = 0.3)+
  scale_fill_manual(breaks=c("Nucleocytoviricota","Uroviricota"),
                    values=c("#6BB5FF","#7AC5CD"))+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_test()+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(hjust =0.5,size=13,colour = 'black'),
        axis.text.y=element_text(size=12,colour = 'black'),
        panel.border = element_rect(size=0.8),
        axis.ticks = element_line(size = 0.4),
        axis.ticks.length = unit(0.08, "cm"),
        strip.text=element_text(size=15),
        legend.position = "none")+xlim(-0.0005,0.0075)
