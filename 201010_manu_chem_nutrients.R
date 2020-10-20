##load libraries
library(gplots)
library(ggplot2)
library(plyr)
library(matrixStats)
library(plotrix)
library(gtools)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(lemon)

chem<-read.delim("191115_waiwera_chem.txt",header = T,row.names = 1)
chem$Distance<-as.character(chem$Distance)

chem2<-chem %>% 
    melt( variable.name = ("variable"), value.name="Values") %>% 
  group_by(group,variable,place,Distance) %>%   
  summarise_each(funs(mean,sd,std.error))

map_dist<-read.delim("map_distance.txt",header = F)

q<-ggplot(filter(chem2,place=="Sed" & variable %in% c("TRP","TS","TN","TOC")), aes(x=group, y=mean, fill=variable)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,  position=position_dodge(.9))+facet_grid(variable~place, scales="free_y")+ theme_bw() + 
    scale_fill_manual(name="Nutrients", values=c("#F8766D", "#E18A00", "#BE9C00", "#8CAB00"))+
    scale_y_continuous(breaks =scales::pretty_breaks(n = 3))+
    theme(legend.position="bottom",legend.title = element_text(size=14),
          legend.text=element_text(size=12),strip.text = element_text(size=13),
          axis.text=element_text(size=11),axis.title = element_text(size=11),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Concentrations (g/m3)") +xlab("Distance")+
    scale_x_discrete(labels = lapply(map_dist[c(1:2,4:11,13:18),],round,1))+
    labs(fill="Nutrients")+guides(fill=guide_legend(ncol=5))

p<-ggplot(filter(chem2,!(variable %in% c("TRP","TS","TN","TOC"))), aes(x=group, y=mean, fill=variable)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,  position=position_dodge(.9))+facet_grid(variable~place, scales="free_y")+ theme_bw() + 
    scale_fill_manual(name="Nutrients", values=c("#24B700", "#00BE70", "#00C1AB", "#00BBDA", "#00ACFC", "#8B93FF", "#D575FE", "#F962DD", "#FF65AC"))+
    scale_y_continuous(breaks =scales::pretty_breaks(n = 3))+
    theme(legend.position="bottom",legend.title = element_text(size=14),
          legend.text=element_text(size=12),strip.text = element_text(size=13),
          axis.text=element_text(size=11),axis.title = element_text(size=11),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Concentrations (g/m3)") +xlab("Distance")+
    scale_x_discrete(labels = lapply(map_dist[c(1:2,4:11,13:18),],round,1))+
    labs(fill="Nutrients")+guides(fill=guide_legend(ncol=5))


pdf("chem_barplot.pdf",width=16,height=10)
grid_arrange_shared_legend(p,arrangeGrob(q,heights = c(0.64,0.36),widths=c(0.65)),ncol = 2,nrow = 1,widths=c(0.64,0.36))
dev.off()

#ZOOM IN AMMONIUM 
ggplot(NH4[9:16,], aes(x=group, y=mean))+ theme_bw() + 
    geom_bar(position=position_dodge(), stat="identity",fill="#2bb45b") +
    geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=.2,  position=position_dodge(.9))+
    theme(legend.position=c(0.9,0.15),legend.title = element_text(size=23),
	legend.text=element_text(size=22),strip.text = element_text(size=24),
	axis.text=element_text(size=19),axis.text.x =element_text(angle=90),axis.title = element_text(size=21)) +
	xlab("Distance")+ylab("Concentrations (g/m3)") 
