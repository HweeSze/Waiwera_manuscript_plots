##load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(hues)
library(lemon)
library(stringr)
library(gridExtra)
library(RColorBrewer)

#read files
map<-read.delim2("mappingcomb.txt")
taxa_wgs_all<-read.delim2("taxa_wgs_all.txt",sep = " ")
taxa_wts_all<-read.delim2("taxa_wts_all.txt",sep = " ")


#convert taxa to character
taxa_wgs_all$taxac<-as.character(taxa_wgs_all$taxac)

### preparation for lineplot

make_gradient <- function(deg = 45, n = 100, cols = blues9) {
    cols <- colorRampPalette(cols)(n + 1)
    rad <- deg / (180 / pi)
    mat <- matrix(
        data = rep(seq(0, 1, length.out = n) * cos(rad), n),
        byrow = TRUE,
        ncol = n
    ) +
        matrix(
            data = rep(seq(0, 1, length.out = n) * sin(rad), n),
            byrow = FALSE,
            ncol = n
        )
    mat <- mat - min(mat)
    mat <- mat / max(mat)
    mat <- 1 + mat * n
    mat <- matrix(data = cols[round(mat)], ncol = n)
    grid::rasterGrob(
        image = mat,
        width = unit(1, "npc"),
        height = unit(1, "npc"), 
        interpolate = TRUE
    )
}

line_wgs <- taxa_wgs_all %>% dplyr::select(X2,X3, f1:S9R3) %>%
    group_by(X2,X3) %>%
    summarise_at(vars(f1:S9R3), ~sum(.)) %>%
    gather(key = "sample", value = "values", f1:S9R3)

line_wts <- taxa_wts_all %>% dplyr::select(X2,X3, f1:S7R3) %>%
    group_by(X2,X3) %>%
    summarise_at(vars(f1:S7R3), ~sum(.)) %>%
    gather(key = "sample", value = "values", f1:S7R3)

line_wgs$X3 <- sapply(line_wgs$X3, function(x) paste( x,' WGS'))
line_wts$X3 <- sapply(line_wts$X3, function(x) paste( x,' WTS'))

#combine wgs wts
line_wgs_wts<-rbind(line_wgs,line_wts)

#log rpm and TPM counts
filter_line_wgs_wts<-left_join(line_wgs_wts,map,by=c("sample"="X")) %>% 
  filter(values !=0) %>% 
  mutate(values = values+1) %>% 
  mutate(log=log10(values))
filter_line_wgs_wts$salinity<-as.numeric(as.character(filter_line_wgs_wts$salinity))

filter_line_wgs_wts <-as.data.frame(filter_line_wgs_wts)  %>%
  mutate(X3= case_when(X3 =="Rubisco WGS" ~ "Calvin cycle WGS", TRUE ~as.character(X3)))
filter_line_wgs_wts <-as.data.frame(filter_line_wgs_wts)  %>%
  mutate(X3= case_when(X3 =="Rubisco WTS" ~ "Calvin cycle WTS", TRUE ~as.character(X3)))

linedata<-separate(data = filter_line_wgs_wts, col = X3, into = c("Specific.pathway", "Data"), sep = " W")
filter_line_wgs_wts$salinity<-as.numeric(as.character(filter_line_wgs_wts$salinity))

linedata$Specific.pathway = str_replace(linedata$Specific.pathway,"_"," ")
linedata$Data = str_replace(linedata$Data,"GS","WGS")
linedata$Data = str_replace(linedata$Data,"TS","WTS")
linedata$Specific.pathway = str_replace(linedata$Specific.pathway,"Sodium−translocating dehydrogenases_nqrF","Na+ translocating dehydrogenase nqrF")
linedata$Specific.pathway = str_replace(linedata$Specific.pathway,"Potassium","K+")


linedata$place = factor(linedata$place, levels=c('Water column','Sediments'))

##lineplot
g <- make_gradient(
  deg = 180, n = 500, cols = alpha(rev(brewer.pal(9, "Blues")),0.1))

d1<-ggplot(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|PhoU phosphate regulon|Specific Transport",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|PhoU phosphate regulon|Specific Transport",linedata$Specific.pathway),],Data=="WGS"),
  color="black")+facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|PhoU phosphate regulon|Specific Transport",linedata$Specific.pathway),],Data=="WTS"),
  color="darkgrey")+
  scale_y_continuous(name = "Log10 TPM") +theme_bw() + 
  ggtitle("Phosphorus")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

d2<-ggplot(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|Thiosulfate",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|Thiosulfate",linedata$Specific.pathway),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|Thiosulfate",linedata$Specific.pathway),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 TPM") +theme_bw() + 
  ggtitle("Sulfur")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d3<-ggplot(linedata[grep("N fixation|Cyanophycin|Ammonia|Complete|DNRA|nitrate|nitric|nitrite|nitrous",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("N fixation|Cyanophycin|mmonia oxidation|DNRA|nitrate|nitric|nitrite|nitrous",linedata$Specific.pathway),],Data=="WGS"),
  color="black")+facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("N fixation|Cyanophycin|mmonia oxidation|DNRA|nitrate|nitric|nitrite|nitrous",linedata$Specific.pathway),],Data=="WTS"),
  color="darkgrey")+
  scale_y_continuous(name = "Log10 TPM") +theme_bw() + ggtitle("Nitrogen")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d4<-ggplot(linedata[grep("glycine|K+|Na+",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("glycine|K+|Na+",linedata$Specific.pathway),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("glycine|K+|Na+",linedata$Specific.pathway),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 TPM") +theme_bw() + 
  ggtitle("Osmoregulation")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d5<-ggplot(linedata[grep("Photosynthesis",linedata$X2), ],col=pal, aes(salinity, log)) +
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Photosynthesis",linedata$X2),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Photosynthesis",linedata$X2),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 TPM") +theme_bw() + 
  ggtitle("Photosynthesis")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),legend.position = "none")

d6<-ggplot(linedata[grep("Carbon fixation",linedata$X2), ],col=pal, aes(salinity, log)) +
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+  
  geom_line(data=filter(linedata[grep("Carbon fixation",linedata$X2),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15),scales="free")+
  geom_line(data=filter(linedata[grep("Carbon fixation",linedata$X2),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 TPM") +theme_bw() + 
  ggtitle("Carbon fixation")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")


pdf("wgs_wts_lineplot.pdf",width=20,height=18)
grid_arrange_shared_legend(d4,arrangeGrob(d6,d1,ncol=1, nrow=2,heights = c(0.35,0.65)),arrangeGrob(d5,d2, nrow=2, 
	ncol=1,heights = c(0.35,0.65)),arrangeGrob(d3, nrow=1, ncol=1,heights = c(1.0)),nrow = 1,ncol = 4)
dev.off()



#### preparation for individual function-taxa barplot
wgs.df2 <- select(taxa_wgs_all1,X2,X3,X4, taxac, 13:30) %>%
    group_by(X2,X3,X4, taxac) %>%
    summarise_at(vars(f1:'9'), ~sum(.)) %>%
    gather(key = "sample", value = "values",f1:'9')
wgs.df2=wgs.df2 %>%mutate(place = case_when(sample %in% c("f1","f2","f3","f4","f5","f6","f7","f8","f9") ~ as.character("Water column"),TRUE ~ as.character("Sediment")))
wgs.df2$sample<-str_replace(wgs.df2$sample,"f","")


wts.df2 <- select(taxa_wts_all1, X2,X3,X4, taxac, f1:"7") %>%
    group_by(X2,X3,X4, taxac) %>%
    summarise_at(vars(f1:"7"), ~sum(as.numeric(.))) %>%
    gather(key = "sample", value = "values", f1:"7")
wts.df2=wts.df2 %>%mutate(place = case_when(sample %in% c("f1","f2","f3","f4","f5","f6","f7","f8","f9") ~ as.character("Water column"),TRUE ~ as.character("Sediment")))
wts.df2$sample<-str_replace(wts.df2$sample,"f","")


##all
colnames(wgs.df2)=c("Function","Specific.pathway","Gene", "taxac",  "sample" ,"values" ,"place")
colnames(wts.df2)=c("Function","Specific.pathway","Gene", "taxac",  "sample" ,"values" ,"place")

wgs.df2$Specific.pathway <- sapply(wgs.df2$Specific.pathway, function(x) paste( x,' WGS'))
wts.df2$Specific.pathway <- sapply(wts.df2$Specific.pathway, function(x) paste( x,' WTS'))


wgs_wts<-rbind(wgs.df2,wts.df2)
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"Rubisco","Calvin cycle")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"wood_lungdahl_pathway","wood lungdahl pathway")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"dehydrogenases_","dehydrogenases ")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"oxide_reduction","oxide reduction")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"Sodium","Na+")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"_"," ")

wgs_wts$place = factor(wgs_wts$place, levels=c('Water column','Sediment'))
wgs_wts$taxac<-str_replace(wgs_wts$taxac,"__$","") 
wgs_wts$taxac<-str_replace(wgs_wts$taxac,"_$","") 
wgs_wts$taxac<-str_replace(wgs_wts$taxac,"unknown$","Unknown") 
wgs_wts$Specific.pathway=str_replace(wgs_wts$Specific.pathway,"^nitr","Nitr")
wgs_wts$Specific.pathway=str_replace(wgs_wts$Specific.pathway,"^Phosphate Regulon","PhoR PhoB Phosphate Regulon")
wgs_wts$Specific.pathway=str_replace(wgs_wts$Specific.pathway,"^PhoR PhoB Inhibitor Protein","PhoU Phosphate Regulon")


### Individual barplot
scaleFUN <- function(x) sprintf("%.1f", x)
set_panel_size <- function(p = NULL, g = ggplotGrob(p),
                           file = NULL, margin = unit(1, "mm"), width = unit(4,
                                                                             "cm"), height = unit(4, "cm")) {
    
    panels <- grep("panel", g$layout$name)
    panel_index_w <- unique(g$layout$l[panels])
    panel_index_h <- unique(g$layout$t[panels])
    nw <- length(panel_index_w)
    nh <- length(panel_index_h)
    
    if (getRversion() < "3.3.0") {
        
        # the following conversion is necessary because
        # there is no `[<-`.unit method so promoting to
        # unit.list allows standard list indexing
        g$widths <- grid:::unit.list(g$widths)
        g$heights <- grid:::unit.list(g$heights)
        
        g$widths[panel_index_w] <- rep(list(width),
                                       nw)
        g$heights[panel_index_h] <- rep(list(height),
                                        nh)
        
    } else {
        
        g$widths[panel_index_w] <- rep(width, nw)
        g$heights[panel_index_h] <- rep(height, nh)
        
    }
    
    if (!is.null(file))
        ggsave(file, g, width = convertWidth(sum(g$widths) +
                                                 margin, unitTo = "in", valueOnly = TRUE),
               height = convertHeight(sum(g$heights) +
                                          margin, unitTo = "in", valueOnly = TRUE))
    
    g
}


pal<-iwanthue(118)
a2<-ggplot(wgs_wts[grep("Osmoregulation",wgs_wts$Function), ],col=pal, aes(sample, values, fill = taxac)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15,color ="black"),axis.title=element_text(size=18,color ="black"),
          legend.text=element_text(size=16,color ="black"),axis.text.y=element_text(size=15,color ="black"),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
          strip.text.x = element_text(size = 16, color ="black" ),strip.text.y = element_text(size = 16, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Osmoregulation", x = "Sites", y = "TPM")
p2=set_panel_size(a2,width = unit(5.5, "cm"),height = unit(4.5, "cm"))


v=wgs_wts[grep("Sulfur",wgs_wts$Function), ]
pal<-iwanthue(88)
a3<-ggplot(v[grep("dsr",v$Gene), ],col=pal, aes(sample, values, fill = taxac)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Gene+Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1, color ="black",vjust=0.5,size =15),axis.title=element_text(size=18, color ="black"),
        legend.text=element_text(size=16, color ="black"),axis.text.y=element_text(size=15, color ="black"),plot.title = element_text(face = 'bold', color ="black", size = 20,hjust=0.5),
        strip.text.x = element_text(size = 16, color ="black" ),strip.text.y = element_text(size = 16, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfite reduction and sulfur oxidation", x = "Sites", y = "TPM")
p3=set_panel_size(a3,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

a1<-ggplot(wgs_wts[grep("Photosynthesis|Carbon fixation",wgs_wts$Function), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( Function+Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
          legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
          strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Photosynthesis", x = "Sites", y = "TPM")
p1=set_panel_size(a1,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

v=wgs_wts[grep("Nitrogen",wgs_wts$Function), ]
pal<-iwanthue(93)
a4<-ggplot(v[grepl("Ammonia|Complete|Nitrite oxidation|Nitrite reduction|DNRA|Nitrate reduction|Nitric oxide|Nitrous oxide|Nitrite reduction",v$Specific.pathway), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( Function+Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=10))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Nitrification", x = "Sites", y = "TPM")
p4=set_panel_size(a4,width = unit(5.5, "cm"),height = unit(5.5, "cm"))


pal<-iwanthue(150)
a6<-ggplot(wgs_wts[grep("Cyanophycin|N fixation|Phosphatase|Phosphate Inorganic Transporter|regulon|Phosphate Regulon|Specific Transport",wgs_wts$Specific.pathway), ],col=pal, 
    aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( Function+Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=3))+ labs(title = "Phosphorus", x = "Sites", y = "TPM")
p6=set_panel_size(a6,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

pal<-iwanthue(88)
a7<-ggplot(wgs_wts[grep("Sulfate reduction|Sulfite reduction|Sulfide oxidation|Sulfur oxidation|Thiosulfate",wgs_wts$Specific.pathway), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( Function+Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
          legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
          strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=3))+ labs(title = "Sulfate and sulfite reduction", x = "Sites", y = "TPM")
p7=set_panel_size(a7,width = unit(5.5, "cm"),height = unit(5.5, "cm"))


pdf("wgs_wts.pdf",width=25,height=50)
marrangeGrob(grobs=list(p2,p1,p4,p6,p7,p3),  nrow = 1,ncol = 1)
dev.off()



##CORRELATIONS BETWEEN THIOSULFATE AND SULFATE CONCENTRATION
thiosulfate_oxidation=taxa_wts_all %>% filter(X3=="Thiosulfate_oxidation") %>%select(f1:S7R3) %>% colSums()
cor(as.numeric(as.character(map[1:9,17])),thiosulfate_oxidation[1:9],method="spearman")
#0.7615129

##CORRELATIONS BETWEEN PHOSPHATE UPTAKE AND PHOSPHATE CONCENTRATION
PHOSPHATE=taxa_wts_all[grep("Phosphate Specific Transport System|Phosphate Inorganic Transporter|Phosphate Regulon|regulon",taxa_wts_all$X3), ] %>%select(f1:S7R3) %>% colSums()
cor(as.numeric(as.character(map[1:9,17])),PHOSPHATE[1:9],method="spearman")
#-0.4937282


# Calculate nrfA and nosZ gene expression
taxa_wts_all1 %>% dplyr::filter(X4=="nrfA") %>%
    group_by(X4) %>%
    summarise_at(vars(f1:'7'), ~sum(.)) 

taxa_wts_all1 %>% dplyr::filter(X4=="nosZ") %>%
    group_by(X4) %>%
    summarise_at(vars(f1:'7'), ~sum(.))     



dev.new()
###use this
pal <- c("#b4ecf9", "#1e79da", "#0a0c93")
n1<-ggplot(data = mds.plot_wgs[["points"]],
mapping = aes(x = NMDS1, y = NMDS2)) +
geom_point(mapping = aes(shape = place, colour = salinity),
size = 4,stroke=2) + theme_bw()+
scale_colour_gradient(low = "#b4ecf9", high = "#0a0c93") +
scale_shape_manual(values=c(1,19))+
labs(shape = "Type", colour = "Salinity")+
ggtitle("nMDS WGS") +
geom_segment(data = test1,
mapping = aes(x = 0, y = 0,
xend = NMDS1* 0.55, yend = NMDS2 * 0.55),
arrow = arrow(length = unit(0.03, "npc"))) +
new_scale("colour") +
geom_text_repel(data = test1,size=5,segment.color = 'transparent', fontface = "bold",
mapping = aes(x = NMDS1 * 0.6, y = NMDS2 * 0.6,
label = rowname,colour = pathway)) +
scale_colour_manual(values = c("orange","#cc99cc", "blue", "#aa4000","purple"))+
labs(colour = "Pathway")+
theme(plot.title = element_text(hjust=0.5,size=14),axis.text=element_text(color = "black",
size=12),axis.title=element_text(color = "black",size=13),legend.text=element_text(size=12),
panel.grid.major = element_blank(), panel.grid.minor = element_blank())+scale_x_continuous(limits=c(-0.8,0.5))

n2<-ggplot(data = mds.plot_wts[["points"]],
mapping = aes(x = NMDS1, y = NMDS2)) +
geom_point(mapping = aes(shape = place, colour = salinity),
size = 4,stroke=2) + theme_bw()+
scale_colour_gradient(low = "#b4ecf9", high = "#0a0c93") +
scale_shape_manual(values=c(1,19))+
ggtitle("nMDS WTS") +
geom_segment(data = test2,
mapping = aes(x = 0, y = 0,
xend = NMDS1* 0.55, yend = NMDS2 * 0.55),
arrow = arrow(length = unit(0.03, "npc"))) +
new_scale("colour") +
geom_text_repel(data = test2,size=5,segment.color = 'transparent',fontface = "bold",
mapping = aes(x = NMDS1 * 0.6, y = NMDS2 * 0.6,
label = rowname,colour = pathway)) +
scale_colour_manual(values = c("orange","#cc99cc", "blue", "#aa4000","green","purple"))+
labs(shape = "Type", colour = "Salinity")+
theme(plot.title = element_text(hjust=0.5,
size=14),axis.text=element_text(color = "black",size=12),axis.title=element_text(color = "black",size=13),
legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
grid_arrange_shared_legend(n1,n2, nrow = 1)

wgs_all<-wgs_function[1:86,c(1,2,4:39)] %>%
group_by(pathway,func) %>%
summarise_at(vars(f1:S9R3), ~sum(.))%>%
column_to_rownames('func') %>%
dplyr::select(-c('pathway'))
wts_all<-wts_function[1:83,c(1,2,4:33)] %>%
group_by(pathway,func) %>%
summarise_at(vars(f1:S7R3), ~sum(.))%>%
column_to_rownames('func') %>%
dplyr::select(-c('pathway'))