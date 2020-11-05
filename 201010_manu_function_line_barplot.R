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
map<-read.delim2("../mappingcomb.txt")
taxa_wgs_all<-read.delim2("../taxa_wgs_all.txt",sep = "\t")
taxa_wts_all<-read.delim2("../taxa_wts_all.txt",sep = "\t")

#Convert to numeric
taxa_wgs_all[,12:47] <- lapply(taxa_wgs_all[,12:47], as.character)
taxa_wts_all[,12:41] <- lapply(taxa_wts_all[,12:41], as.character)

taxa_wgs_all[,12:47] <- lapply(taxa_wgs_all[,12:47], as.numeric)
taxa_wts_all[,12:41] <- lapply(taxa_wts_all[,12:41], as.numeric)

#convert taxa to character
taxa_wgs_all$Class<-as.character(taxa_wgs_all$Class)

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

line_wgs <- taxa_wgs_all %>% dplyr::select(General.function,Specific.pathway, f1:S9R3) %>%
    group_by(General.function,Specific.pathway) %>%
    summarise_at(vars(f1:S9R3), ~sum(.)) %>%
    gather(key = "sample", value = "values", f1:S9R3)

line_wts <- taxa_wgs_all %>% dplyr::select(General.function,Specific.pathway, f1:S9R3) %>%
    group_by(General.function,Specific.pathway) %>%
    summarise_at(vars(f1:S7R3), ~sum(.)) %>%
    gather(key = "sample", value = "values", f1:S7R3)

line_wgs$Specific.pathway <- sapply(line_wgs$Specific.pathway, function(x) paste( x,' WGS'))
line_wts$Specific.pathway <- sapply(line_wts$Specific.pathway, function(x) paste( x,' WTS'))

#combine wgs wts
line_wgs_wts<-rbind(line_wgs,line_wts)

#log rpm and rpkm counts
filter_line_wgs_wts<-left_join(line_wgs_wts,map,by=c("sample"="X")) %>% 
  filter(values !=0) %>% 
  mutate(values = values+1) %>% 
  mutate(log=log10(values))
filter_line_wgs_wts<-filter_line_wgs_wts%>%
  mutate(place = case_when(place=="Filter" ~ as.character("Water Column"),TRUE ~ as.character("Sediments")))
filter_line_wgs_wts$salinity<-as.numeric(as.character(filter_line_wgs_wts$salinity))

filter_line_wgs_wts <-as.data.frame(filter_line_wgs_wts)  %>%
  mutate(Specific.pathway= case_when(Specific.pathway =="Rubisco WGS" ~ "Calvin cycle WGS", TRUE ~as.character(Specific.pathway)))
filter_line_wgs_wts <-as.data.frame(filter_line_wgs_wts)  %>%
  mutate(Specific.pathway= case_when(Specific.pathway =="Rubisco WTS" ~ "Calvin cycle WTS", TRUE ~as.character(Specific.pathway)))

linedata<-separate(data = filter_line_wgs_wts, col = Specific.pathway, into = c("Specific.pathway", "Data"), sep = " W")
filter_line_wgs_wts$salinity<-as.numeric(as.character(filter_line_wgs_wts$salinity))

linedata$Specific.pathway = str_replace(linedata$Specific.pathway,"_"," ")
linedata$Data = str_replace(linedata$Data,"GS","WGS")
linedata$Data = str_replace(linedata$Data,"TS","WTS")
linedata$Specific.pathway = str_replace(linedata$Specific.pathway,"Sodium−translocating dehydrogenases_nqrF","Na+ translocating dehydrogenase nqrF")
linedata$Specific.pathway = str_replace(linedata$Specific.pathway,"Potassium","K+")

typecol<- c("white", "gray88")

linedata$place = factor(linedata$place, levels=c('Water Column','Sediments'))

##lineplot
g <- make_gradient(
  deg = 180, n = 500, cols = alpha(rev(brewer.pal(9, "Blues")),0.1))

d1<-ggplot(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",linedata$Specific.pathway),],Data=="WGS"),
  color="black")+facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",linedata$Specific.pathway),],Data=="WTS"),
  color="darkgrey")+
  scale_y_continuous(name = "Log10 RPKM") +theme_bw() + 
  ggtitle("Phosphorus")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

d2<-ggplot(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|ThioSulfate",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|ThioSulfate",linedata$Specific.pathway),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|ThioSulfate",linedata$Specific.pathway),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 RPKM") +theme_bw() + 
  ggtitle("Sulfur")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d3<-ggplot(linedata[grep("N fixation|Cyanophycin|Ammonia|Complete|DNRA|nitrate|nitric|nitrite|nitrous",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("N fixation|Cyanophycin|mmonia oxidation|DNRA|nitrate|nitric|nitrite|nitrous",linedata$Specific.pathway),],Data=="WGS"),
  color="black")+facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("N fixation|Cyanophycin|mmonia oxidation|DNRA|nitrate|nitric|nitrite|nitrous",linedata$Specific.pathway),],Data=="WTS"),
  color="darkgrey")+
  scale_y_continuous(name = "Log10 RPKM") +theme_bw() + ggtitle("Nitrogen")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d4<-ggplot(linedata[grep("glycine|K+|Na+",linedata$Specific.pathway), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("glycine|K+|Na+",linedata$Specific.pathway),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("glycine|K+|Na+",linedata$Specific.pathway),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 RPKM") +theme_bw() + 
  ggtitle("Osmoregulation")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d5<-ggplot(linedata[grep("Photosynthesis",linedata$General.function), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Photosynthesis",linedata$General.function),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Photosynthesis",linedata$General.function),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 RPKM") +theme_bw() + 
  ggtitle("Photosynthesis")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),legend.position = "none")

d6<-ggplot(linedata[grep("Carbon fixation",linedata$General.function), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=Data),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+  
  geom_line(data=filter(linedata[grep("Carbon fixation",linedata$General.function),],Data=="WGS"),color="black")+
  facet_rep_grid(Specific.pathway~place,labeller = label_wrap_gen(width=15),scales="free")+
  geom_line(data=filter(linedata[grep("Carbon fixation",linedata$General.function),],Data=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log10 RPKM") +theme_bw() + 
  ggtitle("Carbon fixation")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")


pdf("wgs_wts_lineplot.pdf",width=15,height=12)
grid_arrange_shared_legend(d4,arrangeGrob(d6,d1,ncol=1, nrow=2,heights = c(0.35,0.65)),arrangeGrob(d5,d2, nrow=2, 
	ncol=1,heights = c(0.35,0.65)),arrangeGrob(d3, nrow=1, ncol=1,heights = c(1.0)),nrow = 1,ncol = 4)
dev.off()



#### preparation for individual function-taxa barplot

taxa_wgs_all1=cbind((taxa_wgs_all[,1:20]),
                    (t(apply(taxa_wgs_all[,21:47], 1, 
                    tapply, gl(9, 3), mean, na.rm = TRUE))),
                    (taxa_wgs_all[,48:70])) %>% 
                    as.data.frame()
taxa_wts_all1=cbind((taxa_wts_all[,1:20]),
                    t(apply(taxa_wts_all[,21:41], 1, 
                    tapply, gl(7, 3), mean, na.rm = TRUE)),
                    (taxa_wts_all[,42:65]))  %>% 
                    as.data.frame()

wgs.df2 <- select(taxa_wgs_all1,General.function,Specific.pathway,Gene, Class, 12:29) %>%
    group_by(General.function,Specific.pathway,Gene, Class) %>%
    summarise_at(vars(f1:'9'), ~sum(.)) %>%
    gather(key = "sample", value = "values",f1:'9')
wgs.df2=wgs.df2 %>%mutate(place = case_when(sample %in% c("f1","f2","f3","f4","f5","f6","f7","f8","f9") ~ as.character("Water column"),TRUE ~ as.character("Sediments")))
wgs.df2$sample<-str_replace(wgs.df2$sample,"f","")


wts.df2 <- select(taxa_wts_all1, General.function,Specific.pathway,Gene, Class, f1:"7") %>%
    group_by(General.function,Specific.pathway,Gene, Class) %>%
    summarise_at(vars(f1:"7"), ~sum(as.numeric(.))) %>%
    gather(key = "sample", value = "values", f1:"7")
wts.df2=wts.df2 %>%mutate(place = case_when(sample %in% c("f1","f2","f3","f4","f5","f6","f7","f8","f9") ~ as.character("Water column"),TRUE ~ as.character("Sediments")))
wts.df2$sample<-str_replace(wts.df2$sample,"f","")


##all

wgs.df2$Specific.pathway <- sapply(wgs.df2$Specific.pathway, function(x) paste( x,' WGS'))
wts.df2$Specific.pathway <- sapply(wts.df2$Specific.pathway, function(x) paste( x,' WTS'))


wgs_wts<-rbind(wgs.df2,wts.df2)
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"Rubisco","Calvin cycle")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"wood_lungdahl_pathway","wood lungdahl pathway")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"dehydrogenases_","dehydrogenases ")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"oxide_reduction","oxide reduction")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"Sodium","Na+")
wgs_wts$Specific.pathway<-str_replace(wgs_wts$Specific.pathway,"_"," ")

wgs_wts$place = factor(wgs_wts$place, levels=c('Water column','Sediments'))
wgs_wts$Class<-str_replace(wgs_wts$Class,"__$","") 
wgs_wts$Class<-str_replace(wgs_wts$Class,"_$","") 
wgs_wts$Class<-str_replace(wgs_wts$Class,"unknown$","Unknown") 
wgs_wts$Specific.pathway=str_replace(wgs_wts$Specific.pathway,"nitr^","Nitr")
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


#osmoregulation
pal<-iwanthue(112)
a10<-ggplot(wgs_wts[grep("Osmoregulation",wgs_wts$General.function), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Osmoregulation", x = "Sites", y = "RPKM")
p10=set_panel_size(a10,width = unit(5.5, "cm"),height = unit(4.5, "cm"))

#photosynthesis
pal<-iwanthue(24)
a1<-ggplot(wgs_wts[grep("Photosynthesis",wgs_wts$General.function), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Photosynthesis", x = "Sites", y = "RPKM")
p1=set_panel_size(a1,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#carbon
pal<-iwanthue(50)
a2<-ggplot(wgs_wts[grep("Carbon fixation",wgs_wts$General.function), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Carbon fixation", x = "Sites", y = "RPKM")
p2=set_panel_size(a2,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#nitrogen_nitrification
v=wgs_wts[grep("Nitrogen",wgs_wts$General.function), ]
pal<-iwanthue(93)
a3<-ggplot(v[grepl("Cyanophycin|fixation",v$Specific.pathway), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=10))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Nitrogen", x = "Sites", y = "RPKM")
p3=set_panel_size(a3,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

a4<-ggplot(v[grepl("Ammonia|Complete|Nitrite oxidation",v$Specific.pathway), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=10))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Nitrification", x = "Sites", y = "RPKM")
p4=set_panel_size(a4,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

a5<-ggplot(v[grepl("Nitrite reduction|DNRA|Nitrate reduction|Nitric oxide|Nitrous oxide|Nitrite reduction",v$Specific.pathway), ],col=pal, 
           aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=14))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "DNRA/Denitrification", x = "Sites", y = "RPKM")
p5=set_panel_size(a5,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#phosphorus
pal<-iwanthue(107)
a6<-ggplot(wgs_wts[grep("Phosphatase|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",wgs_wts$Specific.pathway), ],col=pal, 
           aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Phosphorus", x = "Sites", y = "RPKM")
p6=set_panel_size(a6,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#sulfur
pal<-iwanthue(88)
a7<-ggplot(wgs_wts[grep("Sulfate reduction|Sulfite reduction",wgs_wts$Specific.pathway), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfate and sulfite reduction", x = "Sites", y = "RPKM")
p7=set_panel_size(a7,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

pal<-iwanthue(88)
a8<-ggplot(wgs_wts[grep("Sulfide oxidation|Sulfur oxidation|ThioSulfate",wgs_wts$Specific.pathway), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfur, sulfide and thiosulfate oxidation", x = "Sites", y = "RPKM")
p8=set_panel_size(a8,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

##dsr rdsr
v=wgs_wts[grep("Sulfur",wgs_wts$General.function), ]
pal<-iwanthue(88)
a9<-ggplot(v[grep("dsr",v$Gene), ],col=pal, aes(sample, values, fill = Class)) +
  theme_bw()  +
  geom_bar(position="stack", stat="identity") +
  facet_grid( Gene+Specific.pathway~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
  scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
  theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
        legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
        strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfite reduction and sulfur oxidation", x = "Sites", y = "RPKM")
p9=set_panel_size(a9,width = unit(5.5, "cm"),height = unit(5.5, "cm"))


pdf("comb_taxa_function_barplot.pdf",width=18,height=32)
marrangeGrob(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10),  nrow = 1,ncol = 1)
dev.off()