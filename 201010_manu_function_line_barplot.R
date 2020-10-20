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

#read files
map<-read.delim2("mappingcomb.txt")
taxa_wgs_all<-read.delim2("taxa_wgs_all.txt",sep = "\t")
taxa_wts_all<-read.delim2("taxa_wts_all.txt",sep = "\t")

#Convert to numeric
taxa_wgs_all[,12:47] <- lapply(taxa_wgs_all[,12:47], as.character)
taxa_wts_all[,12:41] <- lapply(taxa_wts_all[,12:41], as.character)

taxa_wgs_all[,12:47] <- lapply(taxa_wgs_all[,12:47], as.numeric)
taxa_wts_all[,12:41] <- lapply(taxa_wts_all[,12:41], as.numeric)

#convert taxa to character
taxa_wgs_all$taxac<-as.character(taxa_wgs_all$taxac)

### preparation for lineplot
line_wgs <- select(taxa_wgs_all,X2,X3, f1:S9R3) %>%
    group_by(X2,X3) %>%
    summarise_at(vars(f1:S9R3), ~sum(.)) %>%
    gather(key = "sample", value = "values", f1:S9R3)

line_wts <- select(taxa_wts_all,X2,X3, f1:S7R3) %>%
    group_by(X2,X3) %>%
    summarise_at(vars(f1:S7R3), ~sum(.)) %>%
    gather(key = "sample", value = "values", f1:S7R3)

line_wgs$X3 <- sapply(line_wgs$X3, function(x) paste( x,' WGS'))
line_wts$X3 <- sapply(line_wts$X3, function(x) paste( x,' WTS'))

#combine wgs wts
line_wgs_wts<-rbind(line_wgs,line_wts)

#log rpm and rpkm counts
filter_line_wgs_wts<-left_join(line_wgs_wts,map,by=c("sample"="X")) %>% filter(values !=0) %>% mutate(values = values+1) %>% mutate(log=log10(values))
filter_line_wgs_wts<-filter_line_wgs_wts%>%mutate(place = case_when(place=="Filter" ~ as.character("Water Column"),TRUE ~ as.character("Sediments")))
filter_line_wgs_wts$salinity<-as.numeric(as.character(filter_line_wgs_wts$salinity))

filter_line_wgs_wts <-as.data.frame(filter_line_wgs_wts)  %>%mutate(X3= case_when(X3 =="Rubisco WGS" ~ "Calvin cycle WGS", TRUE ~as.character(X3)))
filter_line_wgs_wts <-as.data.frame(filter_line_wgs_wts)  %>%mutate(X3= case_when(X3 =="Rubisco WTS" ~ "Calvin cycle WTS", TRUE ~as.character(X3)))

linedata<-separate(data = filter_line_wgs_wts, col = X3, into = c("X3", "right"), sep = " W")
filter_line_wgs_wts$salinity<-as.numeric(as.character(filter_line_wgs_wts$salinity))

linedata$X3 = str_replace(linedata$X3,"_"," ")
linedata$right = str_replace(linedata$right,"GS","WGS")
linedata$right = str_replace(linedata$right,"TS","WTS")
linedata$X3 = str_replace(linedata$X3,"Sodium−translocating dehydrogenases_nqrF","Na+ translocating dehydrogenase nqrF")
linedata$X3 = str_replace(linedata$X3,"Potassium","K+")

typecol<- c("white", "gray88")

linedata$place = factor(linedata$place, levels=c('Water Column','Sediments'))

##lineplot
g <- make_gradient(
  deg = 180, n = 500, cols = alpha(rev(brewer.pal(9, "Blues")),0.1))

d1<-ggplot(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",linedata$X3), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=right),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",linedata$X3),],right=="WGS"),
  color="black")+facet_rep_grid(X3~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Phosphatase|PhoR PhoB|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",linedata$X3),],right=="WTS"),
  color="darkgrey")+
  scale_y_continuous(name = "Log RPM (WGS) or Log RPKM (WTS)") +theme_bw() + 
  ggtitle("Phosphorus")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none")

d2<-ggplot(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|ThioSulfate",linedata$X3), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=right),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|ThioSulfate",linedata$X3),],right=="WGS"),color="black")+
  facet_rep_grid(X3~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Sulfate reduction|Sulfide|Sulfite|Sulfur|ThioSulfate",linedata$X3),],right=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log RPM (WGS) or Log RPKM (WTS)") +theme_bw() + 
  ggtitle("Sulfur")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d3<-ggplot(linedata[grep("N fixation|Cyanophycin|Ammonia|Complete|DNRA|nitrate|nitric|nitrite|nitrous",linedata$X3), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=right),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("N fixation|Cyanophycin|mmonia oxidation|DNRA|nitrate|nitric|nitrite|nitrous",linedata$X3),],right=="WGS"),
  color="black")+facet_rep_grid(X3~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("N fixation|Cyanophycin|mmonia oxidation|DNRA|nitrate|nitric|nitrite|nitrous",linedata$X3),],right=="WTS"),
  color="darkgrey")+
  scale_y_continuous(name = "Log RPM (WGS) or Log RPKM (WTS)") +theme_bw() + ggtitle("Nitrogen")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d4<-ggplot(linedata[grep("glycine|K+|Na+",linedata$X3), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=right),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("glycine|K+|Na+",linedata$X3),],right=="WGS"),color="black")+
  facet_rep_grid(X3~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("glycine|K+|Na+",linedata$X3),],right=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log RPM (WGS) or Log RPKM (WTS)") +theme_bw() + 
  ggtitle("Osmoregulation")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")

d5<-ggplot(linedata[grep("Photosynthesis",linedata$X2), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=right),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+
  geom_line(data=filter(linedata[grep("Photosynthesis",linedata$X2),],right=="WGS"),color="black")+
  facet_rep_grid(X3~place,labeller = label_wrap_gen(width=15), scales="free")+
  geom_line(data=filter(linedata[grep("Photosynthesis",linedata$X2),],right=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log RPM (WGS) or Log RPKM (WTS)") +theme_bw() + 
  ggtitle("Photosynthesis")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),legend.position = "none")

d6<-ggplot(linedata[grep("Carbon fixation",linedata$X2), ],col=pal, aes(salinity, log)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  geom_jitter(aes(color=right),shape=1,stroke=1.5)+
  scale_colour_manual(values = c("black","darkgrey"))+  
  geom_line(data=filter(linedata[grep("Carbon fixation",linedata$X2),],right=="WGS"),color="black")+
  facet_rep_grid(X3~place,labeller = label_wrap_gen(width=15),scales="free")+
  geom_line(data=filter(linedata[grep("Carbon fixation",linedata$X2),],right=="WTS"),color="darkgrey")+
  scale_y_continuous(name = "Log RPM (WGS) or Log RPKM (WTS)") +theme_bw() + 
  ggtitle("Carbon fixation")+
  theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),axis.ticks.y = element_line(color="black"),
  axis.ticks.x = element_line(color="black"),text = element_text(size=18,color="black"),strip.text.x = element_text(size = 12, color ="black" ),
  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")


pdf("wgs_wts_lineplot.pdf",width=20,height=17)
grid_arrange_shared_legend(d4,arrangeGrob(d6,d1,ncol=1, nrow=2,heights = c(0.35,0.65)),arrangeGrob(d5,d2, nrow=2, 
	ncol=1,heights = c(0.35,0.65)),arrangeGrob(d3, nrow=1, ncol=1,heights = c(1.0)),nrow = 1,ncol = 4)
dev.off()



#### preparation for individual function-taxa barplot
wgs.df2 <- select(taxa_wgs_all1,X2,X3,X4, taxac, 12:29) %>%
    group_by(X2,X3,X4, taxac) %>%
    summarise_at(vars(f1:'9'), ~sum(.)) %>%
    gather(key = "sample", value = "values",f1:'9')
wgs.df2=wgs.df2 %>%mutate(place = case_when(sample %in% c("f1","f2","f3","f4","f5","f6","f7","f8","f9") ~ as.character("Water column"),TRUE ~ as.character("Sediments")))
wgs.df2$sample<-str_replace(wgs.df2$sample,"f","")


wts.df2 <- select(taxa_wts_all1, X2,X3,X4, taxac, f1:"7") %>%
    group_by(X2,X3,X4, taxac) %>%
    summarise_at(vars(f1:"7"), ~sum(as.numeric(.))) %>%
    gather(key = "sample", value = "values", f1:"7")
wts.df2=wts.df2 %>%mutate(place = case_when(sample %in% c("f1","f2","f3","f4","f5","f6","f7","f8","f9") ~ as.character("Water column"),TRUE ~ as.character("Sediments")))
wts.df2$sample<-str_replace(wts.df2$sample,"f","")


##all

wgs.df2$X3 <- sapply(wgs.df2$X3, function(x) paste( x,' WGS'))
wts.df2$X3 <- sapply(wts.df2$X3, function(x) paste( x,' WTS'))


wgs_wts<-rbind(wgs.df2,wts.df2)
wgs_wts$X3<-str_replace(wgs_wts$X3,"Rubisco","Calvin cycle")
wgs_wts$X3<-str_replace(wgs_wts$X3,"wood_lungdahl_pathway","wood lungdahl pathway")
wgs_wts$X3<-str_replace(wgs_wts$X3,"dehydrogenases_","dehydrogenases ")
wgs_wts$X3<-str_replace(wgs_wts$X3,"oxide_reduction","oxide reduction")
wgs_wts$X3<-str_replace(wgs_wts$X3,"Sodium","Na+")
wgs_wts$X3<-str_replace(wgs_wts$X3,"_"," ")

wgs_wts$place = factor(wgs_wts$place, levels=c('Water column','Sediments'))
wgs_wts$taxac<-str_replace(wgs_wts$taxac,"__$","") 
wgs_wts$taxac<-str_replace(wgs_wts$taxac,"_$","") 
wgs_wts$taxac<-str_replace(wgs_wts$taxac,"unknown$","Unknown") 
wgs_wts$X3=str_replace(wgs_wts$X3,"nitr^","Nitr")
wgs_wts$X3=str_replace(wgs_wts$X3,"^Phosphate Regulon","PhoR PhoB Phosphate Regulon")
wgs_wts$X3=str_replace(wgs_wts$X3,"^PhoR PhoB Inhibitor Protein","PhoU Phosphate Regulon")


### Individual barplot
#osmoregulation
pal<-iwanthue(112)
a10<-ggplot(wgs_wts[grep("Osmoregulation",wgs_wts$X2), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Osmoregulation", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p10=set_panel_size(a10,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#photosynthesis
pal<-iwanthue(24)
a1<-ggplot(wgs_wts[grep("Photosynthesis",wgs_wts$X2), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Photosynthesis", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p1=set_panel_size(a1,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#carbon
pal<-iwanthue(50)
a2<-ggplot(wgs_wts[grep("Carbon fixation",wgs_wts$X2), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Carbon fixation", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p2=set_panel_size(a2,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#nitrogen_nitrification
v=wgs_wts[grep("Nitrogen",wgs_wts$X2), ]
pal<-iwanthue(93)
a3<-ggplot(v[grepl("Cyanophycin|fixation",v$X3), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=10))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Nitrogen", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p3=set_panel_size(a3,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

a4<-ggplot(v[grepl("Ammonia|Complete|Nitrite oxidation",v$X3), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=10))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Nitrification", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p4=set_panel_size(a4,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

a5<-ggplot(v[grepl("Nitrite reduction|DNRA|Nitrate reduction|Nitric oxide|Nitrous oxide|Nitrite reduction",v$X3), ],col=pal, 
    aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=14))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "DNRA/Denitrification", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p5=set_panel_size(a5,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#phosphorus
pal<-iwanthue(107)
a6<-ggplot(wgs_wts[grep("Phosphatase|Phosphate Inorganic Transporter|Phosphate Regulon|Specific Transport",wgs_wts$X3), ],col=pal, 
    aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Phosphorus", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p6=set_panel_size(a6,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

#sulfur
pal<-iwanthue(88)
a7<-ggplot(wgs_wts[grep("Sulfate reduction|sulfite reduction",wgs_wts$X3), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfate and sulfite reduction", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p7=set_panel_size(a7,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

pal<-iwanthue(88)
a8<-ggplot(wgs_wts[grep("Sulfide oxidation|Sulfur oxidation|ThioSulfate",wgs_wts$X3), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfur, sulfide and thiosulfate oxidation", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p8=set_panel_size(a8,width = unit(5.5, "cm"),height = unit(5.5, "cm"))

##dsr rdsr
v=wgs_wts[grep("Sulfur",wgs_wts$X2), ]
pal<-iwanthue(88)
a9<-ggplot(v[grep("dsr",v$X4), ],col=pal, aes(sample, values, fill = taxac)) +
    theme_bw()  +
    geom_bar(position="stack", stat="identity") +
    facet_grid( X4+X3~place , scales = "free_y",labeller = label_wrap_gen(width=15))+
    scale_fill_manual(values = pal)+scale_x_discrete(name="Distance (km)", labels=scaleFUN(as.numeric(as.character(map$Distance[1:9]))))+
    theme(legend.position="right",axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size =15),axis.title=element_text(size=18),
    legend.text=element_text(size=16),axis.text.y=element_text(size=15),plot.title = element_text(face = 'bold', size = 20,hjust=0.5),
    strip.text.x = element_text(size = 18, color ="black" ),strip.text.y = element_text(size = 18, color ="black" ),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    guides(fill=guide_legend(ncol=1))+ labs(title = "Sulfite reduction and sulfur oxidation", x = "Sites", y = "RPM (WGS) or RPKM (WTS)")
p9=set_panel_size(a9,width = unit(5.5, "cm"),height = unit(5.5, "cm"))


pdf("comb_taxa_function_barplot.pdf",width=18,height=32)
marrangeGrob(grobs=list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10),  nrow = 1,ncol = 1)
dev.off()