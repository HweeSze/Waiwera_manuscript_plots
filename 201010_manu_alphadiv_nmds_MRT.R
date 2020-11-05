###(Install and) Load libraries:
library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(ggrepel)
library(gtools)
###Read in files:
rar<-read.delim("rar.filt.sintax.tab.otu-table.97.all.cov.norm.emirge.txt", row.names=1, header=T)
map<-read.delim("mappingcomb.txt", row.names=1, header=T)

###Rearrange columns alphabetically:
srar<-cbind(rar[ , mixedsort(names(rar[,1:36]))],rar[,37:43])
###Sample columns:
srar[1,1:36]
srar<-mutate(srar, sum=rowSums(srar[,1:36]),rowname=row.names(srar))

level<-c(rep("Freshwater",2),rep("Brackish",5),rep("Marine",2))
level.comb<-factor(map$state,levels=c("Freshwater", "Brackish",   "Marine"))
map$place_f = factor(map$place, levels=c('Water column','Sediments'))

#filter bacteria,archaea, and microeukaryotes into different dataframe
baconly<-srar %>% filter(Kingdom=="k:Bacteria")
row.names(baconly)<-baconly$rowname

archaea<-srar %>% filter(Kingdom=="k:Archaea")
row.names(archaea)<-archaea$rowname

eukonly<-srar %>% filter(Kingdom=="k:Eukaryota")
row.names(eukonly)<-eukonly$rowname

# Alpha diversity
## Shannon indices for all, microeukaryotes, bacteria and  archaea
shannonrar<- data.frame(x=colnames(srar[1:36]),Y=diversity(t(srar[, 1:36]) ,index="shannon"))
shannonrar_bac<- data.frame(x=colnames(baconly[1:36]),Y=diversity(t(baconly[, 1:36]) ,index="shannon"))
shannonrar_eu<- data.frame(x=colnames(eukonly[1:36]),Y=diversity(t(eukonly[, 1:36]) ,index="shannon"))
shannonrar_arc<- data.frame(x=colnames(archaea[1:36]),Y=diversity(t(archaea[, 1:36]) ,index="shannon"))

#set color
pal <- c("#b4ecf9", "#1e79da", "#0a0c93")

div_all<-ggplot(shannonrar, aes(x = level.comb, y = Y)) +
	geom_boxplot(fill=rep(pal,2),outlier.shape=NA,colour ="black") + scale_x_discrete(name = "Samples") +
	scale_y_continuous(name = "Shannon index") +theme_bw() +facet_wrap( ~ map$place_f ) +
	ggtitle("All")+   geom_jitter()+
	theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),
	text = element_text(size=18,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	strip.text.x = element_text(size = 20, color ="black" ))
div_e<-ggplot(shannonrar_eu, aes(x = level.comb, y = Y)) +
    geom_boxplot(fill=rep(pal,2),outlier.shape=NA,color="black") + scale_x_discrete(name = "Samples") +
    scale_y_continuous(name = "Shannon index") +theme_bw() +facet_wrap( ~ map$place_f ) +
    ggtitle("Microeukaryote")+   geom_jitter()+
    theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),
    text = element_text(size=18,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20, color ="black" ))
div_ar<-ggplot(shannonrar_arc, aes(x = level.comb, y = Y)) +
    geom_boxplot(fill=rep(pal,2),outlier.shape=NA,color="black") + scale_x_discrete(name = "Samples") +
    scale_y_continuous(name = "Shannon index") +theme_bw() +facet_wrap( ~ map$place_f ) +
    ggtitle("Archaea")+   geom_jitter()+
    theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),
    text = element_text(size=18,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20, color ="black" ))
div_b<-ggplot(shannonrar_bac, aes(x = level.comb, y = Y)) +
    geom_boxplot(fill=rep(pal,2),outlier.shape=NA,color="black") + scale_x_discrete(name = "Samples") +
    scale_y_continuous(name = "Shannon index") +theme_bw() +facet_wrap( ~ map$place_f ) +
    ggtitle("Bacteria")+   geom_jitter()+
    theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),
    text = element_text(size=18,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20, color ="black" ))

## Relative OTU richness for all, microeukaryotes, bacteria and  archaea
specrichness<-data.frame(x=colnames(srar[1:36]),y=specnumber(t(srar[, 1:36])))
spec_bac<- data.frame(x=colnames(baconly[1:36]),Y=specnumber(t(baconly[, 1:36])))
spec_eu<- data.frame(x=colnames(eukonly[1:36]),Y=specnumber(t(eukonly[, 1:36])))
spec_arc<- data.frame(x=colnames(archaea[1:36]),Y=specnumber(t(archaea[, 1:36])))

spe_all<-ggplot(specrichness, aes(x = level.comb, y = y)) +
	geom_boxplot(fill=rep(pal,2),outlier.shape=NA,colour="black") + scale_x_discrete(name = "Samples") +
	scale_y_continuous(name = "Relative OTU Richness") +theme_bw() + facet_wrap( ~ map$place_f ) +
	geom_jitter()+ggtitle("All")+
	theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),text = element_text(size=18),
	panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text.x = element_text(size = 20, colour ="black" ))
spe_e<-ggplot(spec_eu, aes(x = level.comb, y = Y)) +
    geom_boxplot(fill=rep(pal,2),outlier.shape=NA,color="black") + scale_x_discrete(name = "Samples") +
    scale_y_continuous(name = "Relative OTU Richness") +theme_bw() +facet_wrap( ~ map$place_f ) +
    ggtitle("Microeukaryote")+   geom_jitter()+
    theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),text = element_text(size=18,color="black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text.x = element_text(size = 20, color ="black" ))
spe_ar<-ggplot(spec_arc, aes(x = level.comb, y = Y)) +
    geom_boxplot(fill=rep(pal,2),outlier.shape=NA,color="black") + scale_x_discrete(name = "Samples") +
    scale_y_continuous(name = "Relative OTU Richness") +theme_bw() +facet_wrap( ~ map$place_f ) +
    ggtitle("Archaea")+   geom_jitter()+
    theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),text = element_text(size=18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text.x = element_text(size = 20, color ="black" ))
spe_b<-ggplot(spec_bac, aes(x = level.comb, y = Y)) +
    geom_boxplot(fill=rep(pal,2),outlier.shape=NA,color="black") + scale_x_discrete(name = "Samples") +
    scale_y_continuous(name = "Relative OTU Richness") +theme_bw() +facet_wrap( ~ map$place_f ) +
    ggtitle("Bacteria")+   geom_jitter()+
    theme(plot.title = element_text(hjust=0.5),axis.text = element_text(color="black"),text = element_text(size=18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.text.x = element_text(size = 20, color ="black" ))
 

pdf("alphadiv_all_boxplot.point.pdf",colormodel="cmyk",width=22,height=7,paper='special')
grid.arrange(div_all, div_ar, div_b, div_e, spe_all, spe_ar, spe_b, spe_e, ncol=4)
dev.off()

##Non-multidimensional analysis
###3. Create ordination
nmdssrar <- metaMDS(t(srar[,1:36]), distance="bray", k=2, trymax=100, autotransform=TRUE, zerodist="add")
nmdssrar
efsrar <- envfit(nmdssrar, map[,c(4,6,11,12,13,15,16,17)], permu = 999,na.rm=T)
efsrar

##Check stress:
stressplot(nmdssrar, pch=1, p.col="blue", l.col="red", lwd=2)
legend("bottomright", legend = paste("Stress: ", round(nmdssrar$stress,digits=4 )), bty = "n", cex = 0.8)

mds.plot <- list()
mds.plot[["points"]] <- scores(x = nmdssrar, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample")
mds.plot[["vector"]] <- scores(x = efsrar, display = "vectors") %>%
    as.data.frame() %>%
    rownames_to_column(var = "vector") %>%
    rename(v1 = NMDS1, v2 = NMDS2)
map$sample <- rownames(map)
mds.plot[["points"]] <- mds.plot[["points"]] %>%
    left_join(select(map, sample, place, salinity), by = "sample")

#select  vector with p-value <0.05
vector1<-mutate(mds.plot$vector,efsrar$vectors$pvals)
vector1<-filter(vector1,efsrar$vectors$pvals<=0.05)
row.names(vector1)<- vector1$vector

pdf("nmds_community.pdf",colormodel="cmyk",width=7,height=5,paper='special')
ggplot(data = mds.plot[["points"]],
       mapping = aes(x = NMDS1, y = NMDS2)) +
    geom_point(mapping = aes(shape = place, colour = salinity),
               size = 4,stroke=2) + theme_bw()+
    ggtitle("NMDS (Rarefied)") +
    geom_segment(data = vector1,
                 mapping = aes(x = 0, y = 0,xend = v1* 0.65, yend = v2 * 0.65),
                 arrow = arrow(length = unit(0.03, "npc"))) +
    geom_text_repel(data = vector1,size=5,segment.color = 'transparent',
                    mapping = aes(x = v1 * 0.7, y = v2 * 0.7,label = vector,)) +
    scale_colour_gradient(low = "#b4ecf9", high = "#0a0c93") +
    scale_shape_manual(values=c(1,19))+
    labs(shape = "Type", colour = "Salinity (ppt)")+ theme(plot.title = element_text(hjust=0.5,size=14),
    axis.text=element_text(color = "black",size=12),axis.title=element_text(color = "black",size=13),
    legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_y_reverse()+ scale_x_continuous(limits = c(-0.9,0.9))
dev.off()

## Multivariate regression tree
library(mvpart)
library(vegan)
library(plyr)
environ <- map[,c(4:6,10:22)]
EnvironRowNames <- rownames(environ)
environ <- lapply(environ,function(x) as.numeric(as.character(x)))
environ <-as.data.frame(environ)
rownames (environ) = EnvironRowNames
environ.df = mutate(environ, place = map$place)
environ.df = mutate(environ.df, state = map$state)
environ.df = mutate(environ.df, comb = map$comb)
environ.df <-as.data.frame(environ.df)
abundance=data.matrix(t(srar[,1:36]))

formula <- abundance ~ salinity + Distance + state + place + Nitrate + DRP + Sulphate + 
  DNPOC + Nitrite + Ammonium + pH + Temp  + DNPOC + Mud.content + TRP + TS +
  TS + TN + TOC 

set.seed(1234)
mvpart_run2 <- mvpart(
  form = formula, 
  data = environ.df,
  rsq = TRUE,  # give "rsq" plot
  pca = TRUE,  # plot PCA of group means and add species and site information
  wgt.ave.pca = TRUE,  # plot weighted averages acorss sites for species
  size=8,bar=F, legend = T, pretty = "full",xv="1se",xval=1000,
)


