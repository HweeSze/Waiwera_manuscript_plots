##install and load libraries
remotes::install_github("adamlilith/enmSdm")
install.packages("remotes")
install.packages("rdist")
install.packages("emmeans")
library(remotes)
library(enmSdm)
library(betapart)
library(ggpmisc)
library(lattice)
library(rdist)
library(ggsci)
library(dplyr)
library(DescTools)
library(hues)
library(emmeans)

##read files
fil_dist<-read.delim("raw_Fil_distance.txt",header = FALSE,row.names = 1)
sed_dist<-read.delim("raw_sed_distance.txt",header = FALSE,row.names = 1)
dissimilarity<-read.csv("bray-curtis-diss.distance.matrix.csv",row.names = 1)
map<-read.delim("mappingcomb.txt", row.names=1, header=T)

### Water column distance decay
#using lat and long to get distance matrix
fil_dist_matrix<-pointDist(fil_dist,longLat = c(2,1))
fil_dist_matrix<-fil_dist_matrix/1000
fil_dist_matrix[is.na(fil_dist_matrix)] <- ""
fil_dist_matrix<-as.dist(fil_dist_matrix)

#dissimilarity matrix
dissimilarity[is.na(dissimilarity)] <- ""
dissimilarity_Filter<- as.matrix(dissimilarity[1:9,1:9])
dissimilarity_Filter<-as.dist(dissimilarity_Filter)

plot(fil_dist_matrix, dissimilarity_Filter, ylim=c(0,1), xlim=c(0, 5.2))

#try decay model
BCI.decay.exp<-decay.model(dissimilarity_Filter, fil_dist_matrix, y.type="dissim", model.type="exp", perm=100)
BCI.decay.pow<-decay.model(dissimilarity_Filter, fil_dist_matrix, y.type="dissim", model.type="pow", perm=100)

new_fil<-BCI.decay.exp$data
colnames(new_fil) <-c("x","y")
fit_fil<-(glm(x ~ y, data=new_fil, family = "gaussian"))
filter_r2<-PseudoR2(fit,which = "Nagelkerke")
summary(fit_fil)
filter_p<-1.752e-07
fil_slope<-BCI.decay.pow$b.slope

### Sediment distance decay
sed_dist_matrix<-pointDist(sed_dist,longLat = c(2,1))
sed_dist_matrix<-sed_dist_matrix/1000
sed_dist_matrix[is.na(sed_dist_matrix)] <- ""
sed_dist_matrix<-as.dist(sed_dist_matrix)

dissimilarity[is.na(dissimilarity)] <- ""
dissimilarity_sed<- as.matrix(dissimilarity[10:36,10:36])
dissimilarity_sed<-as.dist(dissimilarity_sed)

BCI.decay.exp.sed<-decay.model(dissimilarity_sed, sed_dist_matrix, y.type="dissim", model.type="exp", perm=100)
BCI.decay.pow.sed<-decay.model(dissimilarity_sed, sed_dist_matrix, y.type="dissim", model.type="pow", perm=100)

new_sed<-BCI.decay.exp.sed$data
colnames(new_sed) <-c("x","y")
fit<-(glm(x ~ y, data=new_sed, family = "gaussian"))
sed_r2<-PseudoR2(fit,which = "Nagelkerke")
summary(fit)
sed_p<-2e-16
sed_slope<-BCI.decay.pow.sed$b.slope


##combined
new_fil<-mutate(new_fil,state=c("ff",rep("fb",5),rep("fm",2),rep("fb",5),rep("fm",2),rep("bb",4),rep("bm",2),rep("bb",3),rep("bm",2),rep("bb",2),rep("bm",2),rep("bb",1),rep("bm",2),rep("bm",2),"mm"))
new_sed<-mutate(new_sed,state=c(rep("ff",5),rep("fb",15),rep("fm",6),rep("ff",4),rep("fb",15),rep("fm",6),rep("ff",3),rep("fb",15),rep("fm",6),rep("ff",2),rep("fb",15),rep("fm",6),rep("ff",1),rep("fb",15),rep("fm",6),rep("fb",15),rep("fm",6),rep("bb",14),rep("bm",6),rep("bb",13),rep("bm",6),rep("bb",12),rep("bm",6),rep("bb",11),rep("bm",6),rep("bb",10),rep("bm",6),rep("bb",9),rep("bm",6),rep("bb",8),rep("bm",6),rep("bb",7),rep("bm",6),rep("bb",6),rep("bm",6),rep("bb",5),rep("bm",6),rep("bb",4),rep("bm",6),rep("bb",3),rep("bm",6),rep("bb",2),rep("bm",6),rep("bb",1),rep("bm",6),rep("bm",6),rep("mm",5),rep("mm",4),rep("mm",3),rep("mm",2),rep("mm",1)))
##combine both
new_all<-rbind(new_fil,new_sed)
new_all<-mutate(new_all,type=c(rep("Filter",36),rep("Sediment",351)))

##test significance of covariance (analysis of covariance ancova)
anova((lm(y~x*type+x:type,data=new_all)))
inter.lst<-lstrends(inter,"type",var="x")
pairs(inter.lst)

#plot
pal<-iwanthue(6)
ggplot(new_all, aes(x, y))+ ylim(0,1.2)+ 
  geom_point(data=new_all[grep("Sediment",new_all$type), ],colour="black",shape=1,size=4) +
  geom_smooth(data=new_sed, aes(x, y,linetype = "LOESS"),size=1.2,color="black",alpha=0.2,method = "loess")+
  geom_smooth(data=new_sed, aes(x, y,linetype = "GLM"),size=1.2,alpha=0.2,color="black",method = "glm",formula = y~log(x)) +
  geom_smooth(data=new_fil, aes(x, y,linetype = "LOESS"),size=1.2,color="blue",alpha=0.2,method = "loess")+
  geom_smooth(data=new_fil, aes(x, y,linetype = "GLM"),size=1.2,color="blue",alpha=0.2,method = "glm",formula = y~log(x)) +
  theme_bw()+ggtitle("Distance-Decay Plot")+ xlab("Distance (km)")+
  ylab("Bray-curtis dissimilarity")+geom_point(data=new_all[grep("Filter",new_all$type), ],shape=19,mapping = aes(colour = state),size=4) +
  theme(plot.title = element_text(hjust = 0.5,size=15),  legend.title=element_text(size=12), 
        legend.text=element_text(size=11),axis.title=element_text(size=12),axis.text=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = c(0.85, 0.25))+
  geom_text(aes(x=1.3,y=1.1,label = paste("Sediment GLM Adj R2 = ",signif(sed_r2, 5),"\nGLM P-value =",signif(sed_p, 5),"\nGLM slope =",signif(sed_slope, 5)))) +
  geom_text(aes(x=1.3,y=0.1,label = paste("Water Column GLM Adj R2 = ",signif(filter_r2, 5),"\nGLM P-value =",signif(filter_p, 5),"\nGLM slope =",signif(fil_slope, 5))))+
  scale_colour_manual(name="legend", values=pal)+
  labs(shape = "Type", colour = "Line colour")
