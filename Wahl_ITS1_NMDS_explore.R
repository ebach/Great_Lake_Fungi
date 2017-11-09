#Elizabeth Bach, Ryan Williams
#Exploritory NMDS to visualize community differences
#5 June 2015

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session
library(reshape)
library(grid)
library(ggplot2)
library(lme4)
library(lmerTest)
library(bbmle)
library(vegan)
library(gridExtra)

#For this exploritory NMDS, you will not need taxonomic information
#Read in the dataset, Wahl_ITS1_OTUtable_EMB.csv
#we used non-rarified data to explore the commmunity differences
#reading in non-rarified data
data.nosing<-read.csv(file.choose())
head(data.nosing[,1:10])
str(data.nosing)
dim(data.nosing)

#add metadata to summarize results
#Use GL_ITS1_HW_metaData.csv
data.meta<-read.csv(file.choose())
str(data.meta)
#merge metadata with sequence reads
data.nosing3<-merge(data.meta, data.nosing, by="SampleID")
str(data.nosing3)
head(data.nosing3)

#pull out blanks
levels(data.nosing3$SampleID)
data.nosing3b<-droplevels(subset(data.nosing3, data.nosing3$SampleID!="HWblank1"))
data.nosing4<-droplevels(subset(data.nosing3b, data.nosing3$SampleID!="HWblank2"))
str(data.nosing4)
levels(data.nosing4$SampleID)

#Presence/Absence MDS
#First we will look at community structure determined by OTU presence/absence, 
#we will convert the data frame into presnece absence using decostand(), "pa"=presence/absence, again, exclude all metadata so only OTUs are used
#we will run permutations of the multi-dimensional scaling (MDS) algorithm using mteaMDS(), using the Jaccard dissimilarity distance (for presence/absence data)
#we will run MDS ordinations for 6, 5, 4, 3, and 2 dimensions (k), check the stress scores to determine which dimensionality best describes the data, as a rule of thumb, choose the dimensionality that represents a strong decline in stress score, additionally, anything <0.10 is a good representation of the data, 0.10-0.20 may indicate some morphing of true data distances, but if shows the community differences, is acceptable for visualization (you will run statistics seperately anyway), I figure a 2 or 3-D plot with a stress score ~0.19 is better than trying to convey 4 or more dimenisions
#generally 2 or 3 dimensions is all humans can comprehend, if there is a strong reason you need to include 4 or more dimensions, please discuss with collleagues and/or statistician
mds.pa6<-metaMDS(decostand(data.nosing4[,-c(1:20)],"pa" ),distance="jaccard", k=6,autotransform=FALSE, na.rm=TRUE)
#now look at best stress score
mds.pa6
mds.pa5<-metaMDS(decostand(data.nosing4[,-c(1:20)],"pa" ),distance="jaccard", k=5,autotransform=FALSE, na.rm=TRUE)
mds.pa5
mds.pa4<-metaMDS(decostand(data.nosing4[,-c(1:20)],"pa" ),distance="jaccard", k=4,autotransform=FALSE, na.rm=TRUE)
mds.pa4
mds.pa3<-metaMDS(decostand(data.nosing4[,-c(1:20)],"pa" ),distance="jaccard", k=3,autotransform=FALSE, na.rm=TRUE)
mds.pa3
mds.pa2<-metaMDS(decostand(data.nosing4[,-c(1:20)],"pa" ),distance="jaccard", k=2,autotransform=FALSE, na.rm=TRUE)
mds.pa2
#Preview result you want to use
head(mds.pa2)

#co-ordinates for graphing the NMDS must be extracted from the mds objects above using scores()
#coalate mds scors into dataframe for graphing, be sure to adjust mds.pa to the correct number selected above
MDS1<-data.frame(scores(mds.pa2))$NMDS1
MDS2<-data.frame(scores(mds.pa2))$NMDS2
SampleID<-data.nosing$SampleID
NMDS.pa<-data.frame(MDS1,MDS2,SampleID)
head(NMDS.pa)

#Exploritory NMDS figure, this shows all data points
ggplot(NMDS.pa)+geom_point(aes(x=MDS1, y=MDS2, size=2))

#Now we will see if communities group together based on our factors of interest
#this code will produce ovals showing the center of the data spread for each factor of interes and color-code by factor level
#first, run this function, this must be run to produce the graph, but will not produce any output by itself
#graph for publication
ggplot.NMDS<-function(XX,ZZ,COLORS,SHAPES){
	library(ggplot2)
MDS1<-data.frame(scores(XX))$NMDS1
MDS2<-data.frame(scores(XX))$NMDS2
Treatment<-ZZ

NMDS<-data.frame(MDS1,MDS2,Treatment)

NMDS.mean=aggregate(NMDS[,1:2],list(group=Treatment),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }

X1<-ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment,shape=Treatment,),size=3,alpha=0.75) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), show.legend=FALSE, size=2, linetype=5)+theme_bw()+theme(aspect.ratio=1)+scale_color_manual(values=COLORS,breaks=c("CC","P","PF"),labels=c("Maize","Prairie","Fertilized Prairie"))+scale_shape_manual(values=SHAPES,breaks=c("CC","P","PF"),labels=c("Maize","Prairie","Fertilized Prairie"))+theme(axis.text.x=element_text(size=20),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+theme(legend.title=element_text(size=15),legend.text=element_text(size=15))
X1    
}

colors.bw<-c("black","grey")
shapes<-c(15,16)

PA.GL<-ggplot.NMDS(mds.pa2, (data.nosing4$Lake), colors.bw, shapes)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))
PA.GL



#Now lets look at abundance
#This code is the same as above, except we will relatives the abundance of each OTU in each sample using decostand() with "total", agian omit any metadata columns
#we will use Bray-Curtis distance for this NMDS, Bray-Curtis is the default for meatMDS() function, so we do not need to specify it
mds.ab6<-metaMDS(decostand(data.nosing4[,-c(1:20)],"total" ), k=6,autotransform=FALSE, na.rm=TRUE)
#now look at best stress score
mds.ab6
mds.ab5<-metaMDS(decostand(data.nosing4[,-c(1:20)],"total" ), k=5,autotransform=FALSE, na.rm=TRUE)
mds.ab5
mds.ab4<-metaMDS(decostand(data.nosing4[,-c(1:20)],"total" ), k=4,autotransform=FALSE, na.rm=TRUE)
mds.ab4
mds.ab3<-metaMDS(decostand(data.nosing4[,-c(1:20)],"total" ), k=3,autotransform=FALSE, na.rm=TRUE)
mds.ab3
mds.ab2<-metaMDS(decostand(data.nosing4[,-c(1:20)],"total" ), k=2,autotransform=FALSE, na.rm=TRUE)
mds.ab2
#Preview result you want to use
head(mds.ab2)

#Extract co-ordinates and coalate mds scors into dataframe
MDS1.ab<-data.frame(scores(mds.ab2))$NMDS1
MDS2.ab<-data.frame(scores(mds.ab2))$NMDS2
SampleID<-data.nosing4$SampleID
NMDS.ab<-data.frame(MDS1.ab,MDS2.ab,SampleID)
head(NMDS.ab)

#Exploritory NMDS figure, this shows all data points
ggplot(NMDS.ab)+geom_point(aes(x=MDS1, y=MDS2, size=2))

#Now looking at factor groupings, we do not need to re-run the ggplot.NMDS function, it is already saved in R's memory for this working session

AB.GL<-ggplot.NMDS(mds.ab2, (data.nosing4$Lake), colors.bw, shapes)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))
AB.GL

#Depth (only statistically significant for Bray-Curtis similarity)
colors.bw<-c("black","black","black","black")
colors.dp<-c(rgb(199,102,116, max=255), rgb(154,153,69, max=255), rgb(103,121,208, max=255), rgb(151,77,157, max=255))
shapes<-c(15,16,17,18)

AB.DP<-ggplot.NMDS(mds.ab2, (data.nosing4$DepthGroup), colors.dp, shapes)+
theme(axis.line=element_line(size=1.25), aspect.ratio=1, panel.border=element_blank(),axis.ticks=element_line(size=1.25, colour="black"), panel.background=element_blank(), axis.text=element_text(size=10, face="bold", colour="black"), axis.title=element_text(size=12, face="bold", colour="black"))
AB.DP

#If differences in your factor levels are more obvious in the presence/absence figure, that means your communities mostly differ by which OTUs are there or not there
#If differences in factor levels are more obvious in the abundance figure, that means which taxa are present in the community may be similar, but their relative abundances are different

#To determine if your community types are statistically different, see Bach_ITS1_multivariate.R
