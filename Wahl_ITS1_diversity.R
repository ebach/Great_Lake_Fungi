#Elizabeth Bach
#Code to calculate diversity stats for ITS1 sequences from Great Lakes sediment samples
#6 November 2017

#rm() removes all designated objects, this pervents confusion with any previous data you've been working with before running this code
#if you are begining a new R session, this is not necessary, but it won't hurt anything to run it
rm(list=ls())

#libarary() will open the packages necessary for the commands we will be using in this code
#if you do not have these packages intalled, you may intall through the menu or by running install.pacakges()
#You only need to install packages once, they need only be opened each time you start a new R session
library(lme4)
library(lmerTest)
library(bbmle)
library(reshape)

library(vegan)
library(ggplot2)
library(plyr)

#Read in the dataset, file.choose() will open a window from which you can select your data file: Wahl_ITS1_OTUtable_EMB.csv
data.nosing<-read.csv(file.choose())
head(data.nosing[,1:10])
str(data.nosing)
dim(data.nosing)
#remove taxonomic info
data.nosing2<-data.nosing[,-c(51:57)]
str(data.nosing2)

#add metadata to summarize results
data.meta<-read.csv(file.choose())
str(data.meta)
#merge metadata with sequence reads
data.nosing3<-merge(data.meta, data.nosing2, by="SampleID")
str(data.nosing3)
head(data.nosing3)

#pull out blanks
levels(data.nosing3$SampleID)
data.nosing3b<-droplevels(subset(data.nosing3, data.nosing3$SampleID!="HWblank1"))
data.nosing4<-droplevels(subset(data.nosing3b, data.nosing3$SampleID!="HWblank2"))
str(data.nosing4)
levels(data.nosing4$SampleID)

#calculating richness, shannons, and evenness
#Exclude metadata columns from analysis (in this case, columns 1:20), modify for your dataset
richness<-rowSums(data.nosing4[,-c(1:20)]>0)
#Preview
head(richness)
shannons<-diversity(data.nosing4[,-c(1:20)])
#Preview
head(shannons)
#note, other diversity indices can be calculated using diversity(), use ?diversity() in console to see other options

evenness<-shannons/log(richness)
#Preview
head(evenness)

#bind diversity stats into new data frame with just the summary statistics and metadata
div_stats<-data.frame(data.nosing4[,c(1,3,5,8)],richness,shannons,evenness)
#Preview
head(div_stats)
str(div_stats)

#looking at data distribution, this will help determine if transformations would be appropriate
ggplot(div_stats)+geom_histogram(aes(shannons))
#Shannons is best without transformation
ggplot(div_stats)+geom_histogram(aes(log(richness)))
#richness is less skewed if transformed, but not terrible if not transformed
ggplot(div_stats)+geom_histogram(aes(evenness))
#evenness is better withour transformation

#visualize data
#"quick & dirty" graphs
#Richness, depth group
rich<-ddply(div_stats, .(DepthGroup), summarise, .progress="text", mean=mean(richness), SE=sd(richness)/sqrt(length(richness)-1))
rich
ggplot(rich, aes(DepthGroup, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

#Richness, lake, just for visualization, not a valid comaprison as unequal sampling effort between the two lakes
rich<-ddply(div_stats, .(Lake), summarise, .progress="text", mean=mean(richness), SE=sd(richness)/sqrt(length(richness)-1))
rich
ggplot(rich, aes(Lake, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

#Richness, time, for FYI
rich<-ddply(div_stats, .(Date), summarise, .progress="text", mean=mean(richness), SE=sd(richness)/sqrt(length(richness)-1))
rich
ggplot(rich, aes(Date, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

#Evenness
rich<-ddply(div_stats, .(DepthGroup), summarise, .progress="text", mean=mean(evenness), SE=sd(evenness)/sqrt(length(evenness)-1))
ggplot(rich, aes(DepthGroup, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

range(div_stats$evenness)

#Shannon's
rich<-ddply(div_stats, .(DepthGroup), summarise, .progress="text", mean=mean(shannons), SE=sd(shannons)/sqrt(length(shannons)-1))
ggplot(rich, aes(DepthGroup, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

range(div_stats$shannons)


#Statistical tests:  Did not run any of these for the paper...

#general linear model ANOVA, testing main effects of all factors on diversity measures
#Using the * between factors will result in a full model testing each main effect and all possible interactions
#Using a + between factors will run a main effects model, ignoring any interactions
#You can customize to include only significant interactions and main effects by mixing and matching + and : 
#(e.g. Factor1+Factor2+Factor3+Factor1:Factor2 will evaluate main effects of factors 1,2,3 and the interaction of factors 1,2 only)
summary(test<-aov(richness~Factor1*Factor2, data=div_stats))
TukeyHSD(test)

summary(test2<-aov(evenness~Factor1*Factor2, data=div_stats))
TukeyHSD(test2)

summary(test3<-aov(shannons~Factor1*Factor2, data=div_stats))
TukeyHSD(test3)

#Mixed models for studies with blocking
#null model (no effect of any factors)
summary(test.null<-lmer(richness~(1|), data=div_stats, REML=FALSE))
#main effects model, with block
summary(test.main<-lmer(richness~Factor1+Factor2+Factor3+(1|block), data=div_stats, REML=FALSE))
#full model (all interactions)
summary(test.full<-lmer(richness~Factor1*Factor2*Factor3+(1|block), data=div_stats, REML=FALSE))
#you can remove non-significant interactions and include on the significant ones, eg
summary(test.part<-lmer(richness~Factor1+Factor2+Factor3+Factor1:Factor2+(1|block), data=div_stats, REML=FALSE))
#compare tests to see which model best fits the data (AIC comaprison), best fitting will have lowest AIC (by at least 5 usually)
AIC(test.null, test.main, test.full, test.part)
#Can run ANOVA on models to see if one fits "significantly" better than another
anova(test.null, test.main, test.full, test.part)
#Once best model is selected, run the model to see output and get the statistics for the data
#This is what you will report in a manuscript
summary(test.main)

#Nested factor
#main effects model, with block, Factor 2 is nested within Factor 1
summary(test.main<-lmer(richness~Factor1+Factor2+Factor3+(1|block)+(1|Factor1/Factor2), data=div_stats, REML=FALSE))

