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

#Read in the dataset, file.choose() will open a window from which you can select your data file: Wahl_ITS2_OTUtable_EMB.csv
data.nosing<-read.csv(file.choose())
head(data.nosing[,1:10])
str(data.nosing)
dim(data.nosing)

#add metadata to summarize results
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
ggplot(div_stats)+geom_histogram(aes(richness))
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
even<-ddply(div_stats, .(DepthGroup), summarise, .progress="text", mean=mean(evenness), SE=sd(evenness)/sqrt(length(evenness)-1))
ggplot(even, aes(DepthGroup, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

range(div_stats$evenness)

#Shannon's
shan<-ddply(div_stats, .(DepthGroup), summarise, .progress="text", mean=mean(shannons), SE=sd(shannons)/sqrt(length(shannons)-1))
ggplot(shan, aes(DepthGroup, mean))+geom_pointrange(aes(ymax=(mean+SE), ymin=(mean-SE)))

range(div_stats$shannons)

