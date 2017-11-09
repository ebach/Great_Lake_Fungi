#Elizabeth Bach, Ryan Williams
#code for processing QIIME OTU table output so can be used for subsequent analyses
#25 May 2015

#removes all objects previously assigned in R session, ensures no conflicting objects from another analysis
rm(list=ls())

#all useful packages, install.packages(...package name...) if you don't have them
library(labdsv)
library(vegan)
library(plyr)
library(reshape)
library(ggplot2)

# Read in the data, file.choose() will open browser window and allow you to navigate to data file
# QIIME output typcially .txt file rather than csv, so requires read.delim
dataset<-read.delim(file.choose(), skip=2)

#preview dataset to see if it looks like what you expect it should
dataset[1:10,20000:20010]
str(dataset)
dim(dataset)

#First transpose data so samples are rows and OTUs are columns
#create matrix with sample names (column names)
SampleID<-matrix(colnames(dataset[,-52]))
#object with OTU names as characters (R may assume these are numbers)
OTU<-(c("SampleID", as.character(dataset$X.OTU.ID)))

#transpose dataset
#remove collumn 52, which is the taxonomic information, this prevents the taxonomic text from confusing the dataset when transposed
dataset2<-as.data.frame(t(dataset[,-c(1,52)]))
#bind SampleID as first column so you can manipulate by SampleID as needed
dataset3<-cbind(SampleID[-1,], dataset2)
#Set object OTU as column names
colnames(dataset3)<-OTU

#preview transposed dataset3 to ensure it looks correct
dim(dataset3)
str(dataset3[1:10])
dataset3[1:10,1:10]

#let look for singletons
reads.otu<-colSums(dataset3[,-c(1)])
range(reads.otu)
reads.otu
subset(reads.otu, reads.otu<5)
dataset4<-rbind(reads.otu, dataset3)


# First we remove the singletons using the dropspc() function form the labdsv package.  In the line below I bind metadata (column 1)
# and the singleton-removed dataset.  Note that for the dropspc function I only use columns 2 through 7431 as these are all the OTUs (no metadata),
# and I remove species that occur 3 or less times

data.nosing<-cbind(dataset3[1],dropspc(dataset3[,-c(1)],minab=3))
str(data.nosing)
dim(data.nosing)
data.nosing[1:10,1:10]

# Now its time to figure out how many reads to rarefy by...I added a column to our dataset of the total number of reads per sample (row)
reads<-rowSums(data.nosing[,-c(1)])
data.nosing.reads<-cbind(reads,data.nosing)
head(data.nosing.reads[,1:10])
max(data.nosing.reads$reads)
#Sample read totals
Sample.reads<-ddply(data.nosing.reads, .(SampleID), summarise, .progress="text", total_read=sum(reads), avg_read=mean(reads))
Sample.reads
#look at distribution of reads
hist(reads)

#create rarefaction curves for each sample starting around 5000
rared<-rarefy(subset(data.nosing.reads, reads > 4999)[,-c(1:6)],sample=c(1,10,25,50,75,100,250,500,700,1250,2500,5000),se=FALSE )
rared_melt<-melt(rared)
names(rared_melt)<-c("sample","sample_size","OTUs")
rared_melt$sample_size<-c(rep(1,111),rep(10,111),rep(25,111),rep(50,111),rep(75,111),rep(100,111),rep(250,111),rep(500,111),rep(700,111),rep(1250,111),rep(2500,111),rep(5000,111))
head(rared_melt)

ggplot(rared_melt)+geom_line(aes(x=sample_size,y=OTUs,colour=sample,group=sample))+theme(aspect.ratio=1)+theme_bw()

#We ran the above code at several potential rarification cut-offs and decided rarifying to 10,000 reads would elimnate a few low-read samples,
#and retain most OTU data
data.reads2<-subset(data.nosing.reads, data.nosing.reads$reads>10000)
dim(data.nosing.reads)
dim(data.reads2)
min(data.reads2$reads)
mean(data.reads2$reads)

#rarifying to 10000 reads per sample
head(data.nosing[,1:5])
data.nosing.rar<-cbind(subset(data.nosing, reads > 9999)[,1:2],rrarefy(subset(data.nosing,reads > 9999)[,-c(1:2)],10000))
head(data.nosing.rar[,1:10])

#create .csv file of rarified data for subsequent analyses
write.csv(data.nosing.rar, file="COBS_ITS_data_rar.csv")