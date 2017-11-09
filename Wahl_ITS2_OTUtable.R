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
#dataset<-read.delim(file.choose(), skip=5)
dataset<-read.csv(file.choose())

#preview dataset to see if it looks like what you expect it should
dataset[1:10,1:10]
str(dataset)
dim(dataset)

#First transpose data so samples are rows and OTUs are columns
#create matrix with sample names (column labels)
SampleID<-matrix(colnames(dataset[,-c(51:57)]))
#object with OTU names as characters (R may assume these are numbers)
#OTU<-as.character(dataset$OTU_ID)
OTU<-(c("SampleID", as.character(dataset$OTU_ID)))
str(OTU)
dim(OTU)

#transpose dataset
#remove columns 51-57, which is the taxonomic information, this prevents the taxonomic text from confusing the dataset when transposed
dataset2<-as.data.frame(t(dataset[,-c(1,51:57)]))
#bind SampleID as first column so you can manipulate by SampleID as needed
dataset3<-cbind(SampleID[-1,], dataset2)
#Set object OTU as column names
colnames(dataset3)<-OTU

#preview transposed dataset3 to ensure it looks correct
dim(dataset3)
str(dataset3[1:10])
dataset3[1:10,1:10]

#export csv of transposed dataset
write.csv(dataset3, file="Wahl_ITS2_OTUtable_EMB.csv")
