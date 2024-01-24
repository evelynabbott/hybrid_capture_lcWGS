library(dplyr)
library(tidyverse)
library(vegan)
rm(list=ls())
setwd("")
#start with PCA

#load the data - ibs matrix and list of bams ---------------------------
bams_tacc=read.table("bams.nr")
bams_tacc$V1=gsub('_S.*','',bams_tacc$V1) #rename to match metadata

ibsMat=read.table("clad_ngs.ibsMat")
dnames=bams_tacc$V1
dimnames(ibsMat) <- list(dnames,dnames)

#load metadata and remove clones and samples which did not pass the qualRanks filter
load("/Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december/metadata_sampswithmeta.Rdata") #load metadata including all samples
keeps=dnames
metadata=metadata[metadata$id %in% keeps,]

#run capscale --------------------------
pp0_ad=capscale(ibsMat~1)
plot(pp0_ad$CA$eig/sum(pp0_ad$CA$eig))
scores=pp0_ad$CA$u

lab=rownames(ibsMat)

#inital plot to identify outliers
axes2plot=c(1,2)
plot(scores[,axes2plot],pch=19,main="",sub="",cex=0.2)
#ordispider(scores[,axes2plot],lab,label=T,cex=0.5)
ordihull(scores[,axes2plot],lab,label=T,cex=0.8)

#remove outliers from matrix
removenames=c("124","110","120") #example samples
ibsMat = ibsMat[!(rownames(ibsMat)) %in% removenames, ] #remove from matrix
ibsMat = ibsMat[,!(colnames(ibsMat)) %in% removenames ]
metadata = metadata[!metadata$id %in% removenames, ] #remove from metadata

#run capscale again
pp0_ad=capscale(ibsMat~1)
plot(pp0_ad$CA$eig/sum(pp0_ad$CA$eig))
scores=pp0_ad$CA$u

lab=rownames(ibsMat)

#make sure outliers were removed
axes2plot=c(1,2)
plot(scores[,axes2plot],pch=19,main="",sub="",cex=0.2)
#ordispider(scores[,axes2plot],lab,label=T,cex=0.5)
ordihull(scores[,axes2plot],lab,label=T,cex=0.8)

#assign colors to metadata column of interest
metadata$habcol = ifelse(metadata$habitat == "Offshore","green",
                         ifelse(metadata$habitat == "Deep", "darkblue","magenta"))

col=metadata$clustcol
lab=metadata$habitat

axes2plot=c(1,2)
plot(scores[,axes2plot],col=col,pch=19,main="",sub="",cex=1)
#ordispider(scores[,axes2plot],lab,label=T,cex=0.5)
ordihull(scores[,axes2plot],lab,label=T,cex=0.8)

#ADMIXTURE ------------------------
rm(list=ls())
setwd("/Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december/ibs/minmaf_tests")
#load functions
source("/Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december/plot_admixture_v4_function.R")

npops <- 4 #choose number of populations
inName <- paste0('clad_k', npops, '.qopt') #inName corresponds to .qopt files to upload later
#pops
inds=read.table("bams.nr",sep="\t") #bams used
inds=data.frame(lapply(inds, gsub, pattern="_S.*",replacement = "")) #make file names match metadata
load("/Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december/metadata_sampswithmeta.Rdata") #load metadata including all samples

keeps=inds$V1
metadata=metadata[metadata$id %in% keeps,]

coldata=metadata #rename to match this script

pops=data.frame(cbind(coldata$id,coldata$habitat)) #make pops file
write.table(pops, file = "pops",sep = "\t",col.names = c("id","habitat")) #save as tab delimited
pops=read.table("pops") #read as tab delimited table

#---------------
npops=as.numeric(sub("\\D+(\\d+)\\..+","\\1",inName))
tbl=read.table(paste(inName,sep=""),header=F)

#remove outliers identified by PCA
#remove samples from tbl which are not in coldata.
#keeps=bams$qpot_order
#tbl=tbl[rownames(tbl) %in% keeps,]

##########
i2p=read.table("pops")
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
rownames(tbl) <- tbl$ind <- sub("(.*?)\\..*$", "\\1", tbl$ind)
head(tbl,20) # this is how the resulting dataset must look

# putting populations in desired order (edit pop names as needed or skip to plot them alphabetically)
#set factor levels for plotting
tbl$pop=factor(tbl$pop,levels=c("Nearshore","Offshore","Deep"))

plotAdmixture(data=tbl,npops=npops,grouping.method="distance")+
  scale_fill_manual(values=c("#F8766D","#A7055C","#FB761F","#FFC71C"))





