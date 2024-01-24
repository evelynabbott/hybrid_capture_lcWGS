#remove NAs
library(dplyr)
library(tidyverse)

rm(list=ls())
#setwd("/Users/evelynabbott/Dropbox/project/Projects/hybcap_mcav_clad/hybcap_december")
mat=read.table("clad_ngs.ibsMat")
bams <- read.table("bams.qc")[,1] # list of bam files to name matrix
dimnames(mat) <- list(bams,bams)

#print rows with NAs
nasums=data.frame(rowSums(is.na(mat)))
rownames(nasums)=bams
colnames(nasums)="nanum"
#print samples with highest to lowest NAs
#usually its a small number of samples which are generating the NAs in the other samples
nasums=arrange(nasums,desc(nanum))
removenames=rownames(nasums)[1:3]
mat = mat[!(rownames(mat)) %in% removenames, ]
mat = mat[,!(colnames(mat)) %in% removenames ]
#remove samples from bams list
bams <- read.table("bams.qc")
rownames(bams)=bams$V1
bams = bams[!(rownames(bams)) %in% removenames, ]

#did that fix it? make sure NAs are removed
nasums1=data.frame(rowSums(is.na(mat)))
rownames(nasums1)=bams
colnames(nasums1)="nanum"
nasums1=arrange(nasums1,desc(nanum)) #all 0

write.table(mat, file="clad_ngs.ibsMat_na_rem", row.names=FALSE, col.names=FALSE,sep = "\t")
write.table(bams, file="bams_ngs_na_rem", row.names=FALSE, col.names=FALSE,sep="\t")








