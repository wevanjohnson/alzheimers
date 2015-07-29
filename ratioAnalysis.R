setwd("~/collaborations/graves_steven/alzheimers_13_04")

source("gravesFunctions.R")

#### Libraries
library(ROCR) 
library("heatmap.plus")


dat=read.table("adratios.txt", header=T,sep='\t')
samplename=dat[,1]
gender=dat[,2]
status=dat[,3]
date=dat[,4]
markers=dat[,-c(1:4)]


#### Heatmap ####
pdf("results/HeatmapRatio.pdf",width=20,height=15)
hmdat=t(as.matrix(markers))
heatmap.plus(hmdat,col=redgreen(50),ColSideColors=cbind('Disease'=as.character(as.numeric(status)),'Gender'=as.numeric(gender),"Batch"=as.numeric(date)),cexRow=.5,cexCol=.5,main="Alzheimer's Ratio Data")
dev.off()


#### Correlation Analysis  #####
corMat=cor(markers,use="pairwise.complete.obs")
write.table(round(corMat,3),file="results/correlation_matrix.xls",quote=F,sep='\t')

pdf("results/correlation.pdf",width=20,height=15)
heatmap.plus(corMat,cexRow=.5,cexCol=.5,main="Marker Correlations")
dev.off()

cordist=dist(corMat)
pdf(file="results/marker_clust.pdf",height=5,width=30)
plot(hclust(cordist),ann=F,yaxt='n')
dev.off()

##### Single Markers ####

# Gender
genderpreds = singmarktab(markers[-41,],as.factor(as.character(gender[-41])))
write.table(genderpreds,file="results/genderPredictions.xls",sep='\t',quote=F,row.names=F)

# Disease
ADsingle = singmarktab(markers,status,outdir="results/single/")
write.table(ADsingle,file="results/ADPredictions.xls",sep='\t',quote=F,row.names=F)

##### Multiple Markers ####
multmods = markerselect(markers, status)
markertab=markerfreq(multmods,marks=colnames(markers))

markertab
 
##### Selected Markers ####
cut = .1

keep=rownames(markertab)[markertab[,1]>=cut]
keep

markshort=markers[,match(keep,colnames(markers))]

#### Heatmap ####
pdf("results/HeatmapRatioShort.pdf",width=20,height=5)
hmdat=t(as.matrix(markshort))
heatmap.plus(hmdat,col=redgreen(50),ColSideColors=cbind('Disease'=as.character(as.numeric(status)),'Gender'=as.numeric(gender),"Batch"=as.numeric(date)),cexRow=.5,cexCol=.5,main="Alzheimer's Ratio Data")
dev.off()


#### Correlation Analysis  #####
corMat_short=cor(markshort,use="pairwise.complete.obs")
write.table(round(corMat_short,3),file="results/correlation_matrix_short.xls",quote=F,sep='\t')

pdf("results/correlation_short.pdf",width=10,height=10)
heatmap.plus(corMat_short,cexRow=.5,cexCol=.5,main="Marker Correlations")
dev.off()

cordist=dist(corMat_short)
pdf(file="results/marker_clust_short.pdf",height=5,width=15)
plot(hclust(cordist),ann=F,yaxt='n')
dev.off()

##### Single Markers ####

# Gender
genderpreds_short = singmarktab(markshort[-41,],as.factor(as.character(gender[-41])))
write.table(genderpreds_short,file="results/genderPredictions_short.xls",sep='\t',quote=F,row.names=F)

# Disease
ADsingle_short = singmarktab(markshort,status,outdir="results/single/")
ADsingle_short
write.table(ADsingle_short,file="results/ADPredictions_short.xls",sep='\t',quote=F,row.names=F)

##### Multiple Markers ####
multmods_short = markerselect(markshort, status)
markertab_short=markerfreq(multmods_short,marks=colnames(markshort))

markertab_short

### Get models
models= multmods_short$models
models[[14]]=c(1,2,4)
models[[15]]=c(1,3,4)

modelplots(models,mark_dat,status,outdir="finalmarkers/")
 
#Variables: 1 2 4 5 6 ; peptides: marker_531.25 marker_1568.18 marker_804.6 marker_602.30 marker_708.26 ; AUC: 0.9117647 
#Variables: 2 4 1 5 6 ; peptides: marker_1568.18 marker_804.6 marker_531.25 marker_602.30 marker_708.26 ; AUC: 0.9117647 
#Variables: 3 4 1 5 6 ; peptides: marker_1618.20 marker_804.6 marker_531.25 marker_602.30 marker_708.26 ; AUC: 0.9080882 
#Variables: 4 1 2 5 6 ; peptides: marker_804.6 marker_531.25 marker_1568.18 marker_602.30 marker_708.26 ; AUC: 0.9117647 
#Variables: 5 7 2 4 1 6 ; peptides: marker_602.30 marker_701.80 marker_1568.18 marker_804.6 marker_531.25 marker_708.26 ; AUC: 0.8970588 
#Variables: 6 11 2 4 1 5 ; peptides: marker_708.26 marker_874.60 marker_1568.18 marker_804.6 marker_531.25 marker_602.30 ; AUC: 0.890625 
#Variables: 7 2 12 10 13 ; peptides: marker_701.80 marker_1568.18 marker_892.37 marker_660.4 marker_583.3 ; AUC: 0.8549383 
#Variables: 8 2 10 11 ; peptides: marker_989.26 marker_1568.18 marker_660.4 marker_874.60 ; AUC: 0.7771261 
#Variables: 9 1 4 ; peptides: marker_810.615 marker_531.25 marker_804.6 ; AUC: 0.8235294 
#Variables: 10 2 1 4 ; peptides: marker_660.4 marker_1568.18 marker_531.25 marker_804.6 ; AUC: 0.8266254 
#Variables: 11 2 8 10 ; peptides: marker_874.60 marker_1568.18 marker_989.26 marker_660.4 ; AUC: 0.7771261 
#Variables: 12 7 2 10 13 ; peptides: marker_892.37 marker_701.80 marker_1568.18 marker_660.4 marker_583.3 ; AUC: 0.8549383 
#Variables: 13 1 4 2 6 ; peptides: marker_583.3 marker_531.25 marker_804.6 marker_1568.18 marker_708.26 ; AUC: 0.8588235  
 
 
#Variables: 1 2 4 5 6 ; peptides: marker_531.25 marker_1568.18 marker_804.6 marker_602.30 marker_708.26 ; AUC: 0.9117647 
#Variables: 3 4 1 5 6 ; peptides: marker_1618.20 marker_804.6 marker_531.25 marker_602.30 marker_708.26 ; AUC: 0.9080882 
#Variables: 1 2 4 ; peptides: marker_531.25 marker_1568.18 marker_804.6 ; AUC: 0.8421053
#Variables: 1 3 4 ; peptides: marker_531.25 marker_1618.20 marker_804.6 ; AUC: 0.8173375
#Variables: 7 2 12 10 13 ; peptides: marker_701.80 marker_1568.18 marker_892.37 marker_660.4 marker_583.3 ; AUC: 0.8549383 
#Variables: 9 1 4 ; peptides: marker_810.615 marker_531.25 marker_804.6 ; AUC: 0.8235294 


modelplot(c(7,2,12,10,13), models,mark_dat,status)
 
save.image("ratio.RData")


