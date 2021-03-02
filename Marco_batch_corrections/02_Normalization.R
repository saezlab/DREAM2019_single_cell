#heatmap******************************************************************************************************
dir1<- "~/Git/Normalization"
setwd(dir1)
library(data.table)
library(plyr)
library(dplyr)
library(dtplyr)
library(ggplot2)
library(bbRtools)
library(RColorBrewer)
library(ggpubr)
library(utils)
library("flowCore")
source("heatm2_pdf.r")
theme_set(theme_classic())

dir_compensated<- "~/Git/67-cell-lines-compensated"
pp_dir<-file.path(dir1,"plots","DREAM")
dir.create(pp_dir, showWarnings = FALSE)

nms <- read.csv(file.path(paste0(dir1,"/excel"),"prot_info.csv"), sep=';',header = FALSE ) # mapping between metal tags and proteins
nms <- as.character(nms[,2])     # corresponding proteins
nms<-nms[c(18:19,21:36,38:56)]  # selecting "interesting" proteins

dmt <- read.csv(file.path(paste0(dir_compensated,"/excel/info-file_67-cell-lines_DREAM.csv")),sep=";") # info file for the samples 
dmt$fileID<-1:nrow(dmt)
ind_date_corr <- which(!is.na(dmt$associated_cell_line))
for (i in ind_date_corr){
  ind_date_replace <- which(dmt$cell_line == dmt$associated_cell_line[i])
  dmt$date[i] <- dmt$date[ind_date_replace[1]]
}
dmt[,8] <- as.character(dmt[,8])
dmt[,8] <- as.Date(dmt[,8], format="%d.%m.%Y")
dts <- sort(unique(dmt[,8]))

batch13 <-readRDS(file.path(paste0(dir1,"/rds"),"Batch1to3_median_DREAM.rds"))
batch416 <-readRDS(file.path(paste0(dir1,"/rds"),"Batch4to16_median_DREAM.rds"))
batch <-readRDS(file.path(paste0(dir1,"/rds"),"Batch_median_DREAM.rds"))
sd_batch13 <-readRDS(file.path(paste0(dir1,"/rds"),"Batch1to3_sd_DREAM.rds"))
sd_batch416 <-readRDS(file.path(paste0(dir1,"/rds"),"Batch4to16_sd_DREAM.rds"))
sd_batch <-readRDS(file.path(paste0(dir1,"/rds"),"Batch_sd_DREAM.rds"))

pb <- txtProgressBar(min = 1, max = 8, style = 3)
iii<- c(2:3,9:10)
for(i in iii){
DM <- readRDS(file.path(dir1,"rds",paste0("SingleCell_nocontrols_allChannels_fileID",i,".rds")))

#selecting interesting proteins********************************
DM<-DM[which(DM$channel %in% nms),]
DM$time[which(DM$treatment == "full")]<-0
#batch correction********************************
start_id<- DM$fileID[1]
end_id<- DM$fileID[length(DM$fileID)]
dmt_selected <-dmt[dmt$fileID%in%c(start_id:end_id),]
id_firstB <- dmt_selected$fileID[which(dmt_selected$date %in% dts[1:3])]
id_secondB <- dmt_selected$fileID[which(dmt_selected$date %in% dts[4:16])]

if(length(id_firstB)>0){
  for(cc in nms[-17]){
    indd_c<-which(names(batch13)==cc)
    indd<-which(DM$fileID %in% id_firstB&DM$channel==cc)
    DM$value[indd]<-(((DM$value[indd]-batch13[indd_c])/sd_batch13[indd_c])+batch[indd_c])*sd_batch[indd_c]
  }
}

if(length(id_secondB)>0){
  for(cc in nms[-17]){
    indd_c<-which(names(batch416)==cc)
    indd<-which(DM$fileID %in% id_secondB&DM$channel==cc)
    DM$value[indd]<-(((DM$value[indd]-batch416[indd_c])/sd_batch416[indd_c])+batch[indd_c])*sd_batch[indd_c]
  }
}

#data clean up*******************************
if(length(id_firstB)>0){
  DM$value[which(DM$fileID %in% id_firstB&DM$channel=="p-HER2")]<-NA
}
indd <- dmt_selected$fileID[which(dmt_selected$date==dts[4])]
if(length(indd)>0){
  DM$value[which(DM$fileID %in% indd&DM$channel=="p-PLCg2")]<-NA
}

indd <- dmt_selected$fileID[which(dmt_selected$treatment=="singlets")]
if(length(indd)>0){
  DM<-DM[-which(DM$fileID %in% indd),]
}
indd <- dmt_selected$fileID[which(dmt_selected$time=="0short")]
if(length(indd)>0){
  DM$time[which(DM$fileID %in% indd)]<-0
}
DM$treatment<-as.character(unlist(DM$treatment))
DM$cell_line<-as.character(unlist(DM$cell_line))
DM$time<-as.numeric(as.character(unlist(DM$time)))
DM$cellID<-as.numeric(as.character(unlist(DM$cellID)))
DM$fileID<-as.numeric(as.character(unlist(DM$fileID)))
DM$channel<-as.character(unlist(DM$channel))
DM$value<-as.numeric(as.character(unlist(DM$value)))
saveRDS(DM, file.path(dir1,"rds",paste0("SingleCell_fileID",i,".rds")))

#removing control samples**********************************
no_control <- dmt_selected$fileID[which(is.na(dmt_selected$associated_cell_line))]
DM <- DM[which(DM$fileID %in% no_control),]
saveRDS(DM, file.path(dir1,"rds",paste0("SingleCell_nocontrol_fileID",i,".rds")))

#print for quality control********************************
print(unique(DM$treatment))
print(unique(DM$cell_line))
print(unique(DM$time))
print(unique(DM$channel))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)

rm(list=ls())