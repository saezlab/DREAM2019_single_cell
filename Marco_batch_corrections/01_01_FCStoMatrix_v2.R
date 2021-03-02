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

pp_dir<-file.path(dir1,"plots","DREAM")
dir.create(pp_dir, showWarnings = FALSE)

dirFCS<- "~/Git/67-cell-lines-compensated"
frt <- readRDS(paste0(dirFCS, "/rds/data_asinhTransf_67-cell-lines-compensated.rds"))
dmt <- read.csv(file.path(paste0(dirFCS,"/excel/info-file_67-cell-lines.csv")),sep=";") # info file for the samples 
protinfo <- read.csv(file.path(paste0(dir1,"/excel"),"prot_info.csv"), sep=';',header = FALSE ) # mapping between metal tags and proteins
nms <- as.character(protinfo[,2])     # corresponding proteins

ind<-nrow(dmt)

fol_dream = paste0(dir1,"/rds/dream_rds/")
dir.create(fol_dream, showWarnings = FALSE)
prefix_dream = 'dream_'
suffix_rds = '.rds'

pb <- txtProgressBar(min = 1, max = length(ind), style = 3)
for(i in 1:length(ind)){
  DM <- as.data.frame(exprs(frt[[i]]))
  colnames(DM) <- nms
  DM$cell_line <-dmt$cell_line[ind[i]]
  DM$treatment <-dmt$treatment[ind[i]]
  DM$cellID<-c(1:nrow(DM))
  DM$fileID<-i
  DM$time <-dmt$time[ind[i]]
  DM<-melt(DM, id.vars = c("treatment", "cell_line","time","cellID","fileID"), variable.name = "channel", value.name = "value")
  saveRDS(DM, file =  file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)

rm(frt)
print("flowFrame deleted")

pb <- txtProgressBar(min = 1, max = 500, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"1",suffix_rds)))
for(i in 2:500){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID1.rds"))
close(pb)
print("Finished saving 1/10 :) ")
pb <- txtProgressBar(min = 501, max = 1000, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"501",suffix_rds)))
for(i in 502:1000){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID6.rds"))
close(pb)
print("Finished saving 2/10 :) ")

pb <- txtProgressBar(min = 1001, max = 1500, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"1001",suffix_rds)))
for(i in 1002:1500){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID2.rds"))
close(pb)
print("Finished saving 3/10 :) ")


pb <- txtProgressBar(min = 1501, max = 2000, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"1501",suffix_rds)))
for(i in 1502:2000){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID9.rds"))
close(pb)
print("Finished saving 4/10 :) ")

pb <- txtProgressBar(min = 2001, max = 2500, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"2001",suffix_rds)))
for(i in 2002:2500){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID3.rds"))
close(pb)
print("Finished saving 5/10 :) ")

pb <- txtProgressBar(min = 2501, max = 3000, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"2501",suffix_rds)))
for(i in 2502:3000){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID10.rds"))
close(pb)
print("Finished saving 6/10 :) ")

pb <- txtProgressBar(min = 3001, max = 3500, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"3001",suffix_rds)))
for(i in 3002:3500){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID4.rds"))
close(pb)
print("Finished saving 7/10 :) ")

pb <- txtProgressBar(min = 3501, max = 3700, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"3501",suffix_rds)))
for(i in 3502:3700){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID7.rds"))
close(pb)
print("Finished saving 8/10 :) ")

pb <- txtProgressBar(min = 3701, max = 4000, style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"3701",suffix_rds)))
for(i in 3702:4000){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID8.rds"))
close(pb)
print("Finished saving 9/10 :) ")

pb <- txtProgressBar(min = 4000, max = length(ind), style = 3)
plot_mat <- readRDS(file.path(fol_dream,paste0(prefix_dream,"4001",suffix_rds)))
for(i in 4002:length(ind)){
  DM <- readRDS(file.path(fol_dream,paste0(prefix_dream,i,suffix_rds)))
  plot_mat<- rbind(plot_mat,DM)
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
saveRDS(plot_mat, file.path(dir1,"rds","SingleCell_nocontrols_allChannels_fileID5.rds"))
close(pb)
print("Finished saving 10/10 :) ")

rm(list=ls())