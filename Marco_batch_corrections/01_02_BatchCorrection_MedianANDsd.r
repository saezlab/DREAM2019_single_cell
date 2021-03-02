#Batch correction***********************************************************************************************************************************
library("RColorBrewer")
source("heatm2_pdf.r")
dir_compensated<- "~/Git/67-cell-lines-compensated"
dir1 <- "~/Git/Normalization"
setwd(dir1)

Mmi <- readRDS(file.path(paste0(dir_compensated,"/rds/median_DREAM.rds")))
Msd <- readRDS(file.path(paste0(dir_compensated,"/rds/sd_DREAM.rds")))
dmt <- read.csv(file.path(paste0(dir_compensated,"/excel"), "info-file_67-cell-lines.csv"), sep=";")

ind_date_corr <- which(!is.na(dmt$associated_cell_line))
for (i in ind_date_corr){
  ind_date_replace <- which(dmt$cell_line == dmt$associated_cell_line[i])
  dmt$date[i] <- dmt$date[ind_date_replace[1]]
}

dtv <- as.character(dmt[,8])
dtv <- as.Date(dtv, format="%d.%m.%Y")
dts <- sort(unique(dtv))

# Perform batch correction for median
colnames(Mmi) <- paste("col",1:ncol(Mmi),sep="_")
rmeans <- rowMeans(Mmi)            # initial (i.e. before batch correction) 
# means for the protein pairs
indd4 <- which(dtv %in% dts[4:16] )
M12 <- Mmi[,-indd4]
M22 <- Mmi[,indd4]
rmeans_b13 <- rowMeans(M12)
rmeans_b416 <- rowMeans(M22)

#for p-PLCg where measurements are missing
Mmi2<-Mmi[,-which(dtv %in% dts[4])]
rmeansPLC <- rowMeans(Mmi2)            # initial (i.e. before batch correction) 
indd4 <- which(dtv %in% dts[5:16] )
indd13 <- which(dtv %in% dts[1:3] )
M12 <- Mmi[,indd13]
M22 <- Mmi[,indd4]
rmeans_b13PLC <- rowMeans(M12)
rmeans_b416PLC <- rowMeans(M22)

rmeans_b13[30]<-rmeans_b13PLC[30]
rmeans_b416[30]<-rmeans_b416PLC[30]
rmeans[30]<-rmeansPLC[30]

saveRDS(rmeans_b13,file.path(paste0(dir1,"/rds"),"Batch1to3_median_DREAM.rds"))
saveRDS(rmeans_b416,file.path(paste0(dir1,"/rds"),"Batch4to16_median_DREAM.rds"))
saveRDS(rmeans,file.path(paste0(dir1,"/rds"),"Batch_median_DREAM.rds"))

# Perform batch correction for sd
colnames(Msd) <- paste("col",1:ncol(Msd),sep="_")
rmeans <- rowMeans(Msd)            # initial (i.e. before batch correction) 
# means for the protein pairs
indd4 <- which(dtv %in% dts[4:16] )
M12 <- Msd[,-indd4]
M22 <- Msd[,indd4]
rmeans_b13 <- rowMeans(M12)
rmeans_b416 <- rowMeans(M22)

#for p-PLCg where measurements are missing
Msd2<-Msd[,-which(dtv %in% dts[4])]
rmeansPLC <- rowMeans(Msd2)            # initial (i.e. before batch correction) 
indd4 <- which(dtv %in% dts[5:16] )
indd13 <- which(dtv %in% dts[1:3] )
M12 <- Msd[,indd13]
M22 <- Msd[,indd4]
rmeans_b13PLC <- rowMeans(M12)
rmeans_b416PLC <- rowMeans(M22)

rmeans_b13[30]<-rmeans_b13PLC[30]
rmeans_b416[30]<-rmeans_b416PLC[30]
rmeans[30]<-rmeansPLC[30]

saveRDS(rmeans_b13,file.path(paste0(dir1,"/rds"),"Batch1to3_sd_DREAM.rds"))
saveRDS(rmeans_b416,file.path(paste0(dir1,"/rds"),"Batch4to16_sd_DREAM.rds"))
saveRDS(rmeans,file.path(paste0(dir1,"/rds"),"Batch_sd_DREAM.rds"))

rm(list=ls())
