dir1<- "~/Git/67-cell-lines-compensated"
analysis <-"_compensated"
setwd(dir1)
library("flowCore")
library("gplots")
library("RColorBrewer")
source("heatm2_pdf.r")


frt <- readRDS(paste0(dir1, "/rds/data_asinhTransf_67-cell-lines-compensated.rds"))
protinfo <- read.csv(file.path(paste0(dir1,"/excel"),"prot_info.csv"), sep=';',header = FALSE ) # mapping between metal tags and proteins
nms <- as.character(protinfo[,2])     # corresponding proteins


indip <- c(18:19,21:36,38:56)    # indices for "interesting" proteins


medintcomp <- function(fr,ind){
 #Takes a flow frame fr. Extracts the data matrix with 
 #the expression levels, and reduces it to the interesting 
 #proteins, given by the vector of indices ind. 
 #Computes the median intensities of these proteins across 
 #the cells and returns them as a vector.

DM <- exprs(fr)                       # data matrix with expression levels
DMr <- DM[,ind]                       # reduce matrix to the "interesting" proteins
miv <- apply(DMr,2,sd)
return(miv)
}                                     # end function medintcomp



mil <- lapply(frt,function(x){medintcomp(x,indip)}) # list of vectors with median intensities
                                      # for all experiments (samples)
Mmi <- do.call(cbind,mil)             # convert the list into a matrix

rownames(Mmi) <- nms[indip]

saveRDS(Mmi,file.path(paste0(dir1,"/rds"),"sd_DREAM.rds"))