# export  proteomics and transcriptomics

### Proteomics ====
proteomics_raw <- readRDS("./data/proteomics/MSstat_groupComparison_selceted.rds") %>% as_tibble()

proteomics <- proteomics_raw %>% 
	mutate(cell_line = gsub("normal_vs_","",Label)) %>%
	select(Protein,cell_line,Label,log2FC,issue) 
	
	
write_csv(proteomics,path = "./challenge_data/proteomics/proteomics.csv")	
	
### Transcriptomics ====
rnaseq_raw <- readRDS("./data/genomics/dat_RNAseq_Marcotte.rds") %>% as_tibble()

rnaseq <-rnaseq_raw %>% select(Cellline,everything())
 

write_csv(rnaseq,path = "./challenge_data/transcriptomics/rnaseq.csv")	



### Transcriptomics ====
cnv_raw <- readRDS("./data/genomics/CNV_MARCOTTE.rds") %>% as_tibble()




write_csv(cnv_raw,path = "./challenge_data/transcriptomics/CNV_MARCOTTE.csv")	
