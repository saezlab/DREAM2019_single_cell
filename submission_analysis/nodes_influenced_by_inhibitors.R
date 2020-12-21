

# checking median data
# 1. find strong differences between EGF~inhibited
# 1. MEK independent ERK signaling
library(tidyverse)
library(broom)
library(pheatmap)
library(ggplot2)
# median_data = read_csv("./challenge_data/median_phospho/median_phospho_data.csv")
# data_subset = median_data %>% filter(time == 0) %>%
#     gather(marker,value,-cell_line,-treatment,-time) 


median_data = read_rds("./data/median_data/interpolated_median_allsamples_correct_times.rds")

data_subset = median_data %>% 
    #filter(time == 0) %>%
    rename(marker=reporter)


### Which p-sites are significantly different between inhibited and EGF treaments


# Test the change between condition EGF and EGF+iXXX across paired timepoints
# join with itself to compare EGF with other treatments: 
# compute original scale( ie. without the asinh(x/5) scaling) to use log2FC

orig_scale_data = data_subset %>% mutate(value = sinh(5*value))

t_test_scale_results <- orig_scale_data %>%
    filter(treatment =="EGF") %>%
    rename(ref_value = value,ref_treatment= treatment) %>%
    inner_join(filter(orig_scale_data,! treatment %in% c("EGF","full")),by = c("cell_line", "time", "marker")) %>%
    group_by(marker,treatment) %>%
    mutate(log2FC = log2(value/ref_value)) %>%
    nest() %>%
    mutate(ttest  = map(data,function(df){
        tmp <- t.test(x = log2(df$value+1e-3),y=log2(df$ref_value+1e-3),paired = TRUE)
        tmp2 <- t.test(x=df$log2FC)
        ttmp <- tidy(tmp)
        ttmp$mean_inhib = mean(df$value,na.rm = TRUE)
        ttmp$mean_EGF = mean(df$ref_value,na.rm = TRUE)
        return(ttmp)
    })
    ) %>% unnest(ttest) %>%
    mutate(p.adj = p.adjust(p.value,method = "fdr")) %>%
    select(-data) 



t_test_scale_results <- t_test_scale_results %>%
    arrange(p.adj) %>% print(n=40)

t_test_scale_results %>% 
    rename(log2FC = estimate) %>%
    select(marker,treatment,log2FC,p.adj) %>% 
    ggplot() + geom_point(aes(marker,treatment,col=log2FC,size=-log10(p.adj))) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
          legend.position="bottom") + 
    scale_color_gradient2()

ggsave("./publication/figures/inhibitor_influenced_nodes_log2FC.pdf",width = 7, height=3)


t_test_scale_results %>% filter(treatment  =="iMEK",marker=="p.ERK")

# changes in markers used in SC1: 
t_test_scale_results %>% 
    rename(average_log2FC = estimate) %>%
    filter( treatment !="imTOR") %>%
    filter(marker %in% c("p.Akt.Ser473.","p.ERK","p.HER2","p.PLCg2","p.S6")) %>%
    arrange(average_log2FC) %>% select(marker, treatment,average_log2FC, p.adj)


### cell-lines per log2FC in iMEK

orig_scale_data %>%
    filter(treatment =="EGF") %>%
    rename(ref_value = value,ref_treatment= treatment) %>%
    full_join(filter(orig_scale_data,! treatment %in% c("EGF","full")),by = c("cell_line", "time", "marker")) %>%
    group_by(marker,treatment) %>%
    mutate(log2FC = log2(value/ref_value)) %>%
    filter(marker =="p.ERK",treatment =="iMEK") %>%
    ggplot() + geom_density(aes(log2FC))
    



### Vulcano of p-sites

t_test_scale_results %>% 
    rename(log2FC = estimate) %>%
    # calculate the difference on the original space
    mutate(legend = ifelse(p.adj<0.05 && abs(log2FC)>2,paste(marker,treatment,sep = "_"),"" )) %>%
    mutate(signif = ifelse(p.adj<0.05 && abs(log2FC)>1,TRUE,FALSE )) %>%
    ggplot(aes(log2FC,-log10(p.adj))) +
    geom_point(aes(col=signif)) + 
    geom_hline(yintercept = -log10(0.05),linetype = "dotted") +
    geom_vline(xintercept = -1,linetype = "dotted") +
    geom_vline(xintercept = 1,linetype = "dotted") +
    ggrepel::geom_text_repel(aes(label=legend)) +
    scale_color_manual(values = c("TRUE"="red","FALSE"="grey")) +
    theme_bw() + guides(color=FALSE)





### found strong changes: 
# pERK - iMEK
# p.GSK3B - iPKC
# p90RSK - iMEK
# p-S6 - iPKC
# p.S6 - iPI3K
# pERK - iEGFR
# pAKT.Ser.473 - iPI3K
# pGSK3B - iEGFR
# iDU - MEK


#### 
# 1.  MEK indepentend ERK signaling
data_subset_ERK_MEK = median_data %>% filter(time == 0) %>%
    filter(treatment %in% c("EGF","iMEK")) %>% 
    filter(!is.na(p.ERK)) %>%
    select(cell_line, treatment,time,p.ERK) %>%
    filter(!cell_line %in% c("184B5","HCC202","ZR751"))


# order the cell-lines based on difference between the treatments: 
activity_diff <- data_subset_ERK_MEK %>% spread(treatment,p.ERK) %>%
    mutate(reduced_activity = EGF-iMEK) %>% arrange(reduced_activity)

hist(activity_diff$reduced_activity,10)
activity_diff_cl <- activity_diff %>% pull(cell_line)


# show the difference in p-ERK for each cell-line:     
data_subset_ERK_MEK %>% 
    mutate(cell_line = factor(cell_line,levels =activity_diff_cl )) %>%
    
    ggplot() + geom_col(aes(cell_line,p.ERK,fill=treatment),position = "dodge") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))



#### 
# 2. p.GSK3B to iPKC
data_subset_GSK3b_PKC = median_data %>% filter(time == 0) %>%
    filter(treatment %in% c("EGF","iPKC")) %>% 
    filter(!is.na(p.GSK3b)) %>%
    select(cell_line, treatment,time,p.GSK3b)  %>%
     filter(!cell_line %in% c("HCC1428","HCC1806","Hs578T"))


# order the cell-lines based on difference between the treatments: 
activity_diff <- data_subset_GSK3b_PKC %>% spread(treatment,p.GSK3b) %>%
    mutate(reduced_activity = EGF-iPKC) %>% arrange(reduced_activity)

hist(activity_diff$reduced_activity,10)
activity_diff_cl <- activity_diff %>% pull(cell_line)


# show the difference in p-ERK for each cell-line:     
data_subset_GSK3b_PKC %>% 
    mutate(cell_line = factor(cell_line,levels =activity_diff_cl )) %>%
    
    ggplot() + geom_col(aes(cell_line,p.GSK3b,fill=treatment),position = "dodge") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust = 1))


