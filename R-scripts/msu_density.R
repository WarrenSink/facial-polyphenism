setwd("/Users/Warren.Sink/Desktop/R projects/facial-polyphenism/CSVs")

#Needs phenotype data from 

library(tidyverse)
library(readr)

msu_cf = read_csv(file = "msu_cf_nonredundant.csv")
msu_cf_negctrl = read_csv(file = "msu_cf_redundant.csv")

msu_cf$target <- msu_cf$target %>% substr(., 1, 5) 
msu_cf$source <- msu_cf$source %>% substr(., 1, 5) 
msu_cf_negctrl$target <- msu_cf_negctrl$target %>% substr(., 1, 5) 
msu_cf_negctrl$source <- msu_cf_negctrl$source %>% substr(., 1, 5) 

msu_cf <- msu_cf[!grepl("ctrA1", msu_cf$target),]
msu_cf <- msu_cf[!grepl("ctrA2", msu_cf$target),]
msu_cf <- msu_cf[!grepl("ctrB1", msu_cf$target),]
msu_cf <- msu_cf[!grepl("ctrB2", msu_cf$target),]

msu_cf <- msu_cf[!grepl("ctrA1", msu_cf$source),]
msu_cf <- msu_cf[!grepl("ctrA2", msu_cf$source),]
msu_cf <- msu_cf[!grepl("ctrB1", msu_cf$source),]
msu_cf <- msu_cf[!grepl("ctrB2", msu_cf$source),]

msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrA1", msu_cf_negctrl$target),]
msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrA2", msu_cf_negctrl$target),]
msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrB1", msu_cf_negctrl$target),]
msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrB2", msu_cf_negctrl$target),]

msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrA1", msu_cf_negctrl$source),]
msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrA2", msu_cf_negctrl$source),]
msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrB1", msu_cf_negctrl$source),]
msu_cf_negctrl <- msu_cf_negctrl[!grepl("ctrB2", msu_cf_negctrl$source),]

#msu_cf[,2:3] <- lapply(msu_cf[,2:3], function(x) as.numeric(as.character(x)))
#msu_cf_negctrl[,2:3] <- lapply(msu_cf_negctrl[,2:3], function(x) as.numeric(as.character(x)))

#load("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu_dist-ratio_matrices.Rdata")
#msu_all_dist_pheno$msu_id <- as.character(msu_all_dist_pheno$msu_id)

msu_cf = msu_cf %>%
  select(target, source, Similarity) %>%
  left_join(., msu_all_dist_pheno[1:296,1:3], by = c("target" = "msu_id")) %>% 
  left_join(., msu_all_dist_pheno[1:296,1:3], by = c("source" = "msu_id")) %>%
  distinct()

id.List <- list()
for (i in 1:nrow(msu_cf)){
  if (msu_cf[i,4] == msu_cf[i,6]){
    id.List[[i]] = "co-twins"
  }
  else{
    id.List[[i]] = "positive_ctrl"
  }
}
ctrl_id <- do.call(rbind,id.List)
msu_cf <- cbind(ctrl_id, msu_cf)

msu_cf_negctrl = msu_cf_negctrl %>%
  select(target, source, Similarity) %>%
  left_join(., msu_all_dist_pheno[1:296,1:3], by = c("target" = "msu_id")) %>% 
  left_join(., msu_all_dist_pheno[1:296,1:3], by = c("source" = "msu_id")) %>%
  distinct()

neg.List =  list()
for (i in 1:nrow(msu_cf_negctrl)){
  neg.List[[i]] = "neg_ctrl"
}
ctrl_id <- do.call(rbind,neg.List)
msu_cf_negctrl <- cbind(ctrl_id, msu_cf_negctrl)
msu_cf <- rbind(msu_cf, msu_cf_negctrl)

msu_cf <- msu_cf %>% arrange(factor(ctrl_id, levels(ctrl_id)[c(1,2,3)]))


ggplot(subset(similarity_scores_pheno_aldi, status == 'negative_ctrl' |status == 'co-twins'| status == 'positive_ctrl'),
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.5),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 