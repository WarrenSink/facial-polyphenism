setwd("/Users/Warren.Sink/Desktop/R projects/facial-polyphenism/CSVs")

library(tidyverse)
library(readr)

schoeller_cf = read_csv(file = "schoeller_cf_nonredundant.csv")
schoeller_cf_negctrl = read_csv(file = "schoeller_cf_redundant.csv")
schoeller_cf_negctrl_scans = read_csv(file = "schoeller_positive_ctrl.csv")

schoeller_cf_negctrl_scans[1:2,3] <- "Mitchell_Ned"
schoeller_cf_negctrl_scans[3:4,3] <- "Haick_Ava"
schoeller_cf_negctrl_scans <- schoeller_cf_negctrl_scans %>% 
  mutate(ctrl_id = ifelse(status == "positive_ctrl", "neg_ctrl", status)) %>%
  .[,c(21,2:4)] %>%
  left_join(., all_dist_pheno[1:116,c(1:5,9)], by = c("target" = "Names")) %>% 
  left_join(., all_dist_pheno[1:116,c(1:5,9)], by = c("source" = "Names")) %>%
  distinct()

schoeller_cf$target <- schoeller_cf$target %>% str_remove(pattern = ".jpg")
schoeller_cf$source <- schoeller_cf$source %>% str_remove(pattern = ".jpg") 
schoeller_cf_negctrl$target <- schoeller_cf_negctrl$target %>% str_remove(pattern = ".jpg") 
schoeller_cf_negctrl$source <- schoeller_cf_negctrl$source %>% str_remove(pattern = ".jpg") 

#all_dist_pheno$Names <- as.character(all_dist_pheno$Names)
schoeller_cf = schoeller_cf %>%
  select(target, source, Similarity) %>%
  left_join(., all_dist_pheno[1:116,c(1:5,9)], by = c("target" = "Names")) %>% 
  left_join(., all_dist_pheno[1:116,c(1:5,9)], by = c("source" = "Names")) %>%
  distinct()

id.List <- list()
for (i in 1:nrow(schoeller_cf)){
  if (schoeller_cf[i,5] == schoeller_cf[i,10]){
    id.List[[i]] = "co-twins"
  }
  else{
    id.List[[i]] = "positive_ctrl"
  }
}
ctrl_id <- do.call(rbind,id.List)
schoeller_cf <- cbind(ctrl_id, schoeller_cf)

schoeller_cf_negctrl = schoeller_cf_negctrl %>%
  select(target, source, Similarity) %>%
  left_join(., all_dist_pheno[1:116,c(1:5,9)], by = c("target" = "Names")) %>% 
  left_join(., all_dist_pheno[1:116,c(1:5,9)], by = c("source" = "Names")) %>%
  distinct()

neg.List =  list()
for (i in 1:nrow(schoeller_cf_negctrl)){
  neg.List[[i]] = "neg_ctrl"
}
ctrl_id <- do.call(rbind,neg.List)
schoeller_cf_negctrl <- cbind(ctrl_id, schoeller_cf_negctrl)
schoeller_cf <- rbind(schoeller_cf, schoeller_cf_negctrl)
schoeller_cf <- rbind(schoeller_cf, schoeller_cf_negctrl_scans)

schoeller_cf <- schoeller_cf %>% arrange(factor(ctrl_id, levels(ctrl_id)[c(1,2,3)]))

schoeller_cf_sex<-schoeller_cf[!(schoeller_cf$Sex.x=="M" & schoeller_cf$Sex.y=="F"),]
schoeller_cf_sex<-schoeller_cf_sex[!(schoeller_cf_sex$Sex.x=="F" & schoeller_cf_sex$Sex.y=="M"),]
schoeller_cf_sex_female<-schoeller_cf_sex[!(schoeller_cf_sex$Sex.x=="M"),]
schoeller_cf_sex_male<-schoeller_cf_sex[!(schoeller_cf_sex$Sex.x=="F"),]

schoeller_cf_sex_ethnicity<-schoeller_cf_sex[!(schoeller_cf_sex$Ethnicity.x=="caucasian" & schoeller_cf_sex$Ethnicity.y=="hispanic"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="hispanic" & schoeller_cf_sex_ethnicity$Ethnicity.y=="caucasian"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="caucasian" & schoeller_cf_sex_ethnicity$Ethnicity.y=="black"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="black" & schoeller_cf_sex_ethnicity$Ethnicity.y=="caucasian"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="caucasian" & schoeller_cf_sex_ethnicity$Ethnicity.y=="asian"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="asian" & schoeller_cf_sex_ethnicity$Ethnicity.y=="caucasian"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="black" & schoeller_cf_sex_ethnicity$Ethnicity.y=="asian"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="black" & schoeller_cf_sex_ethnicity$Ethnicity.y=="hispanic"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="asian" & schoeller_cf_sex_ethnicity$Ethnicity.y=="black"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="hispanic" & schoeller_cf_sex_ethnicity$Ethnicity.y=="black"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="hispanic" & schoeller_cf_sex_ethnicity$Ethnicity.y=="asian"),]
schoeller_cf_sex_ethnicity<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="asian" & schoeller_cf_sex_ethnicity$Ethnicity.y=="hispanic"),]

schoeller_cf_sex_ethnicity_black<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="caucasian" | schoeller_cf_sex_ethnicity$Ethnicity.x=="hispanic" | schoeller_cf_sex_ethnicity$Ethnicity.x=="asian"),]
schoeller_cf_sex_ethnicity_white<-schoeller_cf_sex_ethnicity[!(schoeller_cf_sex_ethnicity$Ethnicity.x=="black" | schoeller_cf_sex_ethnicity$Ethnicity.x=="hispanic" | schoeller_cf_sex_ethnicity$Ethnicity.x=="asian"),]

schoeller_cf_ethnicity<-schoeller_cf[!(schoeller_cf$Ethnicity.x=="caucasian" & schoeller_cf$Ethnicity.y=="hispanic"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="hispanic" & schoeller_cf_ethnicity$Ethnicity.y=="caucasian"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="caucasian" & schoeller_cf_ethnicity$Ethnicity.y=="black"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="black" & schoeller_cf_ethnicity$Ethnicity.y=="caucasian"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="caucasian" & schoeller_cf_ethnicity$Ethnicity.y=="asian"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="asian" & schoeller_cf_ethnicity$Ethnicity.y=="caucasian"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="black" & schoeller_cf_ethnicity$Ethnicity.y=="asian"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="black" & schoeller_cf_ethnicity$Ethnicity.y=="hispanic"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="asian" & schoeller_cf_ethnicity$Ethnicity.y=="black"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="hispanic" & schoeller_cf_ethnicity$Ethnicity.y=="black"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="hispanic" & schoeller_cf_ethnicity$Ethnicity.y=="asian"),]
schoeller_cf_ethnicity<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Ethnicity.x=="asian" & schoeller_cf_ethnicity$Ethnicity.y=="hispanic"),]

schoeller_cf_ethnicity_black <- schoeller_cf_ethnicity %>%
  filter(Ethnicity.x == "black")

schoeller_cf_ethnicity_white <- schoeller_cf_ethnicity %>%
  filter(Ethnicity.x == "caucasian")

schoeller_cf_sex_ethnicity_white_sex<-schoeller_cf_sex_ethnicity_white[!(schoeller_cf_sex_ethnicity_white$Sex.x=="M" & schoeller_cf_sex_ethnicity_white$Sex.y=="F"),]
schoeller_cf_sex_ethnicity_white_sex<-schoeller_cf_sex_ethnicity_white_sex[!(schoeller_cf_sex_ethnicity_white_sex$Sex.x=="F" & schoeller_cf_sex_ethnicity_white_sex$Sex.y=="M"),]
schoeller_cf_sex_ethnicity_white_sex_male<-schoeller_cf_sex_ethnicity_white_sex[!(schoeller_cf_sex_ethnicity_white_sex$Sex.x=="F"),]
schoeller_cf_sex_ethnicity_white_sex_female<-schoeller_cf_sex_ethnicity_white_sex[!(schoeller_cf_sex_ethnicity_white_sex$Sex.x=="M"),]

schoeller_cf_sex_ethnicity_black_sex<-schoeller_cf_sex_ethnicity_black[!(schoeller_cf_sex_ethnicity_black$Sex.x=="M" & schoeller_cf_sex_ethnicity_black$Sex.y=="F"),]
schoeller_cf_sex_ethnicity_black_sex<-schoeller_cf_sex_ethnicity_black_sex[!(schoeller_cf_sex_ethnicity_black_sex$Sex.x=="F" & schoeller_cf_sex_ethnicity_black$Sex.y=="M"),]
schoeller_cf_sex_ethnicity_black_sex_male<-schoeller_cf_sex_ethnicity_black_sex[!(schoeller_cf_sex_ethnicity_black$Sex.x=="F"),]
schoeller_cf_sex_ethnicity_black_sex_female<-schoeller_cf_sex_ethnicity_black_sex[!(schoeller_cf_sex_ethnicity_black$Sex.x=="M"),]

schoeller_cf_ethnicity_sex<-schoeller_cf_ethnicity[!(schoeller_cf_ethnicity$Sex.x=="M" & schoeller_cf_ethnicity$Sex.y=="F"),]
schoeller_cf_ethnicity_sex<-schoeller_cf_ethnicity_sex[!(schoeller_cf_ethnicity_sex$Sex.x=="F" & schoeller_cf_ethnicity_sex$Sex.y=="M"),]
schoeller_cf_ethnicity_sex_male<-schoeller_cf_ethnicity_sex[!(schoeller_cf_ethnicity_sex$Sex.x=="F"),]
schoeller_cf_ethnicity_sex_female<-schoeller_cf_ethnicity_sex[!(schoeller_cf_ethnicity_sex$Sex.x=="M"),]

library(ggplot2)

sample10_schoeller <- schoeller_cf %>% group_by(ctrl_id) %>% sample_n(size = 10)

similarity_median = sample10_schoeller %>%
  group_by(ctrl_id) %>%
  summarise(similarity.median = median(Similarity))

sample10_schoeller %>%
  ggplot( ., aes(x = Similarity, group = ctrl_id, fill = ctrl_id)) +
  geom_density(adjust=1.5, position="fill")

ggplot(sample10_schoeller, aes(x=Similarity, fill=ctrl_id)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity") +
  geom_density() +
  geom_vline(data=similarity.median, aes(xintercept=similarity.median,  colour=ctrl_id),
             linetype="dashed", size=1)

ggplot(subset(schoeller_cf, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#y = "Density"


ggplot(subset(schoeller_cf_ethnicity, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition w.r.t. Ethnicity", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#y = "Density"

ggplot(subset(schoeller_cf_ethnicity_black, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition w.r.t. African Americans", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#y = "Density"

ggplot(subset(schoeller_cf_ethnicity_white, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition w.r.t. Caucasians", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#y = "Density"

ggplot(subset(schoeller_cf_sex, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition w.r.t. Sex", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(schoeller_cf_sex_male, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition w.r.t. Males", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(subset(schoeller_cf_sex_female, ctrl_id == 'neg_ctrl' |ctrl_id == 'co-twins'| ctrl_id == 'positive_ctrl'),
       aes(x=Similarity, color=ctrl_id)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=ctrl_id), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition w.r.t. Females", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme(legend.position=c(0.5,0.9),
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
