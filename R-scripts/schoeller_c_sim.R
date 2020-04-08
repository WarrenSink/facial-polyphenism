##############script for joining AWS similarity results from that particuler df with dist and ratio df's########

#for test sample
all_dist_pheno_delta_73 <- all_dist_pheno_delta[1:73,]

all_dist_pheno_delta_73$Co_twins_ID <- as.integer(as.character(all_dist_pheno_delta_73$Co_twins_ID)) 

all_dist_pheno_delta_73 <- all_dist_pheno_delta_73 %>% arrange(Co_twins_ID)

schoeller_cf_cotwins <- schoeller_cf %>% filter(ctrl_id == "co-twins") 

schoeller_cf_cotwins$Co_twins_ID.x <- as.integer(as.character(schoeller_cf_cotwins$Co_twins_ID.x))

schoeller_cf_cotwins <- schoeller_cf_cotwins %>% arrange(Co_twins_ID.x)

all_dist_pheno_delta_73 <- cbind(all_dist_pheno_delta_73, schoeller_cf_cotwins[,4]) 

colnames(all_dist_pheno_delta_73)[361] <- "Similarity"

all_dist_pheno_delta_73 <- cbind(all_dist_pheno_delta_73[1:7],all_dist_pheno_delta_73[8:361])

#for neg ctrl
all_dist_pheno_delta_neg_ctrl <- all_dist_pheno_delta[74:78,]

schoeller_cf_neg_ctrl <- schoeller_cf %>% filter(ctrl_id == "neg_ctrl") 

all_dist_pheno_delta_neg_ctrl <- all_dist_pheno_delta_neg_ctrl %>% left_join(schoeller_cf_neg_ctrl[,3:4], by = c("Individual_2" = "source"))


#for pos ctrl
all_dist_pheno_delta_pos_ctrl <- all_dist_pheno_delta[79:83,]

schoeller_cf_pos_ctrl <- schoeller_cf %>% filter(ctrl_id == "positive_ctrl")

all_dist_pheno_delta_pos_ctrl <- all_dist_pheno_delta_pos_ctrl %>% 
  left_join(schoeller_cf_pos_ctrl[,2:4], by = c("Individual_1" = "target", "Individual_2" = "source")) %>%
  left_join(schoeller_cf_pos_ctrl[,2:4], by = c("Individual_1" = "source", "Individual_2" = "target")) %>%
  mutate(Similarity.z = ifelse(!is.na(Similarity.x), Similarity.x, Similarity.y)) %>% select(c(1:360, 363)) %>%
  rename("Similarity" = "Similarity.z")

#combine all 3 dfs together
all_dist_pheno_delta_sim <- rbind(all_dist_pheno_delta_73, all_dist_pheno_delta_neg_ctrl, all_dist_pheno_delta_pos_ctrl)

#general function
add_similarity <- function(feature.df, similarity.df, test.observations.range, feature.descriptors, end.column, ...){
  

  feature.df.2 <- feature.df[test.observations.range,]
  test.observations <- length(test.observations.range)
  feature.df.3 <- feature.df[((test.observations+1):(test.observations+5)),]
  feature.df.4 <- feature.df[((test.observations+6):(test.observations+10)),]
  
  feature.df.2$Co_twins_ID <- as.integer(as.character(feature.df.2$Co_twins_ID)) 
  feature.df.2 <- feature.df.2 %>% arrange(Co_twins_ID)
  similarity.df.cotwins <- similarity.df %>% filter(ctrl_id == "co-twins") 
  similarity.df.cotwins$Co_twins_ID.x <- as.integer(as.character(similarity.df.cotwins$Co_twins_ID.x))
  similarity.df.cotwins <- similarity.df.cotwins %>% arrange(Co_twins_ID.x)
  feature.df.2 <- cbind(feature.df.2, similarity.df.cotwins$Similarity) 
  colnames(feature.df.2)[(end.column+1)] <- "Similarity"
  #feature.df.2 <- cbind(feature.df.2[1:feature.descriptors],feature.df.2[(feature.descriptors+1):end.column])
  
  similarity.df.neg.ctrl <- similarity.df %>% filter(ctrl_id == "neg_ctrl") 
  feature.df.3 <- feature.df.3 %>% left_join(similarity.df.neg.ctrl[,3:4], by = c("Individual_2" = "source"))
  
  similarity.df.pos.ctrl <- similarity.df %>% filter(ctrl_id == "positive_ctrl")
  feature.df.4 <- feature.df.4 %>% 
    left_join(similarity.df.pos.ctrl[,2:4], by = c("Individual_1" = "target", "Individual_2" = "source")) %>%
    left_join(similarity.df.pos.ctrl[,2:4], by = c("Individual_1" = "source", "Individual_2" = "target")) %>%
    mutate(Similarity.z = ifelse(!is.na(Similarity.x), Similarity.x, Similarity.y)) %>% select(c(1:end.column, (end.column+3))) %>%
    rename("Similarity" = "Similarity.z")
  
  #combine all 3 dfs together
  rbind(feature.df.2, feature.df.3, feature.df.4)
  
}









