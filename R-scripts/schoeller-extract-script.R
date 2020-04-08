setwd("/Users/Warren.Sink/github/facial-polyphenism/CSVs")
options(stringsAsFactors = FALSE)
library(tidyverse)
####### Phenotypic data######

#nonpaired

Co_twins_ID <- c('1','1','1','2','2','3','3','4','4','5','6','7','8','9','10','10','11','11','8','12','12','13','14','14','15','15','7','16','16','17','17','18','18','18','18','19',
                 '20','20','21','21','22','22','6','23','23','23','9','24','25','25','26','5','27','27','28','28','29','30','30','31','31','32','32','33','33','34','34',
                 '13','35','35','36','36','37','37','38','38','39','39','40','40','26','29','41','42','42','43','43','44','44','41','45','45','46','46','47','24','24','8','19','24','48',
                 '48','49','49','8','47','50','50','51','51',
                 '52','52','53','53','54','54',
                 #'52_broad','52_broad','53_broad','53_broad',
                 '55','55','56','56','57','57','58','58','59','59',
                 '60','60','61','61','62','62','63','63','64','64')

Individual_ID <- c('1A','1B','1C','2A','2B','3A','3B','4A','4B','5A','6A','7A','8A','9A','10A','10B','11A','11B','8B','12A','12B','13A','14A','14B','15A','15B','7B','16A','16B','17A','17B','18A','18B','18C','18D','19A',
                   '20A','20B','21A','21B','22A','22B','6B','23A','23B','23C','9B','24A','25A','25B','26A','5B','27A','27B','28A','28B','29A','30A','30B','31A','31B','32A','32B','33A','33B','34A','34B',
                   '13B','35A','35B','36A','36B','37A','37B','38A','38B','39A','39B','40A','40B','26B','29B','41A','42A','42B','43A','43B','44A','44B','41B','45A','45B','46A','46B','47A','24B','24C','8C','19B','24D','48A',
                   '48B','49A','49B','8D','47B','50A','50B','51A','51B',
                   '52A','52B','53A','53B','54A','54B',
                   #'52_hdA','52_hdB','53_hdA','53_hdB',
                   '55A','55B','56A','56B','57A','57B','58A','58B','59A','59B',
                   '60A','60B','61A','61B','62A','62B','63A','63B','64A','64B')

Names <- c("Allooh_Jonas","Allooh_Moses","Allooh_Noah","Andrews_Doug","Andrews_Ross",
           "Ayer_Carly","Ayer_Lily","Bowen_David","Bowen_Richard","Bowers_Brenda",
           "Braswell_Orlean","Bridges_Carolyn","Burns_Catherine","Burns_McKenzie","Catalano_Joe",
           "Catalano_Sal","Cepeda_Eurides","Cepeda_Ramon","Claucherty_Marion",
           "Croll-Baehre_Emma","Croll-Baehre_Marta","Davis_June","Demonet_David_2","Demonet_Larry_2",
           "Dwomoh-Piper_Chantelle","Dwomoh-Piper_Danielle","Elder_Marilyn","Estrada_Antonion",
           "Estrada_Jesus","Everingham_Douglas","Everingham_Steve","Furtick_Jason","Furtick_Keith",
           "Furtick_Kevin","Furtick_Victor","Garzilli_Denise","Ginley_Adalynn","Ginley_Breanna",
           "Goldberg_Adam","Goldberg_Peter","Haick_Ava","Haick_Isabel","Hudson_Irene","Humphery_Emma",
           "Humphery_Megan","Humphery_Sarah","Jackson_Sidney","Jacques_Mary","Jensen_Dane","Jensen_Jordan",
           "Kerner_Sharon","Key_Aidan","Kobylski_Cecelia","Kobylski_Elizabeth","Lipfird_Garry","Lipfird_Larry",
           "Maassen_Jamie","Martin_Kate","Martin_emily","McClure_Bryan","McClure_Bryce","McLeod_Keisher",
           "McLeod_Teisher","Mielke_Alyssa","Mielke_Emily","Mitchell_Fred","Mitchell_Ned","Moyer_Jean","Mueller_Kirk_2",
           "Mueller_Nate_2","Nagel_Jeff","Nagel_Steve","Nelson_Audri","Nelson_Erin","Nick_Skyler",
           "Nick_Spencer","Oliver_David","Oliver_Walter","Parks_Katie","Parks_Sarah","Prijatel","Qualkinbush_Jodie","Rabi_Kate",
           "Ramsdell_Cole","Ramsdell_Seth","Ream_Alyssa","Ream_Rachel","Rowan_Hannah","Rowan_Kaci",
           "Scott-Wolf_Adriana","Serra_Ivory","Serra_Shelter","Smith_Connor","Smith_Kyle","Starr_Ellen",
           "Steeves_Carrie","Steeves_Patty","Stewart_Martha","Stillwell_Diane","Thornhill_Jennie","Timotiwu_Jason",
           "Timotiwu_Jordan","Tynski_Daniel","Tynski_Kristen","Ullman_Helen","Underwood_Helen","Whited_Jackie",
           "Whited_Jessica","Widerka_Amanda","Widerka_Beth",
           "Darren_1", "Darren_2", "Karen_1", "Karen_2", "Lois_1", "Lois_2",
           #"Darren_2","Darren_2_broad","Karen_1","Karen_1_broad",
           "Elder_Marilyn(1)","Elder_Marilyn", "Everingham_Steve(1)","Everingham_Steve","Parks_Katie(1)","Parks_Katie", "Ramsdell_Cole(1)","Ramsdell_Cole","Starr_Ellen(1)","Starr_Ellen",
           "Andrews_Doug","Estrada_Antonion","McClure_Bryan","Mitchell_Fred","Nagel_Jeff", "Oliver_David", "Rowan_Hannah", "Starr_Ellen","Timotiwu_Jordan", 'Tynski_Daniel')

Sex <- c('M','M','M','M','M','F','F','M','M','F','F','F','F','F','M','M','M','M','F','F','F','F','M','M','F','F','F','M','M','M','M','M','M','M','M','F',
         'F','F','M','M','F','F','F','F','F','F','F','F','M','M','F','M','F','F','M','M','F','F','F','M','M','F','F','F','F','M','M',
         'F','M','M','M','M','F','F','M','M','M','M','F','F','F','F','F','M','M','F','F','F','F','F','M','M','M','M','F','F','F','F','F','F','M',
         'M','M','F','F','F','F','F','F','F',
         'M','M','F','F','F','F',
         #'M','M','F','F',
         'F','F','M','M','F','F','M','M','F','F',
         'M','M','M',"M",'M','M','F','F','M','M')


Ethnicity <- c('caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','caucasian','caucasian','hispanic','hispanic','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','black','caucasian','caucasian','caucasian','caucasian','caucasian','black','black','black','black','caucasian',
               'black','black','caucasian','caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','black','caucasian','caucasian','black','black',
               'caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','asian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','asian','caucasian','caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','caucasian','caucasian','asian',
               'asian','caucasian','caucasian','caucasian','black','caucasian','caucasian','hispanic','hispanic',
               'caucasian','caucasian','caucasian','caucasian','caucasian','caucasian',
               #'caucasian','caucasian','caucasian','caucasian',
               'caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','caucasian','black','black',
               'caucasian','caucasian','caucasian','black','caucasian','caucasian','caucasian','black','asian','caucasian')

Eye_color <- c('brown','brown','brown','blue','blue','blue','blue','brown','brown','brown','brown','brown','green','brown','brown','brown','brown','brown','green','blue','blue','blue','blue','blue','brown','brown','brown','brown','brown','brown','brown','brown','brown','brown','brown','blue',
               'brown','brown','blue','blue','blue','blue','brown','blue','blue','blue','brown','blue','blue','blue','blue','brown','blue','blue','blue','blue','blue','blue','blue','blue','blue','brown','brown','blue','blue','brown','brown',
               'blue','brown','brown','brown','brown','brown','brown','blue','blue','blue','blue','blue','blue','blue','blue','brown','brown','brown','green','green','blue','blue','brown','brown','brown','green','green','brown','blue','blue','green','blue','blue','brown',
               'brown','brown','brown','brown','brown','blue','blue','brown','brown',
               'green','green','blue','blue','blue','blue',
               #'green','green','blue','blue',
               'brown','brown','brown','brown','blue','blue','brown','brown','brown','brown',
               'blue','brown','blue','brown','brown','blue','blue','brown','brown','brown')

Hair_color <- c('black','black','black','bald','bald','blond','blond','grey','grey','brown','brown','brown','red','black','black','black','black','black','black','red','red','grey','grey','grey','black','black','brown','black','black','grey','grey','black','black','black','black','brown',
                'black','black','bald','bald','brown','brown','brown','blond','blond','blond','blond','blond','blond','blond','blond','grey','brown','brown','grey','grey','brown','blond','blond','blond','blond','black','black','blond','blond','black','black',
                'grey','brown','brown','brown','brown','brown','brown','brown','brown','grey','grey','blond','blond','blond','brown','brown','brown','brown','red','red','blond','blond','brown','black','black','brown','brown','grey','blond','blond','black','brown','blond','black',
                'black','brown','brown','black','grey','blond','blond','brown','brown',
                'brown','brown','blond','blond','blond','blond',
                #'brown','brown','blond','blond',
                'brown','brown','grey','grey','blond','blond','brown','brown','grey','grey',
                'bald','black','blond','black','brown','grey','blond','grey','black','brown')

Birth_year <- c('1985','1985','1985','1965','1965','2006','2006','1957','1957','1964','1945','1944','1958','1991','1943','1943','1977','1977','1958','1996','1996','1939','1940','1940','1987','1987','1944','1988','1988','1947','1947','1990','1990','1990','1990','1953',
                '2005','2005','1970','1970','2007','2007','1945','1998','1998','1998','1991','1982','1986','1986','1949','1964','1996','1996','1954','1954','1981','2002','2002','2000','2000','1977','1977','2000','2000','1950','1950',
                '1939','1984','1984','1967','1967','1977','1977','1992','1992','1936','1936','2001','2001','1949','1981','1982','1988','1988','1997','1997','1992','1992','1982','1972','1972','2001','2001','1942','1982','1982','1958','1953','1982','1994',
                '1994','1984','1984','1958','1942','1991','1991','2005','2005',
                '1992','1992','1972','1972','1952','1952',
                #'1992','1992','1972','1972',
                '1944','1944','1947','1947','2001','2001','1988','1988','1942','1942',
                '1965','1988','2000','1950','1967','1936','1942','1942','1994','1984')


pheno_data <- data.frame(Names, Individual_ID, Co_twins_ID, Sex, Ethnicity, Birth_year, Eye_color, Hair_color)
pheno_data[,1:8] <- lapply(pheno_data[,1:8] , function(x) as.factor(as.character(x)))
#pheno_data <- pheno_data %>% arrange(Co_twins_ID)
pheno_data$Birth_year = as.numeric(as.character(pheno_data$Birth_year))
pheno_data <- pheno_data %>% mutate(., Age = 2012-pheno_data$Birth_year)

rm(list = ls()[!ls() %in% c("pheno_data")])

####### Generate distances and ratios of distances x pheno######

#distance x ratios of distances

rd <- read.csv(file = 'schoeller_extract.csv', header = TRUE, stringsAsFactors = FALSE)
#df1 <- read.csv(file = 'Gimp_features.csv', header = TRUE, stringsAsFactors = FALSE)
df3 <- read.csv(file = 'schoeller_ctrl_coord_nose_norm.csv', header = TRUE, stringsAsFactors = FALSE)
df2 <- read.csv(file = 'schoeller_negctrl_coord_nose_norm.csv', header = TRUE, stringsAsFactors = FALSE)
#df4 <- read.csv(file = 'barack_and_andrew.csv', header = TRUE, stringsAsFactors = FALSE)
#df1 <- slice(df1, 2:31)
df2 <- slice(df2, 2:31)
df3 <- slice(df3, 2:31)
#df <- slice(df4, 2:31)

#rd <- full_join(rd,df,by="X")
#rd <- full_join(rd,df1,by="X")
rd <- full_join(rd,df2,by="X")
rd <- full_join(rd,df3,by="X")
#rd <- full_join(rd,df4,by="X")

rd <- rd[-1,]

vect<-function(n){
  a=strsplit(strsplit(n,split=",")[[1]][1],split="[[]")[[1]][2]
  b=strsplit(strsplit(n,split=",")[[1]][2],split="[]]")[[1]][1]
  return(c(as.numeric(a),as.numeric(b)))
}

eud<-function(x1,y1,x2,y2){
  sqrt((x2-x1)**2+(y2-y1)**2)
}

L=list()
for (i in 1:(nrow(rd)-1)){
  for (j in (i+1):nrow(rd)){
    
    n1=rd[i,1]; n2=rd[j,1] 
    r=as.numeric()
    
    for (v in 2:dim(rd)[2]){
      im.att1=vect(rd[i,v]); im.att2=vect(rd[j,v]) 
      d1= eud(im.att1[1],im.att1[2],im.att2[1],im.att2[2])
      r=c(r,d1)
    }
    
    name=paste(n1,n2,sep="_")
    L[[name]]=r
    #cat(name,"\n")
  }
}

M=do.call(rbind,L)
write.csv(M,"shbang_features_distance.csv")

pairwise_ratio <- function(z){
  tot=choose(length(z),2)
  count=1
  ratio_feature =as.numeric(rep(0,tot))
  
  for (k in 1:(length(z)-1)){
    for (l in (k+1):length(z)){
      
      #ratio_feature = c(ratio_feature, z[k] / z[l])
      ratio_feature[count]=z[k]/z[l]
      
      count=count+1
      if (count%%5000==0){
        print(paste(round((count/tot)*100,1),"%",sep=""))
      }
      
    }
  }
  return(ratio_feature)
}

ratio_names=as.character(rep(0,choose(dim(M)[1],2)))
counter=0
for (k in 1:(length(M[,1])-1)){
  for (l in (k+1):length(M[,1])){
    n3=names(M[k,1]); n4=names(M[l,1]) 
    counter=counter+1
    name1=paste(n3,n4,sep="__")
    ratio_names[counter]=name1
  }
}

M.ratio=matrix(nrow=choose(dim(M)[1],2),ncol=dim(M)[2])
row.names(M.ratio)<-ratio_names

for (j in 1:dim(M)[2]){
  columnz = M[,j]
  column=pairwise_ratio(columnz)
  M.ratio[,j]=column
  cat(paste("column ",j," completed",sep=""))}

#bind pheno and dist/ratio tables

M_t <- t(M)
M.ratio_t <- t(M.ratio)
all_dist_pheno <- cbind(pheno_data,M_t[,1:351])
all_ratio_pheno <- cbind(pheno_data,M.ratio_t[,1:61425])

# arranging the data so greatest value comes first 

all_dist_pheno <- data.frame(all_dist_pheno)
all_ratio_pheno <- data.frame(all_ratio_pheno)

mkshft.dist <- all_dist_pheno[,c(3, 10:360)] %>% group_by(Co_twins_ID) %>% arrange(Co_twins_ID, ChinBottom_LeftEyeDown) 
mkshft.dist <- as.data.frame(mkshft.dist)
mkshft.dist <- left_join(mkshft.dist, all_dist_pheno[,1:10], by = c("Co_twins_ID","ChinBottom_LeftEyeDown"))
all_dist_pheno <- cbind(mkshft.dist[353:354], mkshft.dist[1], mkshft.dist[355:360], mkshft.dist[2:352])
#mkshft.dist$Co_twins_ID <- as.integer(as.character(mkshft.dist$Co_twins_ID)) 
#mkshft.dist <- mkshft.dist %>% arrange(Co_twins_ID)
#all_dist_pheno <- cbind(all_dist_pheno[1:2], mkshft.dist[1], all_dist_pheno[4:9], mkshft.dist[2:352])

mkshft.ratio <- all_ratio_pheno[,c(3, 10:61434)] %>% group_by(Co_twins_ID) %>% arrange_all() 
mkshft.ratio <- as.data.frame(mkshft.ratio)
mkshft.ratio <- left_join(mkshft.ratio, all_ratio_pheno[,1:10], by = c("Co_twins_ID","ChinBottom_LeftEyeDown__ChinBottom_LeftEyeLeft"))
all_ratio_pheno <- cbind(mkshft.ratio[61427:61428], mkshft.ratio[1], mkshft.ratio[61429:61434], mkshft.ratio[2:61426])

#mkshft.ratio <- all_ratio_pheno[,c(3, 10:61434)] %>% group_by(Co_twins_ID) %>% arrange_all() 
#mkshft.ratio <- as.data.frame(mkshft.ratio)
#mkshft.ratio$Co_twins_ID <- as.integer(as.character(mkshft.ratio$Co_twins_ID)) 
#mkshft.ratio <- mkshft.ratio %>% arrange(Co_twins_ID)
#all_ratio_pheno <- cbind(all_ratio_pheno[1:2], mkshft.ratio[1], all_ratio_pheno[4:9], mkshft.ratio[2:61426])


rm(list = ls()[!ls() %in% c("pheno_data","all_dist_pheno","all_ratio_pheno")])

#all_dist_pheno <- all_dist_pheno%>%group_by(Co_twins_ID)%>%arrange(desc(ChinBottom_EyeLeft))
#all_ratio_pheno <- all_ratio_pheno%>%group_by(Co_twins_ID)%>%arrange(desc(ChinBottom_EyeLeft__ChinBottom_EyeRight))

#####Delta/FC###########

#delta
CoIds= all_dist_pheno[,3]


UniqIds=unique(CoIds)

count=1
L.delta<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:(length(idx)-1)){
    for (w in (v+1):length(idx)){
      i=idx[v]
      j=idx[w]
      nam= paste0(all_dist_pheno[i,1],"__",all_dist_pheno[j,1])
      delta<-abs(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]]) - as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]]))
      L.delta[[nam]]<-delta
      #"/n" -> makes new line
      cat(nam,"\n")
    }
  }
}

all_dist_pheno_absdelta<-do.call(rbind,L.delta)
colnames(all_dist_pheno_absdelta)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]


#dist abs(delta)/fc

CoIds= all_dist_pheno[,3]


UniqIds=unique(CoIds)

count=1
L.ratio<-list()
L.delta<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:(length(idx)-1)){
    for (w in (v+1):length(idx)){
      i=idx[v]
      j=idx[w]
      
      nam= paste0(all_dist_pheno[i,1],"__",all_dist_pheno[j,1])
      ratio<-as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]]) /as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]])
      delta<-as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]] - as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]]))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_dist_pheno_fc<-do.call(rbind,L.ratio)
colnames(all_dist_pheno_fc)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]
all_dist_pheno_delta<-do.call(rbind,L.delta)
colnames(all_dist_pheno_delta)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]

#ratio delta

CoIds= all_ratio_pheno[,3]


UniqIds=unique(CoIds)

count=1

L.delta<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:(length(idx)-1)){
    for (w in (v+1):length(idx)){
      i=idx[v]
      j=idx[w]
      
      nam= paste0(all_ratio_pheno[i,1],"__",all_ratio_pheno[j,1])
      delta<-as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]) -as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]])
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_ratio_pheno_delta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_delta)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]

#ratio abs(delta)/fc

CoIds= all_ratio_pheno[,3]


UniqIds=unique(CoIds)

count=1
L.ratio<-list()
L.delta<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:(length(idx)-1)){
    for (w in (v+1):length(idx)){
      i=idx[v]
      j=idx[w]
      
      nam= paste0(all_ratio_pheno[i,1],"__",all_ratio_pheno[j,1])
      ratio<-as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]) /as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]])
      delta<-abs(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]) -as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]]))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}
all_ratio_pheno_fc<-do.call(rbind,L.ratio)
colnames(all_ratio_pheno_fc)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]
all_ratio_pheno_absdelta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_absdelta)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]


###############paired pheno x paired dfs########

paired.rownames <- row.names(all_dist_pheno_delta)

pheno_data_paired = pheno_data[,c(3,5:9)]
Co_twins_ID <- c('1','1','1','10','11','12','13','14','15','16','17','18','18','18','18','18','18','19',
                 '2','20','21','22','23','23','23','24','24','24','24','24','24','25','26','27','28','29',
                 '3','30','31','32','33','34','35','36','37','38','39',
                 '4','40','41','42','43','44','45','46','47','48','49',
                 '5','50','51','52','53','54','55','56','57','58','59',
                 '6','60','61','62','63','64','7','8','8','8','8','8','8','9')

# paired sexes: male = 0; female = 1; both = 2
Sex_Diff <- c("0","0","0","0","0","1","1","0","1",
              "0","0","0","0","0","0","0","0",'1',
              "0","1","0","1","1","1","1","1","1",
              "1","1","1","1","0","1","1","0","1",
              "1","1","0","1","1",'0','0','0','1',
              '0','0','0','1','1','0','1','1','1',
              '0','1','0','2','2','1','1','0','1',
              '1','1','0','1','0','1','1','0','0',
              '0','1','0','1','1','1','1','1','1',
              '1','1')

all_dist_pheno_fc <- data.frame(all_dist_pheno_fc)
all_dist_pheno_fc <- cbind(Co_twins_ID, Sex_Diff, all_dist_pheno_fc)
all_dist_pheno_fc <- full_join(all_dist_pheno_fc, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_fc <- distinct(all_dist_pheno_fc,ChinBottom_LeftPupil,.keep_all= TRUE)
row.names(all_dist_pheno_fc) <- paired.rownames

all_dist_pheno_delta<-data.frame(all_dist_pheno_delta)
all_dist_pheno_delta<- cbind(Co_twins_ID, Sex_Diff, all_dist_pheno_delta)
all_dist_pheno_delta<-left_join(all_dist_pheno_delta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_delta <- distinct(all_dist_pheno_delta,ChinBottom_LeftPupil,.keep_all= TRUE)
row.names(all_dist_pheno_delta) <- paired.rownames

all_dist_pheno_absdelta<-data.frame(all_dist_pheno_absdelta)
all_dist_pheno_absdelta<- cbind(Co_twins_ID, Sex_Diff, all_dist_pheno_absdelta)
all_dist_pheno_absdelta<-left_join(all_dist_pheno_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_absdelta <- distinct(all_dist_pheno_absdelta,ChinBottom_LeftPupil,.keep_all= TRUE)
row.names(all_dist_pheno_absdelta) <- paired.rownames

all_ratio_pheno_fc <- data.frame(all_ratio_pheno_fc)
all_ratio_pheno_fc <- cbind(Co_twins_ID, Sex_Diff, all_ratio_pheno_fc)
all_ratio_pheno_fc <- left_join(all_ratio_pheno_fc, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_fc <- distinct(all_ratio_pheno_fc,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)
row.names(all_ratio_pheno_fc) <- paired.rownames

all_ratio_pheno_delta<-data.frame(all_ratio_pheno_delta)
all_ratio_pheno_delta<- cbind(Co_twins_ID, Sex_Diff, all_ratio_pheno_delta)
all_ratio_pheno_delta<-left_join(all_ratio_pheno_delta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_delta <- distinct(all_ratio_pheno_delta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)
row.names(all_ratio_pheno_delta) <- paired.rownames

all_ratio_pheno_absdelta<-data.frame(all_ratio_pheno_absdelta)
all_ratio_pheno_absdelta<- cbind(Co_twins_ID, Sex_Diff, all_ratio_pheno_absdelta)
all_ratio_pheno_absdelta<-left_join(all_ratio_pheno_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_absdelta <- distinct(all_ratio_pheno_absdelta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)
row.names(all_ratio_pheno_absdelta) <- paired.rownames

#arrange pheno data to be at the front

all_dist_pheno_absdelta <- all_dist_pheno_absdelta[,c(1:2,354:358,3:353)]
all_dist_pheno_delta <- all_dist_pheno_delta[,c(1:2,354:358,3:353)]
all_dist_pheno_fc <- all_dist_pheno_fc[,c(1:2,354:358,3:353)]
all_ratio_pheno_absdelta <- all_ratio_pheno_absdelta[,c(1:2,61428:61432,3:61427)]
all_ratio_pheno_delta <- all_ratio_pheno_delta[,c(1:2,61428:61432,3:61427)]
all_ratio_pheno_fc <- all_ratio_pheno_fc[,c(1:2,61428:61432,3:61427)]

#bring rownames to columns

all_dist_pheno_delta <- all_dist_pheno_delta %>%
  rownames_to_column(var = "Paired_Individuals") %>%
  separate(Paired_Individuals, c('Individual_1', 'Individual_2'), sep = '__') 

all_dist_pheno_absdelta <- all_dist_pheno_absdelta %>%
  rownames_to_column(var = "Paired_Individuals") %>%
  separate(Paired_Individuals, c('Individual_1', 'Individual_2'), sep = '__')

all_dist_pheno_fc <- all_dist_pheno_fc %>%
  rownames_to_column(var = "Paired_Individuals") %>%
  separate(Paired_Individuals, c('Individual_1', 'Individual_2'), sep = '__')

all_ratio_pheno_absdelta <- all_ratio_pheno_absdelta %>%
  rownames_to_column(var = "Paired_Individuals") %>%
  separate(Paired_Individuals, c('Individual_1', 'Individual_2'), sep = '__')

all_ratio_pheno_delta <- all_ratio_pheno_delta %>%
  rownames_to_column(var = "Paired_Individuals") %>%
  separate(Paired_Individuals, c('Individual_1', 'Individual_2'), sep = '__')

all_ratio_pheno_fc <- all_ratio_pheno_fc %>%
  rownames_to_column(var = "Paired_Individuals") %>%
  separate(Paired_Individuals, c('Individual_1', 'Individual_2'), sep = '__')

#arrange data wrt cotwins

all_dist_pheno_absdelta <- all_dist_pheno_absdelta %>% arrange(as.integer(Co_twins_ID)) 
all_dist_pheno_delta <- all_dist_pheno_delta %>% arrange(as.integer(Co_twins_ID)) 
all_dist_pheno_fc <- all_dist_pheno_fc %>% arrange(as.integer(Co_twins_ID)) 
all_ratio_pheno_absdelta <- all_ratio_pheno_absdelta %>% arrange(as.integer(Co_twins_ID)) 
all_ratio_pheno_delta <- all_ratio_pheno_delta %>% arrange(as.integer(Co_twins_ID)) 
all_ratio_pheno_fc <- all_ratio_pheno_fc %>% arrange(as.integer(Co_twins_ID))

#pair each observation (ie cotwin pair) to their similarity score from AWS CompareFaces 
#bring rownames to a column and split character vector in two so as to match with CompareFaces df
##the general function##
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

all_dist_pheno_delta_sim <- add_similarity(feature.df = all_dist_pheno_delta, 
                                           similarity.df = schoeller_cf, 
                                           test.observations.range = c(1:73), 
                                           feature.descriptors = 7, 
                                           end.column = 360,)

all_dist_pheno_fc_sim <- add_similarity(feature.df = all_dist_pheno_fc, 
                                        similarity.df = schoeller_cf, 
                                        test.observations.range = c(1:73), 
                                        feature.descriptors = 7, 
                                        end.column = 360,)

all_dist_pheno_absdelta_sim <- add_similarity(feature.df = all_dist_pheno_absdelta, 
                                              similarity.df = schoeller_cf, 
                                              test.observations.range = c(1:73), 
                                              feature.descriptors = 7, 
                                              end.column = 360,)

all_ratio_pheno_delta_sim <- add_similarity(feature.df = all_ratio_pheno_delta, 
                                            similarity.df = schoeller_cf, 
                                            test.observations.range = c(1:73), 
                                            feature.descriptors = 7, 
                                            end.column = 61434,)

all_ratio_pheno_fc_sim <- add_similarity(feature.df = all_ratio_pheno_fc, 
                                         similarity.df = schoeller_cf, 
                                         test.observations.range = c(1:73), 
                                         feature.descriptors = 7, 
                                         end.column = 61434,)

all_ratio_pheno_absdelta_sim <- add_similarity(feature.df = all_ratio_pheno_absdelta, 
                                               similarity.df = schoeller_cf, 
                                               test.observations.range = c(1:73), 
                                               feature.descriptors = 7, 
                                               end.column = 61434,)

all_dist_pheno_delta$Co_twins_ID <- as.character(all_dist_pheno_delta$Co_twins_ID)
all_dist_pheno_fc$Co_twins_ID <- as.character(all_dist_pheno_fc$Co_twins_ID)
all_dist_pheno_absdelta$Co_twins_ID <- as.character(all_dist_pheno_absdelta$Co_twins_ID)
all_dist_pheno_delta_73$Co_twins_ID <- as.character(all_dist_pheno_delta_73$Co_twins_ID)
all_dist_pheno_fc_73$Co_twins_ID <- as.character(all_dist_pheno_delta_73$Co_twins_ID)
all_dist_pheno_absdelta_73$Co_twins_ID <- as.character(all_dist_pheno_delta_73$Co_twins_ID)
all_ratio_pheno_delta$Co_twins_ID <- as.character(all_ratio_pheno_delta$Co_twins_ID)
all_ratio_pheno_fc$Co_twins_ID <- as.character(all_ratio_pheno_fc$Co_twins_ID)
all_ratio_pheno_absdelta$Co_twins_ID <- as.character(all_ratio_pheno_absdelta$Co_twins_ID)
all_ratio_pheno_delta_73$Co_twins_ID <- as.character(all_ratio_pheno_delta_73$Co_twins_ID)
all_ratio_pheno_fc_73$Co_twins_ID <- as.character(all_ratio_pheno_delta_73$Co_twins_ID)
all_ratio_pheno_absdelta_73$Co_twins_ID <- as.character(all_ratio_pheno_delta_73$Co_twins_ID)

rm(list = ls()[!ls() %in% c("pheno_data","all_dist_pheno","all_ratio_pheno",'all_dist_pheno_delta','all_dist_pheno_absdelta',
                            'all_dist_pheno_fc','all_ratio_pheno_delta','all_ratio_pheno_absdelta','all_ratio_pheno_fc',
                            "all_dist_pheno_log2","all_dist_pheno_abslog2", "all_dist_pheno_log2_fc",'all_dist_pheno_log2_absdelta',
                            'all_dist_pheno_log2_delta','all_dist_pheno_abslog2_fc','all_dist_pheno_abslog2_delta','all_dist_pheno_abslog2_delta',
                            "all_ratio_pheno_log2","all_ratio_pheno_abslog2",'all_ratio_pheno_log2_fc',
                            'all_ratio_pheno_log2_delta','all_ratio_pheno_log2_absdelta',"all_ratio_pheno_abslog2_fc","all_ratio_pheno_abslog2_delta",
                            'all_ratio_pheno_abslog2_absdelta','all_dist_pheno_delta_sim','all_dist_pheno_fc_sim',
                            'all_dist_pheno_absdelta_sim','all_ratio_pheno_delta_sim','all_ratio_pheno_fc_sim',
                            'all_ratio_pheno_absdelta_sim','schoeller_cf')])

######        PCA########

library("FactoMineR")
library("factoextra")
library(ggrepel)
library(ggpmisc)
library(mclust)

#dist

all_dist_pheno$Co_twins_ID <- as.integer(as.character(all_dist_pheno$Co_twins_ID))

all_dist_pheno %<>% arrange(Co_twins_ID)

scaled_all_dist_pheno<-scale(all_dist_pheno[,10:360])
scaled_all_dist_pheno<-cbind(all_dist_pheno[c(1:9)],scaled_all_dist_pheno)
res.pca <- prcomp(scaled_all_dist_pheno[,10:360])

#controls
#pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-PCA/Schoeller-dist-fc/schoeller_pca_dist_fc_controls_arranged.pdf", w=4, h=4)
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = TRUE, geom = "text",
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             col.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             title = "PCA of All Cotwins' Distance between Landmarks") +
  theme(legend.position = "none",
        plot.title = element_text(size=30)) + 
  geom_text(aes(label = all_dist_pheno$Co_twins_ID))
#dev.off()

#sex
groups <- as.factor(all_dist_pheno$Sex) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             col.ind = groups, fill.ind = groups,
             #addEllipses = TRUE,
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             legend.title = "Sex", #theme(legend.position = "none") + 
             title = "PCA of All Cotwins' Distance between Landmarks") +
  theme(plot.title = element_text(size=30))
#theme(legend.position = "none") +   
#geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#age
all_dist_pheno_age <- all_dist_pheno %>%
  mutate(FortyPlus = ifelse(Age > 40, 1, 0)) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             #col.ind = as.character(all_dist_pheno_fc_age$FortyPlus), 
             fill.ind = as.character(all_dist_pheno_age$FortyPlus),
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of All Cotwins' Distance between Landmarks") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  #scale_fill_discrete(name= "",labels = c("cowtins",'negative control','positive control'))+
  theme(plot.title = element_text(size=30),
        legend.position = "none"
  ) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             #col.ind = all_dist_pheno_fc$Age, 
             fill.ind = all_dist_pheno$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of All Cotwins' Distance between Landmarks") +
  scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic() +
  labs(title = "PCA of All Cotwins' Distance between Landmarks", fill = "Age") +
  theme(plot.title = element_text(size=30)) 
#theme(legend.position = "none")

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

fviz_eig(res.pca, choice = c("variance"), geom = c("bar"))

fviz_eig(res.pca, choice = c("eigenvalue"), geom = c("bar"))

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#Mclust w/ first 2 PCs
cluster.pca.adp.1.2 <- Mclust(res.pca[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adp.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.adp.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.adp.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adp.1.2, what = "uncertainty")
dev.off()

# using first 3 PCs

cluster.pca.adp.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.adp.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.adp.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.adp.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adp.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adp.1.2.3, what = "uncertainty")
dev.off()

all_dist_pheno_class12 <- all_dist_pheno
all_dist_pheno_class123 <- all_dist_pheno
all_dist_pheno_class12$mclust_clusters <- cluster.pca.adp.1.2$classification
all_dist_pheno_class123$mclust_clusters <- cluster.pca.adp.1.2.3$classification

all_dist_pheno_class12 %>% 
  group_by(mclust_clusters) %>%
  count(Sex) %>%
  ggplot(aes(x = mclust_clusters, y = n, fill = Sex)) +
  geom_bar(stat="identity", alpha=.6, width=.4) +
  theme_classic() +
  labs(x = "MCLUST Clusters", y = "", title = "Clustering of All Schoeller Individuals") +
  scale_x_discrete(limits=c("1","2")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank())

all_dist_pheno_class12 %>% 
  group_by(mclust_clusters) %>%
  count(Ethnicity) %>%
  group_by(mclust_clusters) %>%
  arrange(n) %>%
  ggplot(aes(x = mclust_clusters, y = n, fill = Ethnicity)) +
  geom_bar(stat="identity", alpha=.6, width=.4) +
  theme_classic() +
  labs(x = "MCLUST Clusters", y = "", title = "Clustering of All Schoeller Individuals") +
  scale_x_discrete(limits=c("1","2")) +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank())


#dist fc

scaled_all_dist_pheno_fc<-scale(all_dist_pheno_fc[,10:360])
scaled_all_dist_pheno_fc<-cbind(all_dist_pheno_fc[c(1:9)],scaled_all_dist_pheno_fc)
res.pca <- prcomp(scaled_all_dist_pheno_fc[,10:360])

#controls
#pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-PCA/Schoeller-dist-fc/schoeller_pca_dist_fc_controls_arranged.pdf", w=4, h=4)
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = TRUE, geom = "text",
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", 
             fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             title = "PCA of FC of Feature Distance") +
  theme(legend.position = "none",
        plot.title = element_text(size=30)) + geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID))
#dev.off()

#sex
groups <- as.factor(scaled_all_dist_pheno_fc$Sex) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups,
             #addEllipses = TRUE,
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             legend.title = "Sex", #theme(legend.position = "none") + 
             title = "PCA of FC of Feature Distance") +
  theme(plot.title = element_text(size=30))
#theme(legend.position = "none") +   
#geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#age
all_dist_pheno_fc_age <- all_dist_pheno_fc %>%
  mutate(FortyPlus = ifelse(Age > 40, 1, 0)) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = as.character(all_dist_pheno_fc_age$FortyPlus), 
             fill.ind = as.character(all_dist_pheno_fc_age$FortyPlus),
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of FC of Feature Distance") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  #scale_fill_discrete(name= "",labels = c("cowtins",'negative control','positive control'))+
  theme(plot.title = element_text(size=30),
        legend.position = "none"
  ) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = all_dist_pheno_fc$Age, 
             fill.ind = all_dist_pheno_fc$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of FC of Feature Distance") +
  scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic() +
  labs(title = "PCA of FC of Feature Distance", fill = "Age") +
  theme(plot.title = element_text(size=30)) 
#theme(legend.position = "none")

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

fviz_eig(res.pca, choice = c("variance"), geom = c("bar"))

fviz_eig(res.pca, choice = c("eigenvalue"), geom = c("bar"))

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(pca.adpfc[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

# using first 3 PCs

cluster.pca.adpfc.1.2.3 <- Mclust(pca.adpfc[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpfc.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2.3, what = "uncertainty")
dev.off()

all_dist_pheno_fc_test <- all_dist_pheno_fc
all_dist_pheno_fc_test$mclust_clusters <- cluster.pca.adpfc.1.2.3$classification

all_dist_pheno_fc_test_1 <- all_dist_pheno_fc_test %>% filter(mclust_clusters == 1)
all_dist_pheno_fc_test_2 <- all_dist_pheno_fc_test %>% filter(mclust_clusters == 2)


all_dist_pheno_fc_test  %>% count(Ethnicity)
#Ethnicity     n
# asian         3
# black        13
# caucasian    65
# hispanic      2
all_dist_pheno_fc_test  %>% count(Sex_Diff)
#Sex_Diff      n
# 0            35
# 1            46
# 2            2
all_dist_pheno_fc_test_1  %>% count(Ethnicity)
#Ethnicity       n
# asian         3
# black         10
# caucasian     46
# hispanic      2
all_dist_pheno_fc_test_1  %>% count(Sex_Diff)
#Sex_Diff     n
# 0           29
# 1           30
# 2            2
all_dist_pheno_fc_test_2  %>% count(Ethnicity)
#Ethnicity         n
# black         3
# caucasian     19
all_dist_pheno_fc_test_2  %>% count(Sex_Diff)
#Ethnicity         n
# black         6
# caucasian     16

#dist delta

scaled_all_dist_pheno_delta<-scale(all_dist_pheno_delta[,10:360])
scaled_all_dist_pheno_delta<-cbind(all_dist_pheno_delta[c(1:9)],scaled_all_dist_pheno_delta)
res.pca <- prcomp(scaled_all_dist_pheno_delta[,10:360])

#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = TRUE, geom = "text",
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", 
             fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             title = "PCA of Delta of Feature Distance") +
  theme(legend.position = "none",
        plot.title = element_text(size=30)) + geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID))

#sex
groups <- as.factor(scaled_all_dist_pheno_fc$Sex_Diff) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups,
             #addEllipses = TRUE,
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             legend.title = "Sex", #theme(legend.position = "none") + 
             title = "PCA of Delta of Feature Distance") +
  theme(plot.title = element_text(size=30))
#theme(legend.position = "none") +   
#geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

all_dist_pheno_delta_age <- all_dist_pheno_delta %>%
  mutate(FortyPlus = ifelse(Age > 40, 1, 0)) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = as.character(all_dist_pheno_fc_age$FortyPlus), 
             fill.ind = as.character(all_dist_pheno_fc_age$FortyPlus),
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Delta of Feature Distance") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  #scale_fill_discrete(name= "",labels = c("cowtins",'negative control','positive control'))+
  theme(plot.title = element_text(size=30),
        legend.position = "none"
  )

#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = all_dist_pheno_delta_age$FortyPlus, 
             fill.ind = all_dist_pheno_delta_age$FortyPlus,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Delta of Feature Distance") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  theme(legend.position = "none")

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = all_dist_pheno_delta$Age, 
             fill.ind = all_dist_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Delta of Feature Distance") +
  scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  theme(plot.title = element_text(size=30),
        legend.position = "none"
  )

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#Mclust w/ first 2 PCs

cluster.pca.adpd.1.2 <- Mclust(pca.adpd[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "uncertainty")
dev.off()

#using first 3 PCs

cluster.pca.adpd.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpd.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2.3, what = "uncertainty")
dev.off()

#dist absdelta

scaled_all_dist_pheno_absdelta<-scale(all_dist_pheno_absdelta[10:360])
scaled_all_dist_pheno_absdelta<-cbind(all_dist_pheno_absdelta[c(1:9)],scaled_all_dist_pheno_absdelta)
res.pca <- prcomp(scaled_all_dist_pheno_absdelta[10:360])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "PCA of Absolute Delta of Feature Distance") + 
  theme(legend.position = "none", plot.title = element_text(size=30)) +
  geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID))


#sex
groups <- as.factor(scaled_all_dist_pheno_absdelta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "PCA of Absolute Delta of Feature Distance") + theme(legend.position = "none",
                                                                          plot.title = element_text(size=30))
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = scaled_all_dist_pheno_absdelta$Age, 
             fill.ind = scaled_all_dist_pheno_absdelta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Absolute Delta of Feature Distance") +
  scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(size=30))


# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "ind", axes = 1:2, top = 120)  +
  theme(axis.text.x = element_text(angle=90))
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 20) +
  theme(axis.text.x = element_text(angle=90))

#Mclust w/ first 2 PCs

cluster.pca.adpdabs.1.2 <- Mclust(res.pca[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpdabs.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.adpdabs.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.adpdabs.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpdabs.1.2, what = "uncertainty")
dev.off()

#using first 3 PCs

cluster.pca.adpdabs.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpdabs.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.adpdabs.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.adpdabs.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.adpdabs.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpdabs.1.2.3, what = "uncertainty")
dev.off()

#ratio

all_ratio_pheno$Co_twins_ID <- as.integer(as.character(all_ratio_pheno$Co_twins_ID))

all_ratio_pheno %<>% arrange(Co_twins_ID)

scaled_all_ratio_pheno<-scale(all_ratio_pheno[,10:(ncol(all_ratio_pheno))])
scaled_all_ratio_pheno<-cbind(all_ratio_pheno[c(1:9)],scaled_all_ratio_pheno)
res.pca <- prcomp(scaled_all_ratio_pheno[,10:(ncol(all_ratio_pheno))])

#controls
#pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-PCA/Schoeller-dist-fc/schoeller_pca_dist_fc_controls_arranged.pdf", w=4, h=4)
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = TRUE, geom = "text",
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             col.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             title = "PCA of All Cotwins' Ratio of Distances between Landmarks") +
  theme(legend.position = "none",
        plot.title = element_text(size=30)) + 
  geom_text(aes(label = all_dist_pheno$Co_twins_ID))
#dev.off()

#sex
groups <- as.factor(all_dist_pheno$Sex) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             col.ind = groups, fill.ind = groups,
             #addEllipses = TRUE,
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             legend.title = "Sex", #theme(legend.position = "none") + 
             title = "PCA of All Cotwins' Ratio of Distances between Landmarks") +
  theme(plot.title = element_text(size=30))
#theme(legend.position = "none") +   
#geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#age
all_dist_pheno_age <- all_dist_pheno %>%
  mutate(FortyPlus = ifelse(Age > 40, 1, 0)) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             #col.ind = as.character(all_dist_pheno_fc_age$FortyPlus), 
             fill.ind = as.character(all_dist_pheno_age$FortyPlus),
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of All Cotwins' Ratio of Distances between Landmarks") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  #scale_fill_discrete(name= "",labels = c("cowtins",'negative control','positive control'))+
  theme(plot.title = element_text(size=30),
        legend.position = "none"
  ) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,24,24,24,24,24,
                            25,25,25,25,25,25,25,25,25,25),
             #col.ind = all_dist_pheno_fc$Age, 
             fill.ind = all_dist_pheno$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of All Cotwins' Ratio of Distances between Landmarks") +
  scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic() +
  labs(title = "PCA of All Cotwins' Ratio of Distances between Landmarks", fill = "Age") +
  theme(plot.title = element_text(size=30)) 
#theme(legend.position = "none")

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

fviz_eig(res.pca, choice = c("variance"), geom = c("bar"))

fviz_eig(res.pca, choice = c("eigenvalue"), geom = c("bar"))

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#Mclust w/ first 2 PCs
cluster.pca.arp.1.2 <- Mclust(res.pca[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adp.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.arp.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.arp.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arp.1.2, what = "uncertainty")
dev.off()

# using first 3 PCs

cluster.pca.arp.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.adp.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.arp.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.arp.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arp.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arp.1.2.3, what = "uncertainty")
dev.off()

#ratio delta

scaled_all_ratio_pheno_delta<-scale(all_ratio_pheno_delta[,10:61434])
scaled_all_ratio_pheno_delta<-cbind(all_ratio_pheno_delta[c(1:9)],scaled_all_ratio_pheno_delta)
res.pca <- prcomp(scaled_all_ratio_pheno_delta[,10:61434])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = TRUE, geom = "text",
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", 
             fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             title = "PCA of Delta of Feature Distance Ratios") +
  theme(legend.position = "none",
        plot.title = element_text(size=30)) + geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID))

#sex
groups <- as.factor(scaled_all_ratio_pheno_delta$Sex_Diff) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups,
             #addEllipses = TRUE,
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             legend.title = "Sex", #theme(legend.position = "none") + 
             title = "PCA of Delta of Feature Distance Ratio") +
  theme(plot.title = element_text(size=30))
#theme(legend.position = "none") +   
#geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_dist_pheno_delta_age$FortyPlus,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Delta of Feature Distance Ratio") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  theme(plot.title = element_text(size=30), 
        legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=75))

#Mclust w/ first 2 PCs

cluster.pca.arpd.1.2 <- Mclust(res.pca[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.arpd.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.arpd.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.arpd.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arpd.1.2, what = "uncertainty")
dev.off()

# using first 3 pca eigenvectors

cluster.pca.arpd.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.arpd.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.arpd.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.arpd.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpd.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arpd.1.2.3, what = "uncertainty")
dev.off()

#ratio fc

scaled_all_ratio_pheno_fc<-scale(all_ratio_pheno_fc[,10:360])
scaled_all_ratio_pheno_fc<-cbind(all_ratio_pheno_fc[c(1:9)],scaled_all_ratio_pheno_fc)
res.pca <- prcomp(scaled_all_ratio_pheno_fc[,10:360])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = TRUE, geom = "text",
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", 
             fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             title = "PCA of FC of Feature Distance Ratio") +
  theme(legend.position = "none",
        plot.title = element_text(size=30)) + geom_text(aes(label = scaled_all_ratio_pheno_fc$Co_twins_ID))

#sex
groups <- as.factor(scaled_all_ratio_pheno_fc$Sex_Diff) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups,
             #addEllipses = TRUE,
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             legend.title = "Sex", #theme(legend.position = "none") + 
             title = "PCA of FC of Feature Distance Ratio") +
  theme(plot.title = element_text(size=30))
#theme(legend.position = "none") +   
#geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = all_dist_pheno_delta_age$FortyPlus, 
             fill.ind = all_dist_pheno_delta_age$FortyPlus,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Delta of Feature Distance Ratio") +
  #scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  theme(plot.title = element_text(size=30), legend.position = "none")

#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             #col.ind = all_dist_pheno_delta$Age, 
             fill.ind = scaled_all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of Delta of Feature Distance Ratio") +
  scale_fill_gradient(heat.colors(5)) +
  #theme_classic()+
  theme(plot.title = element_text(size=30),
        legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

# using first 2 PC eigenvectors

cluster.pca.arpfc.1.2 <- Mclust(res.pca[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.arpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.arpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.arpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arpfc.1.2, what = "uncertainty")
dev.off()

# using first 3 PC eigenvectors

cluster.pca.arpfc.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.arpfc.1.2.3, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.arpfc.1.2.3, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.arpfc.1.2.3, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpfc.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arpfc.1.2.3, what = "uncertainty")
dev.off()

#ratio absdelta
scaled_all_ratio_pheno_absdelta<-scale(all_ratio_pheno_absdelta[10:61434])
scaled_all_ratio_pheno_absdelta-cbind(all_dist_pheno_absdelta[c(1:9)],scaled_all_ratio_pheno_absdelta)
res.pca <- prcomp(scaled_all_ratio_pheno_absdelta[10:61434])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "PCA of Delta of Feature Distance Ratio") + 
  theme(legend.position = "none") +
  geom_text(aes(label=all_ratio_pheno_delta[,1])) 
#sex
groups <- as.factor(all_ratio_pheno_delta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x delta between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x delta between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(pca, choice = "var", axes = 1:2, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#Mclust with first 2 PCs

cluster.pca.arpdabs.1.2 <- Mclust(res.pca[["x"]][,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpdabs.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.arpdabs.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpdabs.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.arpdabs.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpdabs.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.arpdabs.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.pca.arpdabs.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arpdabs.1.2, what = "uncertainty")
dev.off()

# using first 3 PCs 

cluster.pca.arpdabs.1.2.3 <- Mclust(res.pca[["x"]][,1:3], prior = priorControl())
pdf("cluster.pca.arpdabs.1.2.3.BIC.pdf", w=6, h=4)
plot(cluster.pca.arpdabs.1.2.3, what = "BIC")
dev.off()
pdf("cluster.pca.arpdabs.1.2.3.density.pdf", w=4, h=4)
plot(cluster.pca.arpdabs.1.2.3, what = "density")
dev.off()
pdf("cluster.pca.arpdabs.1.2.3.classification.pdf", w=4, h=4)
plot(cluster.pca.arpdabs.1.2.3, what = "classification")
dev.off()
pdf("cluster.pca.arpdabs.1.2.3.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.arpdabs.1.2.3, what = "uncertainty")
dev.off()

###### Sparse PCA ########
##Originally written on the MNP Browser##
library(tidyverse)
options(stringsAsFactors = FALSE)
setwd("/secondary/projects/mnp/warren/sparsity-pca/Plots")
###### Sparse PCA elastic net distance##########

library(elasticnet)

# fc

scaled_adpfc <- scale(all_dist_pheno_fc[,10:360])

sparse.pca.result.all <- 
  elasticnet::spca(scaled_adpfc, K = 2, type = "predictor", sparse = "varnum", para = c(351, 351))
sparse.pca.result.10 <- 
  elasticnet::spca(scaled_adpfc, K = 2, type = "predictor", sparse = "varnum", para = c(10, 10))
sparse.pca.result.30 <- 
  elasticnet::spca(scaled_adpfc, K = 2, type = "predictor", sparse = "varnum", para = c(30, 30))
sparse.pca.result.65 <- 
  elasticnet::spca(scaled_adpfc, K = 2, type = "predictor", sparse = "varnum", para = c(65, 65))
sparse.pca.result.100 <- 
  elasticnet::spca(scaled_adpfc, K = 2, type = "predictor", sparse = "varnum", para = c(100,100))

spca.loadings.all <- unlist(sparse.pca.result.all[["loadings"]])
spca.loadings.10 <- unlist(sparse.pca.result.10[["loadings"]])
spca.loadings.30 <- unlist(sparse.pca.result.30[["loadings"]])
spca.loadings.65 <- unlist(sparse.pca.result.65[["loadings"]])
spca.loadings.100 <- unlist(sparse.pca.result.100[["loadings"]])

mat_spca.loadings.all <- as.matrix(spca.loadings.all)
mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpfc %*% mat_spca.loadings.all
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.all[["pev"]][1], 3) *100
pc2_pev <- round(sparse.pca.result.all[["pev"]][2], 3) *100
pdf("schoeller_spca_dist_fc_controls_arranged_allvar.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.spca.adpfc.1.2.BIC_all.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("cluster.spca.adpfc.1.2.density_all.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.classification_all.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.uncertainty_all.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.10 <- as.matrix(spca.loadings.10)
#mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpfc %*% spca.loadings.10
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.10[["pev"]][1], 3)*100
pc2_pev <- round(sparse.pca.result.10[["pev"]][2], 3)*100
pdf("schoeller_spca_dist_fc_controls_arranged_10var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.spca.adpfc.1.2.BIC_loadings10.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("cluster.spca.adpfc.1.2.density_loadings10.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.classification_loadings10.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.uncertainty_loadings10.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.30 <- as.matrix(spca.loadings.30)
mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpfc %*% mat_spca.loadings.30
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.30[["pev"]][1], 3)*100
pc2_pev <- round(sparse.pca.result.30[["pev"]][2], 3)*100
pdf("schoeller_spca_dist_fc_controls_arranged_30var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.spca.adpfc.1.2.BIC_loadings30.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("cluster.spca.adpfc.1.2.density_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.classification_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.uncertainty_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.65 <- as.matrix(spca.loadings.65)
mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpfc %*% mat_spca.loadings.65
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.65[["pev"]][1], 3)*100
pc2_pev <- round(sparse.pca.result.65[["pev"]][2], 3)*100
pdf("schoeller_spca_dist_fc_controls_arranged_65var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.spca.adpfc.1.2.BIC_loadings65.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("cluster.spca.adpfc.1.2.density_loadings65.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.classification_loadings65.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.uncertainty_loadings65.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.100 <- as.matrix(spca.loadings.100)
#mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpfc %*% mat_spca.loadings.100
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.100[["pev"]][1], 3)*100
pc2_pev <- round(sparse.pca.result.100[["pev"]][2], 3)*100
pdf("schoeller_spca_dist_fc_controls_arranged_100var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.pca.adpfc.1.2.BIC_loadings100.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("cluster.pca.adpfc.1.2.density_loadings100.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("cluster.pca.adpfc.1.2.plot.classification_loadings100.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("cluster.pca.adpfc.1.2.plot.uncertainty_loadings100.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

spca.loadings.30 <- as.data.frame(spca.loadings.30)
spca.loadings.65 <- as.data.frame(spca.loadings.65)
spca.loadings.10 <- as.data.frame(spca.loadings.10)
spca.loadings.100 <- as.data.frame(spca.loadings.100)
spca.loadings.all <- as.data.frame(spca.loadings.all)

sparse.30.PCs <- spca.loadings.30[(which(spca.loadings.30$PC1 != 0 | spca.loadings.30$PC2 != 0 )),]
sparse.30.PC1 <- spca.loadings.30[which(spca.loadings.30$PC1 != 0),]
sparse.30.PC2 <- spca.loadings.30[which(spca.loadings.30$PC2 != 0),]
sparse.65.PCs <- spca.loadings.65[which(spca.loadings.65$PC1 != 0 | spca.loadings.65$PC2 != 0 ),]
sparse.65.PC1 <- spca.loadings.65[which(spca.loadings.65$PC1 != 0),]
sparse.65.PC2 <- spca.loadings.65[which(spca.loadings.65$PC2 != 0),]
sparse.all.PC1 <- spca.loadings.all[which(spca.loadings.all$PC1 != 0),]
sparse.all.PC2 <- spca.loadings.all[which(spca.loadings.all$PC2 != 0),]
sparse.10.PC1 <- spca.loadings.10[which(spca.loadings.10$PC1 != 0),]
sparse.10.PC2 <- spca.loadings.10[which(spca.loadings.10$PC2 != 0),]
sparse.100.PC1 <- spca.loadings.100[which(spca.loadings.100$PC1 != 0),]
sparse.100.PC2 <- spca.loadings.100[which(spca.loadings.100$PC2 != 0),]

sparse.30.PC1 <- as.data.frame(sparse.30.PC1)
sparse.30.PC2 <- as.data.frame(sparse.30.PC2)
sparse.65.PC1 <- as.data.frame(sparse.65.PC1)
sparse.65.PC2 <- as.data.frame(sparse.65.PC2)
sparse.all.PC1 <- as.data.frame(sparse.all.PC1)
sparse.all.PC2 <- as.data.frame(sparse.all.PC2)
sparse.10.PC1 <- as.data.frame(sparse.10.PC1)
sparse.10.PC2 <- as.data.frame(sparse.10.PC2)
sparse.100.PC1 <- as.data.frame(sparse.100.PC1)
sparse.100.PC2 <- as.data.frame(sparse.100.PC2)

# load the library
library(forcats)

# Reorder following the value of another column:
pdf("schoeller_spca_dist_fc_30var_PC1.pdf", w=10, h=8)
sparse.30.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_30var_PC2.pdf", w=10, h=8)
sparse.30.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_65var_PC1.pdf", w=10, h=8)
sparse.65.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_65var_PC2.pdf", w=10, h=8)
sparse.65.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_10var_PC1.pdf", w=10, h=8)
sparse.10.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_10var_PC2.pdf", w=10, h=8)
sparse.10.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_100var_PC1.pdf", w=10, h=8)
sparse.100.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_100var_PC2.pdf", w=10, h=8)
sparse.100.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

#delta

scaled_adpd <- scale(all_dist_pheno_delta[,10:360])

sparse.pca.result.all <- 
  elasticnet::spca(scaled_adpd, K = 2, type = "predictor", sparse = "varnum", para = c(351, 351))
sparse.pca.result.10 <- 
  elasticnet::spca(scaled_adpd, K = 2, type = "predictor", sparse = "varnum", para = c(10, 10))
sparse.pca.result.30 <- 
  elasticnet::spca(scaled_adpd, K = 2, type = "predictor", sparse = "varnum", para = c(30, 30))
sparse.pca.result.65 <- 
  elasticnet::spca(scaled_adpd, K = 2, type = "predictor", sparse = "varnum", para = c(65, 65))
sparse.pca.result.100 <- 
  elasticnet::spca(scaled_adpd, K = 2, type = "predictor", sparse = "varnum", para = c(100,100))

spca.loadings.all <- unlist(sparse.pca.result.all[["loadings"]])
spca.loadings.10 <- unlist(sparse.pca.result.10[["loadings"]])
spca.loadings.30 <- unlist(sparse.pca.result.30[["loadings"]])
spca.loadings.65 <- unlist(sparse.pca.result.65[["loadings"]])
spca.loadings.100 <- unlist(sparse.pca.result.100[["loadings"]])

mat_spca.loadings.all <- as.matrix(spca.loadings.all)
mat_scaled_adpd <- as.matrix(scaled_adpd)

test <- mat_scaled_adpd %*% mat_spca.loadings.all
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.all[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.all[["pev"]][2], 3) * 100
pdf("schoeller_spca_dist_delta_controls_arranged_allvar.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpd.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.pca.adpd.1.2.BIC_loadingsAll.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2, what = "BIC")
dev.off()
pdf("cluster.pca.adpd.1.2.density_loadingsAll.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "density")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.classification_loadingsAll.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "classification")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.uncertainty_loadingsAll.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.10 <- as.matrix(spca.loadings.10)
#mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpd %*% spca.loadings.10
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.10[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.10[["pev"]][2], 3) * 100
pdf("schoeller_spca_dist_delta_controls_arranged_10var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpd.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.pca.adpd.1.2.BIC_loadings10.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2, what = "BIC")
dev.off()
pdf("cluster.pca.adpd.1.2.density_loadings10.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "density")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.classification_loadings10.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "classification")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.uncertainty_loadings10.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.30 <- as.matrix(spca.loadings.30)

test <- mat_scaled_adpd %*% mat_spca.loadings.30
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.30[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.30[["pev"]][2], 3) * 100
pdf("schoeller_spca_dist_delta_controls_arranged_30var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpd.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.pca.adpd.1.2.BIC_loadings30.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2, what = "BIC")
dev.off()
pdf("cluster.pca.adpd.1.2.density_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "density")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.classification_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "classification")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.uncertainty_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.65 <- as.matrix(spca.loadings.65)

test <- mat_scaled_adpd %*% mat_spca.loadings.65
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.65[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.65[["pev"]][2], 3) * 100
pdf("schoeller_spca_dist_delta_controls_arranged_65var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpd.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.pca.adpd.1.2.BIC_loadings65.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2, what = "BIC")
dev.off()
pdf("cluster.pca.adpd.1.2.density_loadings65.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "density")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.classification_loadings65.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "classification")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.uncertainty_loadings65.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.100 <- as.matrix(spca.loadings.100)

test <- mat_scaled_adpd %*% mat_spca.loadings.100
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.100[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.100[["pev"]][2], 3) * 100
pdf("schoeller_spca_dist_delta_controls_arranged_100var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpd.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.pca.adpd.1.2.BIC_loadings100.pdf", w=6, h=4)
plot(cluster.pca.adpd.1.2, what = "BIC")
dev.off()
pdf("cluster.pca.adpd.1.2.density_loadings100.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "density")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.classification_loadings100.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "classification")
dev.off()
pdf("cluster.pca.adpd.1.2.plot.uncertainty_loadings100.pdf", w=4, h=4)
plot(cluster.pca.adpd.1.2, what = "uncertainty")
dev.off()

spca.loadings.30 <- as.data.frame(spca.loadings.30)
spca.loadings.65 <- as.data.frame(spca.loadings.65)
spca.loadings.10 <- as.data.frame(spca.loadings.10)
spca.loadings.100 <- as.data.frame(spca.loadings.100)
spca.loadings.all <- as.data.frame(spca.loadings.all)

sparse.30.PCs <- spca.loadings.30[(which(spca.loadings.30$PC1 != 0 | spca.loadings.30$PC2 != 0 )),]
sparse.30.PC1 <- spca.loadings.30[which(spca.loadings.30$PC1 != 0),]
sparse.30.PC2 <- spca.loadings.30[which(spca.loadings.30$PC2 != 0),]
sparse.65.PCs <- spca.loadings.65[which(spca.loadings.65$PC1 != 0 | spca.loadings.65$PC2 != 0 ),]
sparse.65.PC1 <- spca.loadings.65[which(spca.loadings.65$PC1 != 0),]
sparse.65.PC2 <- spca.loadings.65[which(spca.loadings.65$PC2 != 0),]
sparse.all.PC1 <- spca.loadings.all[which(spca.loadings.all$PC1 != 0),]
sparse.all.PC2 <- spca.loadings.all[which(spca.loadings.all$PC2 != 0),]
sparse.10.PC1 <- spca.loadings.10[which(spca.loadings.10$PC1 != 0),]
sparse.10.PC2 <- spca.loadings.10[which(spca.loadings.10$PC2 != 0),]
sparse.100.PC1 <- spca.loadings.100[which(spca.loadings.100$PC1 != 0),]
sparse.100.PC2 <- spca.loadings.100[which(spca.loadings.100$PC2 != 0),]

sparse.30.PC1 <- as.data.frame(sparse.30.PC1)
sparse.30.PC2 <- as.data.frame(sparse.30.PC2)
sparse.65.PC1 <- as.data.frame(sparse.65.PC1)
sparse.65.PC2 <- as.data.frame(sparse.65.PC2)
sparse.all.PC1 <- as.data.frame(sparse.all.PC1)
sparse.all.PC2 <- as.data.frame(sparse.all.PC2)
sparse.10.PC1 <- as.data.frame(sparse.10.PC1)
sparse.10.PC2 <- as.data.frame(sparse.10.PC2)
sparse.100.PC1 <- as.data.frame(sparse.100.PC1)
sparse.100.PC2 <- as.data.frame(sparse.100.PC2)

# load the library
library(forcats)

# Reorder following the value of another column:
pdf("schoeller_spca_dist_delta_30var_PC1.pdf", w=10, h=8)
sparse.30.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_30var_PC2.pdf", w=10, h=8)
sparse.30.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_65var_PC1.pdf", w=10, h=8)
sparse.65.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_65var_PC2.pdf", w=10, h=8)
sparse.65.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_10var_PC1.pdf", w=10, h=8)
sparse.10.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_10var_PC2.pdf", w=10, h=8)
sparse.10.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_100var_PC1.pdf", w=10, h=16)
sparse.100.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_delta_100var_PC2.pdf", w=10, h=16)
sparse.100.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

###### Sparse PCA elastic net ratio#############

#delta

scaled_arpd <- scale(all_ratio_pheno_delta[,10:61434])

sparse.pca.result.all <- 
  elasticnet::spca(scaled_arpd, K = 2, type = "predictor", sparse = "varnum", para = c(61434, 61434))
sparse.pca.result.10 <- 
  elasticnet::spca(scaled_arpd, K = 2, type = "predictor", sparse = "varnum", para = c(10, 10))
sparse.pca.result.30 <- 
  elasticnet::spca(scaled_arpd, K = 2, type = "predictor", sparse = "varnum", para = c(30, 30))
sparse.pca.result.65 <- 
  elasticnet::spca(scaled_arpd, K = 2, type = "predictor", sparse = "varnum", para = c(65, 65))
sparse.pca.result.100 <- 
  elasticnet::spca(scaled_arpd, K = 2, type = "predictor", sparse = "varnum", para = c(100,100))

spca.loadings.all <- unlist(sparse.pca.result.all[["loadings"]])
spca.loadings.10 <- unlist(sparse.pca.result.10[["loadings"]])
spca.loadings.30 <- unlist(sparse.pca.result.30[["loadings"]])
spca.loadings.65 <- unlist(sparse.pca.result.65[["loadings"]])
spca.loadings.100 <- unlist(sparse.pca.result.100[["loadings"]])

mat_spca.loadings.all <- as.matrix(spca.loadings.all)
mat_scaled_arpd <- as.matrix(scaled_arpd)

test <- mat_scaled_arpd %*% mat_spca.loadings.all
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.all[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.all[["pev"]][2], 3) * 100
pdf("schoeller_spca_ratio_delta_controls_arranged_allvar.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

mat_spca.loadings.10 <- as.matrix(spca.loadings.10)
#mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_adpd %*% spca.loadings.10
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.10[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.10[["pev"]][2], 3) * 100
pdf("schoeller_spca_ratio_delta_controls_arranged_10var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

mat_spca.loadings.30 <- as.matrix(spca.loadings.30)

test <- mat_scaled_arpd %*% mat_spca.loadings.30
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.30[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.30[["pev"]][2], 3) * 100
pdf("schoeller_spca_ratio_delta_controls_arranged_30var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.65 <- as.matrix(spca.loadings.65)

test <- mat_scaled_arpd %*% mat_spca.loadings.65
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.65[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.65[["pev"]][2], 3) * 100
pdf("schoeller_spca_ratio_delta_controls_arranged_65var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.BIC.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/MSU-mclust/cluster.pca.adpfc.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

mat_spca.loadings.100 <- as.matrix(spca.loadings.100)

test <- mat_scaled_adpd %*% mat_spca.loadings.100
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.100[["pev"]][1], 3) * 100
pc2_pev <- round(sparse.pca.result.100[["pev"]][2], 3) * 100
pdf("schoeller_spca_ratio_delta_controls_arranged_100var.pdf", w=12, h=10)
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

#Mclust w/ first 2 PCs
cluster.pca.adpfc.1.2 <- Mclust(test[,1:2], prior = priorControl())
pdf("cluster.spca.adpfc.1.2.BIC_loadings30.pdf", w=6, h=4)
plot(cluster.pca.adpfc.1.2, what = "BIC")
dev.off()
pdf("cluster.spca.adpfc.1.2.density_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "density")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.classification_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "classification")
dev.off()
pdf("cluster.spca.adpfc.1.2.plot.uncertainty_loadings30.pdf", w=4, h=4)
plot(cluster.pca.adpfc.1.2, what = "uncertainty")
dev.off()

spca.loadings.30 <- as.data.frame(spca.loadings.30)
spca.loadings.65 <- as.data.frame(spca.loadings.65)
spca.loadings.10 <- as.data.frame(spca.loadings.10)
spca.loadings.100 <- as.data.frame(spca.loadings.100)
spca.loadings.all <- as.data.frame(spca.loadings.all)

sparse.30.PCs <- spca.loadings.30[(which(spca.loadings.30$PC1 != 0 | spca.loadings.30$PC2 != 0 )),]
sparse.30.PC1 <- spca.loadings.30[which(spca.loadings.30$PC1 != 0),]
sparse.30.PC2 <- spca.loadings.30[which(spca.loadings.30$PC2 != 0),]
sparse.65.PCs <- spca.loadings.65[which(spca.loadings.65$PC1 != 0 | spca.loadings.65$PC2 != 0 ),]
sparse.65.PC1 <- spca.loadings.65[which(spca.loadings.65$PC1 != 0),]
sparse.65.PC2 <- spca.loadings.65[which(spca.loadings.65$PC2 != 0),]
sparse.all.PC1 <- spca.loadings.all[which(spca.loadings.all$PC1 != 0),]
sparse.all.PC2 <- spca.loadings.all[which(spca.loadings.all$PC2 != 0),]
sparse.10.PC1 <- spca.loadings.10[which(spca.loadings.10$PC1 != 0),]
sparse.10.PC2 <- spca.loadings.10[which(spca.loadings.10$PC2 != 0),]
sparse.100.PC1 <- spca.loadings.100[which(spca.loadings.100$PC1 != 0),]
sparse.100.PC2 <- spca.loadings.100[which(spca.loadings.100$PC2 != 0),]

sparse.30.PC1 <- as.data.frame(sparse.30.PC1)
sparse.30.PC2 <- as.data.frame(sparse.30.PC2)
sparse.65.PC1 <- as.data.frame(sparse.65.PC1)
sparse.65.PC2 <- as.data.frame(sparse.65.PC2)
sparse.all.PC1 <- as.data.frame(sparse.all.PC1)
sparse.all.PC2 <- as.data.frame(sparse.all.PC2)
sparse.10.PC1 <- as.data.frame(sparse.10.PC1)
sparse.10.PC2 <- as.data.frame(sparse.10.PC2)
sparse.100.PC1 <- as.data.frame(sparse.100.PC1)
sparse.100.PC2 <- as.data.frame(sparse.100.PC2)

# load the library
library(forcats)

# Reorder following the value of another column:
pdf("schoeller_spca_ratio_delta_30var_PC1.pdf", w=10, h=8)
sparse.30.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_ratio_delta_30var_PC2.pdf", w=10, h=8)
sparse.30.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_ratio_delta_65var_PC1.pdf", w=10, h=8)
sparse.65.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_ratio_delta_65var_PC2.pdf", w=10, h=8)
sparse.65.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_ratio_delta_10var_PC1.pdf", w=10, h=8)
sparse.10.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_ratio_delta_10var_PC2.pdf", w=10, h=8)
sparse.10.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

# fc

scaled_arpfc <- scale(all_ratio_pheno_fc[,10:360])

sparse.pca.result.10 <- 
  elasticnet::spca(scaled_arpfc, K = 2, type = "predictor", sparse = "varnum", para = c(10, 10))
sparse.pca.result.30 <- 
  elasticnet::spca(scaled_arpfc, K = 2, type = "predictor", sparse = "varnum", para = c(30, 30))
sparse.pca.result.65 <- 
  elasticnet::spca(scaled_arpfc, K = 2, type = "predictor", sparse = "varnum", para = c(65, 65))
sparse.pca.result.100 <- 
  elasticnet::spca(scaled_arpfc, K = 2, type = "predictor", sparse = "varnum", para = c(100,100))

spca.loadings.10 <- unlist(sparse.pca.result.10[["loadings"]])
spca.loadings.30 <- unlist(sparse.pca.result.30[["loadings"]])
spca.loadings.65 <- unlist(sparse.pca.result.65[["loadings"]])
spca.loadings.100 <- unlist(sparse.pca.result.100[["loadings"]])

mat_spca.loadings.all <- as.matrix(spca.loadings.all)
mat_scaled_arpfc <- as.matrix(scaled_arpfc)

test <- mat_scaled_arpfc %*% mat_spca.loadings.all
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.all[["pev"]][1], 3)
pc2_pev <- round(sparse.pca.result.all[["pev"]][2], 3)
pdf("schoeller_spca_ratio_fc_controls_arranged_allvar.pdf")
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

mat_spca.loadings.10 <- as.matrix(spca.loadings.10)
#mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_arpfc %*% spca.loadings.10
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.10[["pev"]][1], 3)
pc2_pev <- round(sparse.pca.result.10[["pev"]][2], 3)
pdf("schoeller_spca_ratio_fc_controls_arranged_10var.pdf")
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

mat_spca.loadings.30 <- as.matrix(spca.loadings.30)
#mat_scaled_arpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_arpfc %*% mat_spca.loadings.30
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.30[["pev"]][1], 3)
pc2_pev <- round(sparse.pca.result.30[["pev"]][2], 3)
pdf("schoeller_spca_ratio_fc_controls_arranged_30var.pdf")
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

mat_spca.loadings.65 <- as.matrix(spca.loadings.65)
#mat_scaled_arpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_arpfc %*% mat_spca.loadings.65
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.65[["pev"]][1], 3)
pc2_pev <- round(sparse.pca.result.65[["pev"]][2], 3)
pdf("schoeller_spca_ratio_fc_controls_arranged_65var.pdf")
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

mat_spca.loadings.100 <- as.matrix(spca.loadings.100)
#mat_scaled_adpfc <- as.matrix(scaled_adpfc)

test <- mat_scaled_arpfc %*% mat_spca.loadings.100
test <- as.data.frame(test)

pc1_pev <- round(sparse.pca.result.100[["pev"]][1], 3)
pc2_pev <- round(sparse.pca.result.100[["pev"]][2], 3)
pdf("schoeller_spca_ratio_fc_controls_arranged_100var.pdf")
ggplot(test) + geom_point(aes(x=test$PC1, y=test$PC2, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                              "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                              fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                     "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                          alpha = 0.3, size = 10, shape=21) + 
  theme_minimal() + # alpha is the opacity
  geom_text(aes(x=PC1, y=PC2, label=all_dist_pheno_fc$Co_twins_ID)) +
  labs( x=bquote("sPCA dim1 (" ~ .(pc1_pev) ~ "%)"), y=bquote("sPCA dim2 (" ~ .(pc2_pev) ~ "%)"), color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))
dev.off()

spca.loadings.30 <- as.data.frame(spca.loadings.30)
spca.loadings.65 <- as.data.frame(spca.loadings.65)
spca.loadings.10 <- as.data.frame(spca.loadings.10)
spca.loadings.100 <- as.data.frame(spca.loadings.100)
spca.loadings.all <- as.data.frame(spca.loadings.all)

sparse.30.PCs <- spca.loadings.30[(which(spca.loadings.30$PC1 != 0 | spca.loadings.30$PC2 != 0 )),]
sparse.30.PC1 <- spca.loadings.30[which(spca.loadings.30$PC1 != 0),]
sparse.30.PC2 <- spca.loadings.30[which(spca.loadings.30$PC2 != 0),]
sparse.65.PCs <- spca.loadings.65[which(spca.loadings.65$PC1 != 0 | spca.loadings.65$PC2 != 0 ),]
sparse.65.PC1 <- spca.loadings.65[which(spca.loadings.65$PC1 != 0),]
sparse.65.PC2 <- spca.loadings.65[which(spca.loadings.65$PC2 != 0),]
sparse.all.PC1 <- spca.loadings.30[which(spca.loadings.all$PC1 != 0),]
sparse.all.PC2 <- spca.loadings.30[which(spca.loadings.all$PC2 != 0),]
sparse.10.PC1 <- spca.loadings.30[which(spca.loadings.10$PC1 != 0),]
sparse.10.PC2 <- spca.loadings.30[which(spca.loadings.10$PC2 != 0),]
sparse.100.PC1 <- spca.loadings.30[which(spca.loadings.100$PC1 != 0),]
sparse.100.PC2 <- spca.loadings.30[which(spca.loadings.100$PC2 != 0),]

sparse.30.PC1 <- as.data.frame(sparse.30.PC1)
sparse.30.PC2 <- as.data.frame(sparse.30.PC2)
sparse.65.PC1 <- as.data.frame(sparse.65.PC1)
sparse.65.PC2 <- as.data.frame(sparse.65.PC2)
sparse.all.PC1 <- as.data.frame(sparse.all.PC1)
sparse.all.PC2 <- as.data.frame(sparse.all.PC2)
sparse.10.PC1 <- as.data.frame(sparse.10.PC1)
sparse.10.PC2 <- as.data.frame(sparse.10.PC2)
sparse.100.PC1 <- as.data.frame(sparse.100.PC1)
sparse.100.PC2 <- as.data.frame(sparse.100.PC2)

# load the library
library(forcats)

# Reorder following the value of another column:
pdf("schoeller_spca_dist_fc_30var_PC1.pdf", w=10, h=8)
sparse.30.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_30var_PC2.pdf", w=10, h=8)
sparse.30.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_65var_PC1.pdf", w=10, h=8)
sparse.65.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_65var_PC2.pdf", w=10, h=8)
sparse.65.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_10var_PC1.pdf", w=10, h=8)
sparse.10.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_10var_PC2.pdf", w=10, h=8)
sparse.10.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_100var_PC1.pdf", w=10, h=8)
sparse.100.PC1 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC1)) %>%
  ggplot( aes(x=variables, y=PC1)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()

pdf("schoeller_spca_dist_fc_100var_PC1.pdf", w=10, h=8)
sparse.100.PC2 %>%
  rownames_to_column(var = "variables") %>%
  mutate(variables = fct_reorder(variables,PC2)) %>%
  ggplot( aes(x=variables, y=PC2)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()
dev.off()





###### t-SNE dist#######
library(tidyverse)
library(Rtsne)
library(ggplot2)

setwd("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-tSNE/")

#all dist pheno fc

scaled_all_dist_pheno_fc<-scale(all_dist_pheno_fc[10:360])
scaled_all_dist_pheno_fc<-cbind(all_dist_pheno_fc[1:9],scaled_all_dist_pheno_fc)

tsne <- Rtsne(scaled_all_dist_pheno_fc[,10:360],check_duplicates = FALSE, perplexity=12, max_iter = 20000, theta = 0.0) 

#plot t-sne results (you put it in a plot)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
pdf("schoeller_rtsne_dist_fc_controls_arranged_pplx12.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_fc$Co_twins_ID)) +
  labs( x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_dist_fc_sex_arranged_pplx12.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_fc$Sex, 
                                   fill=scaled_all_dist_pheno_fc$Sex), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_fc$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_dist_fc_age-gradient_arranged_pplx12.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_fc$Age, 
                                   fill=scaled_all_dist_pheno_fc$Age), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_fc$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#, 
#legend.position = "none")
dev.off()

#Mclust w/ first tSNE x and y dimensions
cluster.tsne.adpfc.1.2 <- Mclust(tsne_plot[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpfc.1.2.BIC.pdf", w=6, h=4)
plot(cluster.tsne.adpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.tsne.adpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.tsne.adpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpfc.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.tsne.adpfc.1.2, what = "uncertainty")
dev.off()

#dist delta


scaled_all_dist_pheno_delta<-scale(all_dist_pheno_delta[10:360])
scaled_all_dist_pheno_delta<-cbind(all_dist_pheno_delta[1:9],scaled_all_dist_pheno_delta)

tsne <- Rtsne(scaled_all_dist_pheno_delta[,10:360],check_duplicates = FALSE, perplexity=20, max_iter = 20000, theta = 0.0) 

#plot t-sne results (you put it in a plot)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
pdf("schoeller_rtsne_dist_delta_controls_arranged_pplx20.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_dist_delta_age-gradient_arranged_pplx20.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_delta$Sex, 
                                   fill=scaled_all_dist_pheno_delta$Sex), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  #geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_dist_delta_sex_arranged_pplx20.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_delta$Age, 
                                   fill=scaled_all_dist_pheno_delta$Age), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

#Mclust w/ tSNE dimensions

cluster.tsne.adpd.1.2 <- Mclust(tsne_plot[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpd.1.2.BIC.pdf", w=6, h=4)
plot(cluster.tsne.adpd.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpd.1.2.density.pdf", w=4, h=4)
plot(cluster.tsne.adpd.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpd.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.tsne.adpd.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.adpd.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.tsne.adpd.1.2, what = "uncertainty")
dev.off()


###### t-SNE ratio#############


#import of tsnes from server and plots of ratio matrices 

arpd_tsne_tab <- read.table(sfile = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_all_ratio_delta.tab", sep = "\t", header = F)
arpfc_tsne_tab <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_all_ratio_fc.tab", sep = "\t", header = F)
arpd_arranged_tsne_tab <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_delta.tab", sep = "\t", header = F)
arpfc_arranged_tsne_tab <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_fc.tab", sep = "\t", header = F)
arpdabs_arranged_tsne_tab <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_absdelta.tab", sep = "\t", header = F)
arpd_arranged_tsne_tab_pplx10 <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_delta_pplx10.tab", sep = "\t", header = F)
arpfc_arranged_tsne_tab_pplx10 <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_fc_pplx10.tab", sep = "\t", header = F)
arpd_arranged_tsne_tab_pplx20 <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_delta_pplx20.tab", sep = "\t", header = F)
arpfc_arranged_tsne_tab_pplx20 <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_fc_pplx20.tab", sep = "\t", header = F)
arpd_arranged_tsne_tab_pplx40 <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_delta_pplx40.tab", sep = "\t", header = F)
arpfc_arranged_tsne_tab_pplx40 <- read.table(file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/tsne_schoeller_all_ratio_pheno_fc_pplx40.tab", sep = "\t", header = F)


arpd_tsne <- cbind(all_dist_pheno_delta[,1:9], arpd_arranged_tsne_tab)
arpfc_tsne_pplx10 <- cbind(all_dist_pheno_fc[,1:9], arpfc_arranged_tsne_tab_pplx10)
arpd_tsne_pplx10 <- cbind(all_dist_pheno_delta[,1:9], arpd_arranged_tsne_tab_pplx10)
arpfc_tsne_pplx20 <- cbind(all_dist_pheno_fc[,1:9], arpfc_arranged_tsne_tab_pplx20)
arpd_tsne_pplx20 <- cbind(all_dist_pheno_delta[,1:9], arpd_arranged_tsne_tab_pplx20)
arpfc_tsne_pplx40 <- cbind(all_dist_pheno_fc[,1:9], arpfc_arranged_tsne_tab_pplx40)
arpd_tsne_pplx40 <- cbind(all_dist_pheno_delta[,1:9], arpd_arranged_tsne_tab_pplx40)
arpfc_tsne <- cbind(all_dist_pheno_fc[,1:9], arpfc_arranged_tsne_tab)arpdabs_tsne <- cbind(all_dist_pheno_absdelta[,1:9], arpdabs_arranged_tsne_tab)


arpd_tsne_male <- arpd_tsne%>%filter(.,Sex=="M")
arpd_tsne_female <- arpd_tsne%>%filter(.,Sex=="F")
#arpfc_tsne<-arpd_tsne%>%filter(.,Sex.x=='M')
#arpfc_tsne<-arpd_tsne%>%filter(.,Sex.x=='F')

#ratio delta tsne
pdf("schoeller_rtsne_ratio_delta_controls_arranged_pplx30.pdf", w=10, h=8)
ggplot(arpd_tsne[,10:11]) + 
  geom_point(aes(x=arpd_tsne[,10], 
                 y=arpd_tsne[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne[,7], fill=arpd_tsne[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_delta_sex_arranged_pplx30.pdf", w=10, h=8)
ggplot(arpd_tsne[,10:11]) + 
  geom_point(aes(x=arpd_tsne[,10], 
                 y=arpd_tsne[,11], 
                 color=arpd_tsne$Sex, 
                 fill=arpd_tsne$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne[,7], fill=arpd_tsne[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# age
pdf("schoeller_rtsne_ratio_delta_age-gradient_arranged_pplx30.pdf", w=10, h=8)
ggplot(arpd_tsne[,10:11]) + 
  geom_point(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], 
                 color=arpd_tsne$Age, fill=arpd_tsne$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne$Age, fill=arpd_tsne$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_ratio_delta_controls_arranged_pplx10.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx10[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx10[,10], 
                 y=arpd_tsne_pplx10[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpd_tsne_pplx10[,10], y=arpd_tsne_pplx10[,11], label=arpd_tsne_pplx10$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx10[,7], fill=arpd_tsne_pplx10[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_delta_sex_arranged_pplx10.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx10[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx10[,10], 
                 y=arpd_tsne_pplx10[,11], 
                 color=arpd_tsne_pplx10$Sex, 
                 fill=arpd_tsne_pplx10$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx10[,7], fill=arpd_tsne_pplx10[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# age
pdf("schoeller_rtsne_ratio_delta_age-gradient_arranged_pplx10.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx10[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx10[,10], y=arpd_tsne_pplx10[,11], 
                 color=arpd_tsne_pplx10$Age, fill=arpd_tsne_pplx10$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx10$Age, fill=arpd_tsne_pplx10$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_ratio_delta_controls_arranged_pplx20.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx20[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx20[,10], 
                 y=arpd_tsne_pplx20[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpd_tsne_pplx20[,10], y=arpd_tsne_pplx20[,11], label=arpd_tsne_pplx20$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx20[,7], fill=arpd_tsne_pplx20[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_delta_sex_arranged_pplx20.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx20[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx20[,10], 
                 y=arpd_tsne_pplx20[,11], 
                 color=arpd_tsne_pplx20$Sex, 
                 fill=arpd_tsne_pplx20$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx20[,7], fill=arpd_tsne_pplx20[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# age
pdf("schoeller_rtsne_ratio_delta_age-gradient_arranged_pplx20.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx20[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx20[,10], y=arpd_tsne_pplx20[,11], 
                 color=arpd_tsne_pplx20$Age, fill=arpd_tsne_pplx20$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx20$Age, fill=arpd_tsne_pplx20$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_ratio_delta_controls_arranged_pplx40.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx40[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx40[,10], 
                 y=arpd_tsne_pplx40[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpd_tsne_pplx40[,10], y=arpd_tsne_pplx40[,11], label=arpd_tsne_pplx40$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx40[,7], fill=arpd_tsne_pplx40[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_delta_sex_arranged_pplx40.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx40[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx40[,10], 
                 y=arpd_tsne_pplx40[,11], 
                 color=arpd_tsne_pplx40$Sex, 
                 fill=arpd_tsne_pplx40$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx40[,7], fill=arpd_tsne_pplx40[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()
# age
pdf("schoeller_rtsne_ratio_delta_age-gradient_arranged_pplx40.pdf", w=10, h=8)
ggplot(arpd_tsne_pplx40[,10:11]) + 
  geom_point(aes(x=arpd_tsne_pplx40[,10], y=arpd_tsne_pplx40[,11], 
                 color=arpd_tsne_pplx40$Age, fill=arpd_tsne_pplx40$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx40$Age, fill=arpd_tsne_pplx40$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

cluster.tsne.arpd.1.2 <- Mclust(arpd_tsne[,10:11], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpd.1.2.BIC.pdf", w=6, h=4)
plot(cluster.tsne.arpd.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpd.1.2.density.pdf", w=4, h=4)
plot(cluster.tsne.arpd.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpd.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.tsne.arpd.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpd.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.tsne.arpd.1.2, what = "uncertainty")
dev.off()


#all ratio pheno fc coloring just the controls
ggplot(arpfc_tsne[,10:11]) + 
  geom_point(aes(x=arpfc_tsne[,10], 
                 y=arpfc_tsne[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpfc_tsne[,10], y=arpfc_tsne[,11], label=arpfc_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of FC of Ratio Distances", x="tSNE dim1", y="tSNE dim2") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))

#ratio delta sex
ggplot(arpfc_tsne[,10:11]) + 
  geom_point(aes(x=arpfc_tsne[,10], 
                 y=arpfc_tsne[,11], 
                 color=arpfc_tsne$Sex, 
                 fill=arpfc_tsne$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of FC of Ratio Distances", x="tSNE dim1", y="tSNE dim2") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

# age
ggplot(arpfc_tsne[,10:11]) + 
  geom_point(aes(x=arpfc_tsne[,10], y=arpfc_tsne[,11], 
                 color=arpfc_tsne$Age, fill=arpfc_tsne$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpfc_tsne[,10], y=arpfc_tsne[,11], label=arpfc_tsne[,3])) +
  labs(title = "t-SNE of FC of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne$Age, fill=arpfc_tsne$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

pdf("schoeller_rtsne_ratio_fc_controls_arranged_pplx10.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx10[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx10[,10], 
                 y=arpfc_tsne_pplx10[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpfc_tsne_pplx10[,10], y=arpfc_tsne_pplx10[,11], label=arpfc_tsne_pplx10$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx10[,7], fill=arpfc_tsne_pplx10[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_fc_sex_arranged_pplx10.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx10[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx10[,10], 
                 y=arpfc_tsne_pplx10[,11], 
                 color=arpfc_tsne_pplx10$Sex, 
                 fill=arpfc_tsne_pplx10$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx10[,7], fill=arpd_tsne_pplx10[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# age
pdf("schoeller_rtsne_ratio_fc_age-gradient_arranged_pplx10.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx10[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx10[,10], y=arpfc_tsne_pplx10[,11], 
                 color=arpfc_tsne_pplx10$Age, fill=arpfc_tsne_pplx10$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx10$Age, fill=arpfc_tsne_pplx10$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_ratio_fc_controls_arranged_pplx20.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx20[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx20[,10], 
                 y=arpfc_tsne_pplx20[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpfc_tsne_pplx20[,10], y=arpfc_tsne_pplx20[,11], label=arpfc_tsne_pplx20$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx20[,7], fill=arpfc_tsne_pplx20[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_fc_sex_arranged_pplx20.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx20[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx20[,10], 
                 y=arpfc_tsne_pplx20[,11], 
                 color=arpfc_tsne_pplx20$Sex, 
                 fill=arpfc_tsne_pplx20$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne_pplx20[,7], fill=arpd_tsne_pplx20[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# age
pdf("schoeller_rtsne_ratio_fc_age-gradient_arranged_pplx20.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx20[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx20[,10], y=arpfc_tsne_pplx20[,11], 
                 color=arpfc_tsne_pplx20$Age, fill=arpfc_tsne_pplx20$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx20$Age, fill=arpfc_tsne_pplx20$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

pdf("schoeller_rtsne_ratio_fc_controls_arranged_pplx40.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx40[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx40[,10], 
                 y=arpfc_tsne_pplx40[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpfc_tsne_pplx40[,10], y=arpfc_tsne_pplx40[,11], label=arpfc_tsne_pplx40$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx40[,7], fill=arpfc_tsne_pplx40[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))
dev.off()

#ratio delta sex
pdf("schoeller_rtsne_ratio_fc_sex_arranged_pplx40.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx40[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx40[,10], 
                 y=arpfc_tsne_pplx40[,11], 
                 color=arpfc_tsne_pplx40$Sex, 
                 fill=arpfc_tsne_pplx40$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx40[,7], fill=arpfc_tsne_pplx40[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()
# age
pdf("schoeller_rtsne_ratio_fc_age-gradient_arranged_pplx40.pdf", w=10, h=8)
ggplot(arpfc_tsne_pplx40[,10:11]) + 
  geom_point(aes(x=arpfc_tsne_pplx40[,10], y=arpfc_tsne_pplx40[,11], 
                 color=arpfc_tsne_pplx40$Age, fill=arpfc_tsne_pplx40$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,10], y=arpd_tsne[,11], label=arpd_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpfc_tsne_pplx40$Age, fill=arpfc_tsne_pplx40$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

#Mclust w/ tSNE dimensions
cluster.tsne.arpfc.1.2 <- Mclust(arpd_tsne[,10:11], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpfc.1.2.BIC.pdf", w=6, h=4)
plot(cluster.tsne.arpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.tsne.arpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.tsne.arpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpfc.1.2.uncertainty.pdf", w=4, h=4)
plot(cluster.tsne.arpfc.1.2, what = "uncertainty")
dev.off()

# ratio absdelta 

ggplot(arpdabs_tsne[,10:11]) + 
  geom_point(aes(x=arpdabs_tsne[,10], 
                 y=arpdabs_tsne[,11],
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=arpdabs_tsne[,10], y=arpdabs_tsne[,11], label=arpdabs_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpdabs_tsne[,7], fill=arpdabs_tsne[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))

#ratio delta sex
ggplot(arpdabs_tsne[,10:11]) + 
  geom_point(aes(x=arpdabs_tsne[,10], 
                 y=arpdabs_tsne[,11], 
                 color=arpdabs_tsne$Sex, 
                 fill=arpdabs_tsne$Sex), 
             alpha = 0.8, 
             size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpd_tsne[,8], y=arpd_tsne[,9], label=arpd_tsne$Co_twins_ID)) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpd_tsne[,7], fill=arpd_tsne[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

# age
ggplot(arpdabs_tsne[,10:11]) + 
  geom_point(aes(x=arpdabs_tsne[,10], y=arpdabs_tsne[,11], 
                 color=arpdabs_tsne$Age, fill=arpdabs_tsne$Age), 
             alpha = 0.8, size = 15, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                               24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  #geom_text(aes(x=arpdabs_tsne[,10], y=arpdabs_tsne[,11], label=arpdabs_tsne[,3])) +
  labs(title = "t-SNE of Delta of Ratio Distances", x="tSNE dim1", y="tSNE dim2", color=arpdabs_tsne$Age, fill=arpdabs_tsne$Age) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

cluster.tsne.arpdabs.1.2 <- Mclust(arpdabs_tsne[,10:11], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpdabs.1.2.BIC.pdf", w=6, h=4)
plot(cluster.tsne.arpdabs.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpdabs.1.2.density.pdf", w=4, h=4)
plot(cluster.tsne.arpdabs.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpdabs.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.tsne.arpdabs.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.tsne.arpdabs.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.tsne.arpdabs.1.2, what = "uncertainty")
dev.off()

###### UMAP #######

library(umap)

library(M3C)

library(uwot)

t <- theme(plot.title = element_text(hjust = 0.5, size = 25), 
           axis.title.x = element_text(size = 20),
           axis.title.y = element_text(size = 20),  
           axis.text.x= element_text(size = 20), 
           axis.text.y= element_text(size = 20),
           legend.position = c(0.15,0.80),
           legend.title = element_text(size = 15), 
           legend.text = element_text(size = 15))

# dist delta

umap_adpd <- umap(all_dist_pheno_delta[,c(10:360)])
umap_adpd<-as.data.frame(umap_adpd)

pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_delta_controls_arranged.pdf", w=10, h=8)
ggplot(umap_adpd) + 
  geom_point(aes(x=V1, y=V2, 
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_delta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t+
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()

pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_delta_sex_arranged.pdf", w=10, h=8)
ggplot(umap_adpd) + 
  geom_point(aes(x=V1, y=V2, 
                 color=all_dist_pheno_delta$Sex_Diff,
                 fill=all_dist_pheno_delta$Sex_Diff),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_delta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t +
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()



pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_delta_age-gradient_arranged.pdf", w=10, h=8)
ggplot(umap_adpd) + 
  geom_point(aes(x=V1, y=V2, 
                 color=all_dist_pheno_delta$Age,
                 fill=all_dist_pheno_delta$Age),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_delta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + 
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) +
  t + theme(legend.position = "none")
dev.off()

cluster.umap.adpd.1.2 <- Mclust(umap_adpd[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpd.1.2.BIC.pdf", w=6, h=4)
plot(cluster.umap.adpd.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpd.1.2.density.pdf", w=4, h=4)
plot(cluster.umap.adpd.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpd.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.umap.adpd.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpd.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.umap.adpd.1.2, what = "uncertainty")
dev.off()

# dist absdelta

umap_adpdabs <- umap(all_dist_pheno_absdelta[,c(10:360)])
umap_adpdabs<-as.data.frame(umap_adpdabs)
setwd("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/")
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_absdelta_controls_arranged.pdf", w=10, h=8)
ggplot(umap_adpdabs) + 
  geom_point(aes(x=V1, y=V2, 
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_absdelta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t+
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()

pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_absdelta_sex_arranged.pdf", w=10, h=8)
ggplot(umap_adpdabs) + 
  geom_point(aes(x=V1, y=V2, 
                 color=all_dist_pheno_absdelta$Sex_Diff,
                 fill=all_dist_pheno_absdelta$Sex_Diff),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_absdelta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t +
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()



pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_absdelta_age-gradient_arranged.pdf", w=10, h=8)
ggplot(umap_adpdabs) + 
  geom_point(aes(x=V1, y=V2, 
                 color=all_dist_pheno_absdelta$Age,
                 fill=all_dist_pheno_absdelta$Age),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_absdelta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + 
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) +
  t + theme(legend.position = "none")
dev.off()

cluster.umap.adpdabs.1.2 <- Mclust(umap_adpdabs[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpdabs.1.2.BIC.pdf", w=6, h=4)
plot(cluster.umap.adpdabs.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpdabs.1.2.density.pdf", w=4, h=4)
plot(cluster.umap.adpdabs.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpdabs.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.umap.adpdabs.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpdabs.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.umap.adpdabs.1.2, what = "uncertainty")
dev.off()

# dist fc
umap_adpfc <- umap(all_dist_pheno_fc[,c(10:360)])
umap_adpfc<-as.data.frame(umap_adpfc)

pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_fc_controls_arranged.pdf", w=10, h=8)
ggplot(umap_adpfc) + 
  geom_point(aes(x=V1, y=V2, 
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")),
             alpha = 0.8, size = 15, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_fc$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t+
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()

pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_fc_sex_arranged.pdf", w=10, h=8)
ggplot(umap_adpfc) + 
  geom_point(aes(x=V1, y=V2, 
                 color=all_dist_pheno_fc$Sex_Diff,
                 fill=all_dist_pheno_fc$Sex_Diff),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_fc$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t +
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()



pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-UMAP/schoeller_umap_dist_fc_age-gradient_arranged.pdf", w=10, h=8)
ggplot(umap_adpfc) + 
  geom_point(aes(x=V1, y=V2, 
                 color=all_dist_pheno_fc$Age,
                 fill=all_dist_pheno_fc$Age),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              25,25,25,25,25)) + 
  geom_text(aes(x=V1, y=V2, label=all_dist_pheno_fc$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + 
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) +
  t + theme(legend.position = "none")
dev.off()

cluster.umap.adpfc.1.2 <- Mclust(umap_adpfc[,1:2], prior = priorControl())
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpfc.1.22.BIC.pdf", w=6, h=4)
plot(cluster.umap.adpfc.1.2, what = "BIC")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpfc.1.2.density.pdf", w=4, h=4)
plot(cluster.umap.adpfc.1.2, what = "density")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpfc.1.2.plot.classification.pdf", w=4, h=4)
plot(cluster.umap.adpfc.1.2, what = "classification")
dev.off()
pdf("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-mclust/cluster.umap.adpfc.1.2.plot.uncertainty.pdf", w=4, h=4)
plot(cluster.umap.adpfc.1.2, what = "uncertainty")
dev.off()

#for umap plots of ratio dfs, go to lab server and use node 

goodnames=colnames(dist_matrix_tot_scaled)[4:8]
plot_list = list()
for (i in 1:length(goodnames)){
  featurename=goodnames[i]
  p <- ggplot(umap) + 
    geom_point(aes(x=V1, y=V2,color=dist_matrix_tot_scaled[,3+i], fill=dist_matrix_tot_scaled[,3+i]),
               alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                24,24,24,24,24,24,24,24,24,24,
                                                22,22,22,22,22,22,22,22,22,22)) + 
    #geom_text(aes(x=V1, y=V2, label=dist_table$Individual_ID), size = 7) +
    labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y", color=featurename) +
    theme_classic() +t + guides(fill = FALSE)
  plot_list[[i]]=p
}

filename="UMAP_total_individuals_"
for (i in 1:length(goodnames)){
  featurename=goodnames[i]
  newfilename=paste(filename,featurename, sep="")
  newfilename=paste(newfilename,".pdf", sep="")
  ggsave(newfilename, plot= plot_list[[i]], device='pdf', width = 10, height= 8)
}


pdf("UMAP_total_individuals_Age.pdf", w=10, h=8)
ggplot(umap) + 
  geom_point(aes(x=V1, y=V2,color=as.numeric(dist_matrix_tot_scaled$Age)),
             alpha = 0.8, size = 20, shape= c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,
                                              19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,
                                              19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,
                                              19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,
                                              19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,
                                              19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,
                                              17,17,17,17,17,17,17,17,17,17,
                                              15,15,15,15,15,15,15,15,15,15)) + 
  #geom_text(aes(x=V1, y=V2, label=dist_table$Individual_ID), size = 7) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic()+t+
  scale_color_gradient2(low="white", mid="blue",high="red", name = "AGE")
dev.off()

#ratio matrix tsne unranked/uncorrelated
ggplot(umap_arpd[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_delta_nops_tsne[,8], 
                 y=all_ratio_pheno_delta_nops_tsne[,9], 
                 color=all_ratio_pheno_delta_nops_tsne[,7], 
                 fill=all_ratio_pheno_delta_nops_tsne[,7]), 
             alpha = 0.3, 
             size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                  21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,
                                  24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_delta_nops_tsne[,8], y=all_ratio_pheno_delta_nops_tsne[,9], label=all_ratio_pheno_delta_nops_tsne[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_delta_nops_tsne[,7], fill=all_ratio_pheno_delta_nops_tsne[,7]) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 



umap_arpfc <- umap(all_ratio_pheno_fc[,c(2:61426)])
umap_arpfc<-as.data.frame(umap_arpfc)
pdf("UMAP_arpfc.pdf", w=10, h=8)
ggplot(umap_arpfc) + 
  geom_point(aes(x=V1, y=V2, 
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none")),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              22,22,22,22,22,23,23)) + 
  geom_text(aes(x=V1, y=V2, label=all_ratio_pheno_fc$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t+
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()

pdf("umap_arpfc_noName.pdf", w=10, h=8)
ggplot(umap_arpfc) + 
  geom_point(aes(x=V1, y=V2, 
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none")),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              22,22,22,22,22,23,23)) + 
  #geom_text(aes(x=V1, y=V2, label=dist_matrix_tot_scaled$Individual_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() +t+
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()

pdf("umap_arpfc_noName_TWINS.pdf", w=10, h=8)
ggplot(umap_arpfc) + 
  geom_point(aes(x=V1, y=V2, 
                 color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "none","none","none","none","none","none","none","none","none","none","none","none","none",
                         "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                         "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"),
                 fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "none","none","none","none","none","none","none","none","none","none","none","none","none",
                        "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                        "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none")),
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              22,22,22,22,22,23,23)) + 
  #geom_text(aes(x=V1, y=V2, label=dist_matrix_tot_scaled$Individual_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() +t+
  scale_color_manual(values=c("#0727f5","#76EEC6" ,"#FF0000", "#FF7F24" , "#615292", "#AEC960", "#EE1289", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c("#0727f5","#76EEC6","#FF0000" , "#FF7F24" ,"#615292", "#AEC960", "#EE1289","white"), guide=FALSE)
dev.off()


##########Seurat########
library(Seurat)

scaled_all_dist_pheno_delta <- scale(all_dist_pheno_delta[,10:360])

pca_all_dist_pheno_delta <- RunPCA(scaled_all_dist_pheno_delta, verbose = FALSE)

#pca_all_dist_pheno_delta <- prcomp(scaled_all_dist_pheno_delta)

pca_all_dist_pheno_delta <- as.data.frame(pca_all_dist_pheno_delta[["sdev"]])

ElbowPlot(pca_all_dist_pheno_delta@stdev)



######### PCA 2 tSNE #######

library(Rtsne)
library(tidyverse)

# dist delta

scaled_all_dist_pheno_delta <- scale(all_dist_pheno_delta[,10:360])

pca_all_dist_pheno_delta <- prcomp(scaled_all_dist_pheno_delta)

pca_all_dist_pheno_delta_scores <- as.data.frame(pca_all_dist_pheno_delta[["x"]])

pca_all_dist_pheno_delta_scores <- cbind(all_dist_pheno_delta[1:2], pca_all_dist_pheno_delta_scores)

tsne <- Rtsne(pca_all_dist_pheno_delta_scores[,3:85],check_duplicates = FALSE, perplexity=20, max_iter = 20000, theta = 0.0)

tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
setwd("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-tSNE")
pdf("schoeller_rtsne_pca_dist_delta_controls_arranged_pplx20.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# dist fc

scaled_all_dist_pheno_fc <- scale(all_dist_pheno_fc[,10:360])

pca_all_dist_pheno_fc <- prcomp(scaled_all_dist_pheno_fc)

pca_all_dist_pheno_fc_scores <- as.data.frame(pca_all_dist_pheno_fc[["x"]])

pca_all_dist_pheno_fc_scores <- cbind(all_dist_pheno_fc[1:2], pca_all_dist_pheno_fc_scores)

# pplx 20

tsne <- Rtsne(pca_all_dist_pheno_fc_scores[,3:85],check_duplicates = FALSE, perplexity=20, max_iter = 20000, theta = 0.0)

tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
setwd("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-tSNE")
pdf("schoeller_rtsne_pca_dist_fc_controls_arranged_pplx20.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

#pplx 25

tsne <- Rtsne(pca_all_dist_pheno_fc_scores[,3:85],check_duplicates = FALSE, perplexity=25, max_iter = 20000, theta = 0.0)

tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
setwd("/Users/Warren.Sink/github/facial-polyphenism/Plots/Schoeller-tSNE")
pdf("schoeller_rtsne_pca_dist_fc_controls_arranged_pplx25.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

# ratio delta

scaled_all_ratio_pheno_delta <- scale(all_ratio_pheno_delta[,10:61434])

pca_all_ratio_pheno_delta <- prcomp(scaled_all_ratio_pheno_delta)

pca_all_ratio_pheno_delta_scores <- as.data.frame(pca_all_ratio_pheno_delta[["x"]])

pca_all_ratio_pheno_delta_scores <- cbind(all_dist_pheno_delta[1:2], pca_all_ratio_pheno_delta_scores)

tsne <- Rtsne(pca_all_ratio_pheno_delta_scores[,3:85],check_duplicates = FALSE, perplexity=20, max_iter = 20000, theta = 0.0)

tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
pdf("schoeller_rtsne_pca_ratio_delta_controls_arranged_pplx20.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()

tsne <- Rtsne(pca_all_ratio_pheno_delta_scores[,3:85],check_duplicates = FALSE, perplexity=15, max_iter = 20000, theta = 0.0)

tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
pdf("schoeller_rtsne_pca_ratio_delta_controls_arranged_pplx15.pdf", w=10, h=8)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                                     "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"),
                                   fill=c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                          "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
                                          "#668d3c","#668d3c","#668d3c","#668d3c","#668d3c")), 
                               alpha = 0.3, size = 10, shape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                                                 24,24,24,24,24,25,25,25,25,25)) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=all_dist_pheno_delta$Co_twins_ID)) +
  labs(x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 
dev.off()



