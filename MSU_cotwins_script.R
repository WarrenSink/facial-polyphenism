setwd("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute")
library(tidyverse)
####### Phenotypic data######



####### Generate distances and ratios of distances x pheno######

#distance x ratios of distances

msu <- read.csv(file = 'msu_coord_nose_norm.csv', header = TRUE, stringsAsFactors = FALSE)
msu_controls <- read.csv(file = 'msu_ctrl_coord_nose_norm.csv', header = TRUE, stringsAsFactors = FALSE)

msu <- full_join(msu,msu_controls,by="X")

msu[1,] <- msu[1,] %>% substr(., 1, 5)

#msu_neg <- (msu[,c(149,149,151,151,153,153,155,155,157,157)])

msu_neg <- (msu[,c(149,149,151,151,153,153,155,155,157,157)]) %>% rename(.,X147.2 = X147,
                                                                         X149.2 = X149,
                                                                         X151.2 = X151,
                                                                         X153.2 = X153,
                                                                         X155.2 = X155)
msu <- cbind(msu, msu_neg)

msu_id <- msu[1,]

msu <- slice(msu, -1)

vect<-function(n){
  a=strsplit(strsplit(n,split=",")[[1]][1],split="[[]")[[1]][2]
  b=strsplit(strsplit(n,split=",")[[1]][2],split="[]]")[[1]][1]
  return(c(as.numeric(a),as.numeric(b)))
}

eud<-function(x1,y1,x2,y2){
  sqrt((x2-x1)**2+(y2-y1)**2)
}

L=list()
for (i in 1:(nrow(msu)-1)){
  for (j in (i+1):nrow(msu)){
    
    n1=msu[i,1]; n2=msu[j,1] 
    r=as.numeric()
    
    for (v in 2:dim(msu)[2]){
      im.att1=vect(msu[i,v]); im.att2=vect(msu[j,v]) 
      d1= eud(im.att1[1],im.att1[2],im.att2[1],im.att2[2])
      r=c(r,d1)
    }
    
    name=paste(n1,n2,sep="_")
    L[[name]]=r
    #cat(name,"\n")
  }
}

M=do.call(rbind,L)
write.csv(M,"msu_features_distance.csv")

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

Co_twins_ID <- NULL
for (i in 1:nrow(M_t)){
  counter = 0 + i
  counter2 = 0 + (i-300)
  counter3 = 0 + (i-310)
  if (counter < 301){
    if ((counter %% 2) == 0) Co_twins_ID[[i]] = (counter/2)
    else Co_twins_ID[[i]] = ((counter/2)+0.5)
  }
  else if (counter > 300){
    if (counter < 311){
      if ((counter %% 2) == 0){
        Co_twins_ID[[i]] = paste("pos", (counter2/2),sep="_")
      } 
      else Co_twins_ID[[i]] = paste("pos", ((counter2/2)+0.5), sep = "_")
    }
    else{
      if ((counter %% 2) == 0) Co_twins_ID[[i]] = paste("neg", (counter3/2),sep="_")
      else Co_twins_ID[[i]] = paste("neg", ((counter3/2)+0.5), sep = "_")
    }
  }
}

Individual_ID <- NULL
for (i in 1:nrow(M_t)){
  counter = 0 + i
  counter2 = 0 + (i-300)
  counter3 = 0 + (i-310)
  if (counter < 301){
    if ((counter %% 2) == 0) Individual_ID[[i]] = paste((counter/2),"B",sep="")
    else Individual_ID[[i]] = paste(((counter/2)+0.5), "A", sep = "")
  }
  else if (counter > 300){
    if (counter < 311){
      if ((counter %% 2) == 0) Individual_ID[[i]] = paste("pos",(counter2/2),"B",sep="_")
      else Individual_ID[[i]] = paste("pos",((counter2/2)+0.5), "A", sep = "_")
    }
    else{
      if ((counter %% 2) == 0) Individual_ID[[i]] = paste("neg",(counter3/2),"B",sep="_")
      else Individual_ID[[i]] = paste("neg", ((counter3/2)+0.5),"A", sep = "_")
    }
  }
}

all_dist_pheno <- data.frame(cbind(t(msu_id[-1]),Co_twins_ID,Individual_ID,M_t))
all_dist_pheno <- remove_rownames(all_dist_pheno)
colnames(all_dist_pheno)[1] <- "msu_id"
msu_all_dist_pheno <- all_dist_pheno
msu_all_dist_pheno[,1] <- as.character(msu_all_dist_pheno[,1])
msu_all_dist_pheno[,c(4:354)] <- lapply(msu_all_dist_pheno[,c(4:354)], function(x) as.numeric(as.character(x)))

all_ratio_pheno <- data.frame(cbind(t(msu_id[-1]),Co_twins_ID,Individual_ID,M.ratio_t))
all_ratio_pheno <- remove_rownames(all_ratio_pheno)
colnames(all_ratio_pheno)[1] <- "msu_id"
msu_all_ratio_pheno <- all_ratio_pheno
msu_all_ratio_pheno[,1] <- as.character(msu_all_ratio_pheno[,1])
msu_all_ratio_pheno[,c(4:354)] <- lapply(msu_all_ratio_pheno[,c(4:354)], function(x) as.numeric(as.character(x)))


rm(list = ls()[!ls() %in% c("msu_all_dist_pheno","msu_all_ratio_pheno")])

#####Delta/FC###########

#absdelta
CoIds= msu_all_dist_pheno[,2]


UniqIds=unique(CoIds)

count=1
L.delta<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:(length(idx)-1)){
    for (w in (v+1):length(idx)){
      i=idx[v]
      j=idx[w]
      nam= paste0(msu_all_dist_pheno[i,1],"__",msu_all_dist_pheno[j,1])
      nam2=paste0(msu_all_dist_pheno)
      delta<-abs((as.numeric(msu_all_dist_pheno[i,4:dim(msu_all_dist_pheno)[2]]) - as.numeric(msu_all_dist_pheno[j,4:dim(msu_all_dist_pheno)[2]])))
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

msu_all_dist_pheno_absdelta<-do.call(rbind,L.delta)
colnames(msu_all_dist_pheno_absdelta)<-colnames(msu_all_dist_pheno)[10:ncol(msu_all_dist_pheno)]


#dist delta/fc

CoIds= msu_all_dist_pheno[,2]


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
      
      nam= paste0(msu_all_dist_pheno[i,1],"__",msu_all_dist_pheno[j,1])
      ratio<-as.numeric(msu_all_dist_pheno[i,4:dim(msu_all_dist_pheno)[2]]) /as.numeric(msu_all_dist_pheno[j,4:dim(msu_all_dist_pheno)[2]])
      delta<-as.numeric(msu_all_dist_pheno[i,4:dim(msu_all_dist_pheno)[2]] - as.numeric(msu_all_dist_pheno[j,4:dim(msu_all_dist_pheno)[2]]))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

msu_all_dist_pheno_fc<-do.call(rbind,L.ratio)
colnames(msu_all_dist_pheno_fc)<-colnames(msu_all_dist_pheno)[10:ncol(msu_all_dist_pheno)]
msu_all_dist_pheno_delta<-do.call(rbind,L.delta)
colnames(msu_all_dist_pheno_delta)<-colnames(msu_all_dist_pheno)[10:ncol(msu_all_dist_pheno)]


#ratio delta

CoIds= msu_all_ratio_pheno[,2]


UniqIds=unique(CoIds)

count=1

L.delta<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:(length(idx)-1)){
    for (w in (v+1):length(idx)){
      i=idx[v]
      j=idx[w]
      
      nam= paste0(msu_all_ratio_pheno[i,1],"__",msu_all_ratio_pheno[j,1])
      delta<-abs(as.numeric(msu_all_ratio_pheno[i,4:dim(msu_all_ratio_pheno)[2]]) -as.numeric(msu_all_ratio_pheno[j,4:dim(msu_all_ratio_pheno)[2]]))
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

msu_all_ratio_pheno_delta<-do.call(rbind,L.delta)
colnames(msu_all_ratio_pheno_delta)<-colnames(msu_all_ratio_pheno)[4:ncol(msu_all_ratio_pheno)]

#ratio abs(delta)/fc

CoIds= msu_all_ratio_pheno[,2]


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
      
      nam= paste0(msu_all_ratio_pheno[i,1],"__",msu_all_ratio_pheno[j,1])
      ratio<-as.numeric(msu_all_ratio_pheno[i,4:dim(msu_all_ratio_pheno)[2]]) / as.numeric(msu_all_ratio_pheno[j,4:dim(msu_all_ratio_pheno)[2]])
      delta<-abs(as.numeric(all_ratio_pheno[i,4:dim(msu_all_ratio_pheno)[2]]) - as.numeric(msu_all_ratio_pheno[j,4:dim(msu_all_ratio_pheno)[2]]))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}
all_ratio_pheno_fc<-do.call(rbind,L.ratio)
colnames(all_ratio_pheno_fc)<-colnames(msu_all_ratio_pheno)[4:ncol(msu_all_ratio_pheno)]
all_ratio_pheno_absdelta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_absdelta)<-colnames(msu_all_ratio_pheno)[4:ncol(msu_all_ratio_pheno)]

#######Log2 and FC/Delta#######

#Log2FC twins tables
all_dist_pheno_log2 <- all_dist_pheno %>% group_by(Co_twins_ID)  %>%
  mutate_at(vars(names(.[,c(10:360)])), list(log2 = ~ log2(.))) %>%
  ungroup() %>% subset(., select = -c(10:360))
all_dist_pheno_abslog2 <- all_dist_pheno_log2 %>% group_by(Co_twins_ID)  %>%
  mutate_at(vars(names(.[,c(10:360)])), list(abslog2 = ~ abs(.))) %>%
  ungroup() %>% subset(., select = -c(10:360))
all_ratio_pheno_log2 <- all_ratio_pheno %>% group_by(Co_twins_ID)  %>%
  mutate_at(vars(names(.[,c(10:61434)])), list(log2 = ~ log2(.))) %>%
  ungroup() %>% subset(., select = -c(10:61434))
all_ratio_pheno_abslog2 <- all_ratio_pheno_log2 %>% group_by(Co_twins_ID)  %>%
  mutate_at(vars(names(.[,c(10:61434)])), list(abslog2 = ~ abs(.))) %>%
  ungroup() %>% subset(., select = -c(10:61434))

#log2 dist delta, fc, and absdelta

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
      ratio<-log2(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]])) / log2(as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]]))
      delta<-log2(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]])) - log2(as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]]))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_dist_pheno_log2_fc<-do.call(rbind,L.ratio)
colnames(all_dist_pheno_log2_fc)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]
all_dist_pheno_log2_delta<-do.call(rbind,L.delta)
colnames(all_dist_pheno_log2_delta)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]

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
      delta<-abs(log2(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]])) - log2(as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]])))
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_dist_pheno_log2_absdelta<-do.call(rbind,L.delta)
colnames(all_dist_pheno_log2_absdelta)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]

#abslog2 dist delta, fc, and absdelta

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
      ratio<-abs(log2(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]]))) / abs(log2(as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]])))
      delta<-abs(log2(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]]))) - abs(log2(as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]])))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_dist_pheno_abslog2_fc<-do.call(rbind,L.ratio)
colnames(all_dist_pheno_abslog2_fc)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]
all_dist_pheno_abslog2_delta<-do.call(rbind,L.delta)
colnames(all_dist_pheno_abslog2_delta)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]

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
      delta<-abs(abs(log2(as.numeric(all_dist_pheno[i,10:dim(all_dist_pheno)[2]])) - abs(log2(as.numeric(all_dist_pheno[j,10:dim(all_dist_pheno)[2]])))))
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_dist_pheno_abslog2_absdelta<-do.call(rbind,L.delta)
colnames(all_dist_pheno_abslog2_absdelta)<-colnames(all_dist_pheno)[10:ncol(all_dist_pheno)]

#log2 and abslog2 ratio

CoIds= all_ratio_pheno[,3]


UniqIds=unique(CoIds)

count=1
L.logtwo<-list()
L.abslogtwo<-list()

for (UniqId in UniqIds){
  
  idx<-as.numeric(which(CoIds==UniqId))
  
  for (v in 1:length(idx)){
    i=idx[v]
    nam= paste0(all_ratio_pheno[i,1])
    logtwo<-log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]))
    abslogtwo<-abs(log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]])))
    L.logtwo[[nam]]<-logtwo
    L.abslogtwo[[nam]]<-abslogtwo
    cat(nam,"\n")
  }
}
all_ratio_pheno_log2<-do.call(rbind,L.logtwo)
colnames(all_ratio_pheno_log2)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]
all_ratio_pheno_abslog2<-do.call(rbind,L.abslogtwo)
colnames(all_ratio_pheno_abslog2)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]

L.logtwo<-list()
L.abslogtwo<-list()

for (i in 1:length(all_ratio_pheno[,1])){
  nam= paste0(all_ratio_pheno[i,1])
  logtwo<-log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]))
  abslogtwo<-abs(log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]])))
  L.logtwo[[nam]]<-logtwo
  L.abslogtwo[[nam]]<-abslogtwo
  cat(nam,"\n")
}
all_ratio_pheno_log2<-do.call(rbind,L.logtwo)
colnames(all_ratio_pheno_log2)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]
all_ratio_pheno_abslog2<-do.call(rbind,L.abslogtwo)
colnames(all_ratio_pheno_abslog2)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]


#log2 ratio delta, fc

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
      ratio<-log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]])) / log2(as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]]))
      delta<-log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]])) - log2(as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]]))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_ratio_pheno_log2_fc<-do.call(rbind,L.ratio)
colnames(all_ratio_pheno_log2_fc)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]
all_ratio_pheno_log2_delta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_log2_delta)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]

#log2 absdelta

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
      delta<-abs(log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]])) - log2(as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]])))
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_ratio_pheno_log2_absdelta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_log2_absdelta)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]


#abslog2 ratio delta, fc, and absdelta

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
      ratio<-abs(log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]))) / abs(log2(as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]])))
      delta<-abs(log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]]))) - abs(log2(as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]])))
      L.ratio[[nam]]<-ratio
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_ratio_pheno_abslog2_fc<-do.call(rbind,L.ratio)
colnames(all_ratio_pheno_abslog2_fc)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]
all_ratio_pheno_abslog2_delta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_abslog2_delta)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]

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
      delta<-abs(abs(log2(as.numeric(all_ratio_pheno[i,10:dim(all_ratio_pheno)[2]])) - log2(as.numeric(all_ratio_pheno[j,10:dim(all_ratio_pheno)[2]]))))
      L.delta[[nam]]<-delta
      cat(nam,"\n")
    }
  }
}

all_ratio_pheno_abslog2_absdelta<-do.call(rbind,L.delta)
colnames(all_ratio_pheno_abslog2_absdelta)<-colnames(all_ratio_pheno)[10:ncol(all_ratio_pheno)]


####### getting rid of NAs + paired pheno x paired dfs########

#getting rid of NAs

all_dist_pheno_fc[sapply(all_dist_pheno_fc, is.infinite)] <- NA
all_dist_pheno_fc<-all_dist_pheno_fc[, colSums(is.na(all_dist_pheno_fc)) != nrow(all_dist_pheno_fc)]

all_dist_pheno_log2[sapply(all_dist_pheno_log2, is.infinite)] <- NA
all_dist_pheno_log2<-all_dist_pheno_log2[, colSums(is.na(all_dist_pheno_log2)) != nrow(all_dist_pheno_log2)]

all_dist_pheno_log2_fc[sapply(all_dist_pheno_log2_fc, is.infinite)] <- NA
all_dist_pheno_log2_fc<-all_dist_pheno_log2_fc[, colSums(is.na(all_dist_pheno_log2_fc)) != nrow(all_dist_pheno_log2_fc)]

all_dist_pheno_log2_delta[sapply(all_dist_pheno_log2_delta, is.infinite)] <- NA
all_dist_pheno_log2_delta<-all_dist_pheno_log2_delta[, colSums(is.na(all_dist_pheno_log2_delta)) != nrow(all_dist_pheno_log2)]

all_dist_pheno_log2_absdelta[sapply(all_dist_pheno_log2_absdelta, is.infinite)] <- NA
all_dist_pheno_log2_absdelta<-all_dist_pheno_log2_absdelta[, colSums(is.na(all_dist_pheno_log2_absdelta)) != nrow(all_dist_pheno_log2_absdelta)]

all_dist_pheno_abslog2[sapply(all_dist_pheno_abslog2, is.infinite)] <- NA
all_dist_pheno_abslog2<-all_dist_pheno_abslog2[, colSums(is.na(all_dist_pheno_abslog2)) != nrow(all_dist_pheno_abslog2)]

all_dist_pheno_abslog2_fc[sapply(all_dist_pheno_abslog2_fc, is.infinite)] <- NA
all_dist_pheno_abslog2_fc<-all_dist_pheno_abslog2_fc[, colSums(is.na(all_dist_pheno_abslog2_fc)) != nrow(all_dist_pheno_abslog2_fc)]

all_dist_pheno_abslog2_delta[sapply(all_dist_pheno_abslog2_delta, is.infinite)] <- NA
all_dist_pheno_abslog2_delta<-all_dist_pheno_abslog2_delta[, colSums(is.na(all_dist_pheno_abslog2_delta)) != nrow(all_dist_pheno_abslog2)]

all_dist_pheno_abslog2_absdelta[sapply(all_dist_pheno_abslog2_absdelta, is.infinite)] <- NA
all_dist_pheno_abslog2_absdelta<-all_dist_pheno_abslog2_absdelta[, colSums(is.na(all_dist_pheno_abslog2_absdelta)) != nrow(all_dist_pheno_abslog2_absdelta)]

all_ratio_pheno[sapply(all_ratio_pheno, is.infinite)] <- NA
all_ratio_pheno<-all_ratio_pheno[, colSums(is.na(all_ratio_pheno)) != nrow(all_ratio_pheno)]

all_ratio_pheno_delta[sapply(all_ratio_pheno_delta, is.infinite)] <- NA
all_ratio_pheno_delta<-all_ratio_pheno_delta[, colSums(is.na(all_ratio_pheno_delta)) != nrow(all_ratio_pheno_delta)]

all_ratio_pheno_absdelta[sapply(all_ratio_pheno_absdelta, is.infinite)] <- NA
all_ratio_pheno_absdelta<-all_ratio_pheno_absdelta[, colSums(is.na(all_ratio_pheno_absdelta)) != nrow(all_ratio_pheno_absdelta)]

all_ratio_pheno_fc[sapply(all_ratio_pheno_fc, is.infinite)] <- NA
all_ratio_pheno_fc<-all_ratio_pheno_fc[, colSums(is.na(all_ratio_pheno_fc)) != nrow(all_ratio_pheno_fc)]

all_ratio_pheno_log2[sapply(all_ratio_pheno_log2, is.infinite)] <- NA
all_ratio_pheno_log2<-all_ratio_pheno_log2[, colSums(is.na(all_ratio_pheno_log2)) != nrow(all_ratio_pheno_log2)]

all_ratio_pheno_log2_fc[sapply(all_ratio_pheno_log2_fc, is.infinite)] <- NA
all_ratio_pheno_log2_fc<-all_ratio_pheno_log2_fc[, colSums(is.na(all_ratio_pheno_log2_fc)) != nrow(all_ratio_pheno_log2_fc)]

all_ratio_pheno_log2_delta[sapply(all_ratio_pheno_log2_delta, is.infinite)] <- NA
all_ratio_pheno_log2_delta<-all_ratio_pheno_log2_delta[, colSums(is.na(all_ratio_pheno_log2_delta)) != nrow(all_ratio_pheno_log2_delta)]

all_ratio_pheno_log2_absdelta[sapply(all_ratio_pheno_log2_absdelta, is.infinite)] <- NA
all_ratio_pheno_log2_absdelta<-all_ratio_pheno_log2_absdelta[, colSums(is.na(all_ratio_pheno_log2_absdelta)) != nrow(all_ratio_pheno_log2_absdelta)]

all_ratio_pheno_abslog2[sapply(all_ratio_pheno_abslog2, is.infinite)] <- NA
all_ratio_pheno_abslog2<-all_ratio_pheno_abslog2[, colSums(is.na(all_ratio_pheno_abslog2)) != nrow(all_ratio_pheno_abslog2)]

all_ratio_pheno_abslog2_fc[sapply(all_ratio_pheno_abslog2_fc, is.infinite)] <- NA
all_ratio_pheno_abslog2_fc<-all_ratio_pheno_abslog2_fc[, colSums(is.na(all_ratio_pheno_abslog2_fc)) != nrow(all_ratio_pheno_abslog2_fc)]

all_ratio_pheno_abslog2_delta[sapply(all_ratio_pheno_abslog2_delta, is.infinite)] <- NA
all_ratio_pheno_abslog2_delta<-all_ratio_pheno_abslog2_delta[, colSums(is.na(all_ratio_pheno_abslog2_delta)) != nrow(all_ratio_pheno_abslog2_delta)]

all_ratio_pheno_abslog2_absdelta[sapply(all_ratio_pheno_abslog2_absdelta, is.infinite)] <- NA
all_ratio_pheno_abslog2_absdelta<-all_ratio_pheno_abslog2_absdelta[, colSums(is.na(all_ratio_pheno_abslog2_absdelta)) != nrow(all_ratio_pheno_abslog2_absdelta)]

#paired pheno x paired dfs#

pheno_data_paired = pheno_data[,c(3:9)]
Co_twins_ID <- c('1','1','1','10','11','12','13','14','15','16','17','18','18','18','18',
                 '18','18','19','2','20','21','22','23','23','23','24','24','24','24','24',
                 '24','25','26','27','28','29','3','30','31','32','33','34','35','36','37','38',
                 '39','4','40','41','42','43','44','45','46','47','48','49','5','50','51','52','53',
                 '54','6','7','8','8','8','8','8','8','9',
                 'pos1','pos2','pos3','pos4','pos5','neg1','neg2','neg3','neg4','neg5',"55","56")

all_dist_pheno_fc <- data.frame(all_dist_pheno_fc)
all_dist_pheno_fc <- cbind(Co_twins_ID, all_dist_pheno_fc)
all_dist_pheno_fc <- full_join(all_dist_pheno_fc, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_fc <- distinct(all_dist_pheno_fc,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_delta<-data.frame(all_dist_pheno_delta)
all_dist_pheno_delta<- cbind(Co_twins_ID, all_dist_pheno_delta)
all_dist_pheno_delta<-left_join(all_dist_pheno_delta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_delta <- distinct(all_dist_pheno_delta,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_absdelta<-data.frame(all_dist_pheno_absdelta)
all_dist_pheno_absdelta<- cbind(Co_twins_ID, all_dist_pheno_absdelta)
all_dist_pheno_absdelta<-left_join(all_dist_pheno_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_absdelta <- distinct(all_dist_pheno_absdelta,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_log2<-data.frame(all_dist_pheno_log2)
all_dist_pheno_abslog2<-data.frame(all_dist_pheno_abslog2)

all_dist_pheno_log2_delta<-data.frame(all_dist_pheno_log2_delta)
all_dist_pheno_log2_delta<- cbind(Co_twins_ID, all_dist_pheno_log2_delta)
all_dist_pheno_log2_delta<-left_join(all_dist_pheno_log2_delta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_log2_delta <- distinct(all_dist_pheno_log2_delta,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_log2_absdelta<-data.frame(all_dist_pheno_log2_absdelta)
all_dist_pheno_log2_absdelta<- cbind(Co_twins_ID, all_dist_pheno_log2_absdelta)
all_dist_pheno_log2_absdelta<-left_join(all_dist_pheno_log2_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_log2_absdelta <- distinct(all_dist_pheno_log2_absdelta,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_log2_fc<-data.frame(all_dist_pheno_log2_fc)
all_dist_pheno_log2_fc<- cbind(Co_twins_ID, all_dist_pheno_log2_fc)
all_dist_pheno_log2_fc<-left_join(all_dist_pheno_log2_fc, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_log2_fc <- distinct(all_dist_pheno_log2_fc,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_abslog2_delta<-data.frame(all_dist_pheno_abslog2_delta)
all_dist_pheno_abslog2_delta<- cbind(Co_twins_ID, all_dist_pheno_abslog2_delta)
all_dist_pheno_abslog2_delta<-left_join(all_dist_pheno_abslog2_delta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_abslog2_delta <- distinct(all_dist_pheno_abslog2_delta,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_abslog2_absdelta<-data.frame(all_dist_pheno_abslog2_absdelta)
all_dist_pheno_abslog2_absdelta<- cbind(Co_twins_ID, all_dist_pheno_abslog2_absdelta)
all_dist_pheno_abslog2_absdelta<-left_join(all_dist_pheno_abslog2_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_abslog2_absdelta <- distinct(all_dist_pheno_abslog2_absdelta,ChinBottom_LeftPupil,.keep_all= TRUE)

all_dist_pheno_abslog2_fc<-data.frame(all_dist_pheno_abslog2_fc)
all_dist_pheno_abslog2_fc<- cbind(Co_twins_ID, all_dist_pheno_abslog2_fc)
all_dist_pheno_abslog2_fc<-left_join(all_dist_pheno_abslog2_fc, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_abslog2_fc <- distinct(all_dist_pheno_abslog2_fc,ChinBottom_LeftPupil,.keep_all= TRUE)

all_ratio_pheno_fc <- data.frame(all_ratio_pheno_fc)
all_ratio_pheno_fc <- cbind(Co_twins_ID, all_ratio_pheno_fc)
all_ratio_pheno_fc <- left_join(all_ratio_pheno_fc, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_fc <- distinct(all_ratio_pheno_fc,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_delta<-data.frame(all_ratio_pheno_delta)
all_ratio_pheno_delta<- cbind(Co_twins_ID, all_ratio_pheno_delta)
all_ratio_pheno_delta<-left_join(all_ratio_pheno_delta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_delta <- distinct(all_ratio_pheno_delta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_absdelta<-data.frame(all_ratio_pheno_absdelta)
all_ratio_pheno_absdelta<- cbind(Co_twins_ID, all_ratio_pheno_absdelta)
all_ratio_pheno_absdelta<-left_join(all_ratio_pheno_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_absdelta <- distinct(all_ratio_pheno_absdelta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_log2<-data.frame(all_ratio_pheno_log2)
all_ratio_pheno_abslog2<-data.frame(all_ratio_pheno_abslog2)

all_ratio_pheno_log2_delta<-data.frame(all_ratio_pheno_log2_delta)
all_ratio_pheno_log2_delta<- cbind(Co_twins_ID, all_ratio_pheno_log2_delta)
all_ratio_pheno_log2_delta<-left_join(all_ratio_pheno_log2_delta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_log2_delta <- distinct(all_ratio_pheno_log2_delta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_log2_absdelta<-data.frame(all_ratio_pheno_log2_absdelta)
all_ratio_pheno_log2_absdelta<- cbind(Co_twins_ID, all_ratio_pheno_log2_absdelta)
all_ratio_pheno_log2_absdelta<-left_join(all_ratio_pheno_log2_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_log2_absdelta <- distinct(all_ratio_pheno_log2_absdelta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_log2_fc<-data.frame(all_ratio_pheno_log2_fc)
all_ratio_pheno_log2_fc<- cbind(Co_twins_ID, all_ratio_pheno_log2_fc)
all_ratio_pheno_log2_fc<-left_join(all_ratio_pheno_log2_fc, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_log2_fc <- distinct(all_ratio_pheno_log2_fc,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_abslog2_fc<-data.frame(all_ratio_pheno_abslog2_fc)
all_ratio_pheno_abslog2_fc<- cbind(Co_twins_ID, all_ratio_pheno_abslog2_fc)
all_ratio_pheno_abslog2_fc<-left_join(all_ratio_pheno_abslog2_fc, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_abslog2_fc <- distinct(all_ratio_pheno_abslog2_fc,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_abslog2_delta<-data.frame(all_ratio_pheno_abslog2_delta)
all_ratio_pheno_abslog2_delta<-cbind(Co_twins_ID, all_ratio_pheno_abslog2_delta)
all_ratio_pheno_abslog2_delta<-left_join(all_ratio_pheno_abslog2_delta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_abslog2_delta<-distinct(all_ratio_pheno_abslog2_delta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

all_ratio_pheno_abslog2_absdelta<-data.frame(all_ratio_pheno_abslog2_absdelta)
all_ratio_pheno_abslog2_absdelta<-cbind(Co_twins_ID, all_ratio_pheno_abslog2_absdelta)
all_ratio_pheno_abslog2_absdelta<-left_join(all_ratio_pheno_abslog2_absdelta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_abslog2_absdelta<-distinct(all_ratio_pheno_abslog2_absdelta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)

rm(list = ls()[!ls() %in% c("pheno_data","all_dist_pheno","all_ratio_pheno",'all_dist_pheno_delta','all_dist_pheno_absdelta',
                            'all_dist_pheno_fc','all_ratio_pheno_delta','all_ratio_pheno_absdelta','all_ratio_pheno_fc',
                            "all_dist_pheno_log2","all_dist_pheno_abslog2", "all_dist_pheno_log2_fc",'all_dist_pheno_log2_absdelta',
                            'all_dist_pheno_log2_delta','all_dist_pheno_abslog2_fc','all_dist_pheno_abslog2_delta','all_dist_pheno_abslog2_delta',
                            "all_ratio_pheno_log2","all_ratio_pheno_abslog2",'all_ratio_pheno_log2_fc',
                            'all_ratio_pheno_log2_delta','all_ratio_pheno_log2_absdelta',"all_ratio_pheno_abslog2_fc","all_ratio_pheno_abslog2_delta",
                            'all_ratio_pheno_abslog2_absdelta')])

#save only the numeric columns for ipython tsne

all_ratio_pheno_delta_nops<-all_ratio_pheno_delta[,2:61426]
all_ratio_pheno_absdelta_nops<-all_ratio_pheno_absdelta[,2:61426]
all_ratio_pheno_fc_nops<-all_ratio_pheno_fc[,2:61426]
all_ratio_pheno_abslog2_delta_nops<-all_ratio_pheno_abslog2_delta[,2:61426]
all_ratio_pheno_abslog2_absdelta_nops<-all_ratio_pheno_abslog2_absdelta[,2:61426]
all_ratio_pheno_abslog2_fc_nops<-all_ratio_pheno_abslog2_fc[,2:61426]
all_ratio_pheno_log2_delta_nops<-all_ratio_pheno_log2_delta[,2:61426]
all_ratio_pheno_log2_fc_nops<-all_ratio_pheno_log2_fc[,2:61426]

save(all_ratio_pheno_delta_nops, file="all_ratio_pheno_delta_nops.Rdata")
save(all_ratio_pheno_absdelta_nops, file="all_ratio_pheno_absdelta_nops.Rdata")
save(all_ratio_pheno_fc_nops, file="all_ratio_pheno_fc_nops.Rdata")
save(all_ratio_pheno_abslog2_delta_nops, file="all_ratio_pheno_abslog2_delta_nops.Rdata")
save(all_ratio_pheno_abslog2_absdelta_nops, file="all_ratio_pheno_abslog2_absdelta_nops.Rdata")
save(all_ratio_pheno_abslog2_fc_nops, file="all_ratio_pheno_abslog2_fc_nops.Rdata")
save(all_ratio_pheno_log2_delta_nops, file="all_ratio_pheno_log2_delta_nops.Rdata")
save(all_ratio_pheno_log2_fc_nops, file="all_ratio_pheno_log2_fc_nops.Rdata")

write_csv(all_ratio_pheno_delta_nops,path = "Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/alpd_nops.csv")

########hpc-ratio-tsne trampoline########

tab_all_ratio_pheno_absdelta<-read.table(file="/Volumes/projects_secondary/mnp/warren/test/tsne_all_ratio_delta.tab", sep ="\t", header=F)
tab_all_ratio_pheno_fc<-read.table(file="/Volumes/projects_secondary/mnp/warren/test/tsne_all_ratio_fc.tab", sep ="\t", header=F)

all_ratio_pheno_delta_tsne<-cbind(all_dist_pheno_fc[c(1,406:411)],tab_all_ratio_pheno_absdelta)
all_ratio_pheno_fc_tsne<-cbind(all_dist_pheno_fc[c(1,407:413)],tab_all_ratio_pheno_fc)

all_ratio_pheno_fc_tsne_male<-all_ratio_pheno_fc_tsne%>%filter(.,Sex.y=="M")
all_ratio_pheno_fc_tsne_female<-all_ratio_pheno_fc_tsne%>%filter(.,Sex.y=="F")
all_ratio_pheno_delta_tsne_male<-all_ratio_pheno_delta_tsne%>%filter(.,Sex.y=='M')
all_ratio_pheno_delta_tsne_female<-all_ratio_pheno_delta_tsne%>%filter(.,Sex.y=='F')

rm(list = ls()[!ls() %in% c("pheno_data","all_dist_pheno","all_ratio_pheno",'all_dist_pheno_delta',
                            'all_dist_pheno_fc','all_ratio_pheno_delta','all_ratio_pheno_fc',
                            "all_dist_pheno_log2","all_dist_pheno_abslog2", "all_dist_pheno_log2_fc",
                            'all_dist_pheno_log2_delta','all_dist_pheno_abslog2_fc','all_dist_pheno_abslog2_delta',
                            'all_dist_pheno_abslog2_absdelta','all_dist_pheno_log2_absdelta',
                            "all_ratio_pheno_log2","all_ratio_pheno_abslog2",'all_ratio_pheno_log2_fc',
                            'all_ratio_pheno_log2_delta',"all_ratio_pheno_abslog2_fc","all_ratio_pheno_abslog2_delta",
                            "all_ratio_pheno_abslog2_absdelta","all_ratio_pheno_delta_tsne","all_ratio_pheno_fc_tsne",
                            "all_ratio_pheno_fc_tsne_male","all_ratio_pheno_fc_tsne_female","all_ratio_pheno_delta_tsne_male",
                            "all_ratio_pheno_delta_tsne_female","schoeller_similarity_scores_pheno","rd_0.Rdata",
                            "similarity_scores_pheno_aldi","all_ratio_pheno_delta_nops_tsne","all_ratio_pheno_fc_nops_tsne",
                            "all_ratio_pheno_delta_nops_tsne_male","all_ratio_pheno_delta_nops_tsne_female",
                            "all_ratio_pheno_fc_nops_tsne_male","all_ratio_pheno_fc_nops_tsne_female")])

rm(list = ls()[!ls() %in% c("all_ratio_pheno_delta_nops_tsne","all_ratio_pheno_fc_nops_tsne",
                            "all_ratio_pheno_delta_nops_tsne_male","all_ratio_pheno_delta_nops_tsne_female",
                            "all_ratio_pheno_fc_nops_tsne_male","all_ratio_pheno_fc_nops_tsne_female")])
####Elimination of highly correlative or retention of most covariate elements#####
##correlation matrix##
write.csv(all_ratio_pheno_delta[,2:82102],"schoeller_all_ratio_pheno_delta.csv")
write.csv(all_ratio_pheno_fc[,2:82407],"schoeller_all_ratio_pheno_fc.csv")
write.csv(all_ratio_pheno_log2_delta[,2:81407],"schoeller_all_ratio_pheno_log2_delta.csv")
write.csv(all_ratio_pheno_log2_fc[,2:81351],"schoeller_all_ratio_pheno_log2_fc.csv")
write.csv(all_ratio_pheno_abslog2_delta[,2:81407],"schoeller_all_ratio_pheno_abslog2_delta.csv")
write.csv(all_ratio_pheno_abslog2_fc[,2:81346],"schoeller_all_ratio_pheno_abslog2_fc.csv")


library(corrplot)
correlation_matrix_ratio_delta <- cor(all_ratio_pheno_delta, all_ratio_pheno_delta, method = "pearson" )
correlation_matrix_ratio_fc <- cor(all_ratio_pheno_fc, all_ratio_pheno_fc, method = "pearson" )
#correlation_matrix <- cor(data[c(6)] , data[c(11:50)] , use="complete.obs" , method = "pearson" ) ###this is to skip NA data###

library("Hmisc") ##this is to obtain pval associated to correlation matrix##
all_ratio_pheno_delta_mat<-as.matrix(all_ratio_pheno_delta)
all_ratio_pheno_fc_mat<-as.matrix(all_ratio_pheno_fc)
correlation_matrix_ratio_delta_2 <- rcorr(all_ratio_pheno_delta_mat , all_ratio_pheno_delta_mat , type  = "pearson" )
correlation_matrix_ratio_fc_2 <- rcorr(as.matrix(all_ratio_pheno_fc) , as.matrix(data[all_ratio_pheno_fc]) , type  = "pearson" )
correlation_matrix_ratio_delta_2$r ##Extract the correlation coefficients##
correlation_matrix_ratio_delta_2$P ##Extract p-values##
correlation_matrix_ratio_fc_2$r
correlation_matrix_ratio_fc_2$P

##this is to format the correlation matrix##

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
} 
flattenCorrMatrix(correlation_matrix_2$r, correlation_matrix_2$P)
flatten_corr_matrix <- flattenCorrMatrix(correlation_matrix_2$r, correlation_matrix_2$P)
View(flatten_corr_matrix)
##the following is for correlation matrix with corrplot##
library(corrplot)
correlation_matrix_total <- cor(data[c(6,11:50)] , data[c(6,11:50)] , use="complete.obs" , method = "pearson" )
corrplot(correlation_matrix_total, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, tl.cex = 0.7)
##this is to add pval information##
corrplot(correlation_matrix_2$r, type="upper", order="hclust",
         p.mat = correlation_matrix_2$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.srt = 45, tl.cex = 0.7)

###########PCA########

library("FactoMineR")
library("factoextra")

#all dist pheno fc

#Andrew(84) drives data dispersion too much
all_dist_pheno_fc_no84<-slice(all_dist_pheno_fc, c(1:83,85))
scaled_all_dist_pheno_fc<-scale(all_dist_pheno_fc_no84[,2:352])
scaled_all_dist_pheno_fc<-cbind(all_dist_pheno_fc_no84[c(1,353:358)],scaled_all_dist_pheno_fc)
res.pca <- prcomp(scaled_all_dist_pheno_fc[8:358])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","none"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "PCA of FC of Feature Distance") + 
  theme(legend.position = "none") +
  geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#sex
groups <- as.factor(scaled_all_dist_pheno_fc$Sex) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "PCA of FC of Feature Distance") + 
  theme(legend.position = "none") +   
  geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,
                            24,24,24,24,24,25,25,25,25,25,22),
             col.ind = scaled_all_dist_pheno_fc$Age, 
             fill.ind = scaled_all_dist_pheno_fc$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of FC of Feature Distance") +
  scale_fill_gradientn(colours = heat.colors(5))+
  #theme_classic()+
  theme(legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#PCA fc with Andrew(84)
scaled_all_dist_pheno_fc<-scale(all_dist_pheno_fc[,2:352])
scaled_all_dist_pheno_fc<-cbind(all_dist_pheno_fc[c(1,353:358)],scaled_all_dist_pheno_fc)
res.pca <- prcomp(scaled_all_dist_pheno_fc[8:358])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c"), 
             alpha.ind = 0.4 , 
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "PCA of FC of Feature Distance") + 
  theme(legend.position = "none") +
  geom_text(aes(label = scaled_all_dist_pheno_fc$Co_twins_ID)) 

#sex
groups <- as.factor(all_dist_pheno_fc$Sex) # group by Sex
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,
                            24,24,24,24,24,25,25,25,25,25),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "PCA of FC of Feature Distance") + theme(legend.position = "none") +   geom_text(aes(label = Co_twins_ID)) 

#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_dist_pheno_fc$Age.x, 
             fill.ind = all_dist_pheno_fc$Age.x,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "PCA of FC of Feature Distance") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#all dist pheno delta

all_dist_pheno_delta_no84<-slice(all_dist_pheno_delta, c(1:83,85))
scaled_all_dist_pheno_delta<-scale(all_dist_pheno_delta_no84[2:352])
scaled_all_dist_pheno_delta<-cbind(all_dist_pheno_delta_no84[c(1,353:358)],scaled_all_dist_pheno_delta)
res.pca <- prcomp(scaled_all_dist_pheno_delta[8:358])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Feature Distances x Delta between Cotwins") + theme(legend.position = "none") +
  geom_text(aes(label = scaled_all_dist_pheno_delta$Co_twins_ID)) 

#sex
groups <- as.factor(scaled_all_dist_pheno_delta$Sex.x)
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x Delta between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x Delta between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#all dist pheno absdelta
scaled_all_dist_pheno_absdelta<-scale(all_dist_pheno_absdelta[2:352])
scaled_all_dist_pheno_absdelta<-cbind(all_dist_pheno_absdelta[c(1,353:358)],scaled_all_dist_pheno_absdelta)
res.pca <- prcomp(scaled_all_dist_pheno_absdelta[8:358])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x absDelta between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(scaled_all_dist_pheno_absdelta$Sex.x) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x absDelta between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = scaled_all_dist_pheno_absdelta$Age.x, 
             fill.ind = scaled_all_dist_pheno_absdelta$Age.x,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x absDelta between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "ind", axes = 1:2, top = 120)  +
  theme(axis.text.x = element_text(angle=90))
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 20)

#all ratio pheno delta
res.pca <- prcomp(all_ratio_pheno_delta[2:61426])
#controls
groups <- as.factor(all_ratio_pheno_delta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x Delta between Cotwins") + 
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
                            24,24,24,24,24,25,25,25,25,25,22,22),
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
                            24,24,24,24,24,25,25,25,25,25,22,22),
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


#all ratio pheno absdelta
res.pca <- prcomp(all_ratio_pheno_absdelta[2:61426])
#controls
groups <- as.factor(all_ratio_pheno_absdelta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x Delta between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_absdelta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
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
                            24,24,24,24,24,25,25,25,25,25,22,22),
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

#all ratio pheno fc
res.pca <- prcomp(all_ratio_pheno_fc[2:61426])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_fc$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x FC between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(pca, choice = "var", axes = 1:2, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#all ratio pheno log2fc
res.pca <- prcomp(all_ratio_pheno_log2_fc[2:61426])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_fc$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x FC between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

#all ratio pheno abslog2fc
res.pca <- prcomp(all_ratio_pheno_abslog2_fc[2:61426])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_abslog2_fc$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x FC between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

#all ratio pheno log2delta
res.pca <- prcomp(all_ratio_pheno_log2_delta[2:61426])
#controls
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_log2_delta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x FC between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#all ratio pheno log2absdelta
res.pca <- prcomp(all_ratio_pheno_log2_absdelta[2:61426])
#controls
groups <- as.factor(all_ratio_pheno_log2_absdelta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_log2_delta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x FC between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(pca, choice = "var", axes = 1:2, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

#all ratio pheno abslog2 absdelta
res.pca <- prcomp(all_ratio_pheno_abslog2_absdelta[2:61426])
#controls
groups <- as.factor(all_ratio_pheno_abslog2_absdelta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = "black", fill.ind = c("none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "none","none","none","none","none","none","none","none","none","none","none","none","none",
                                             "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#668d3c","#668d3c","#668d3c","#668d3c","#668d3c","none","none"), 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#sex
groups <- as.factor(all_ratio_pheno_log2_absdelta$Sex) # group by Twin_ID
fviz_pca_ind(res.pca, pointsize = 10, axes = c(1, 2), repel = T, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = groups, fill.ind = groups, 
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             palette = c(), 
             #legend.title = "Sex" ) #+ theme(legend.position = "none") 
             title = "Ratio of Distances x FC between Cotwins") + theme(legend.position = "none")
#age
fviz_pca_ind(res.pca, pointsize = 10, 
             axes = c(1, 2), repel = F, 
             geom.ind = "point", 
             mean.point = F, 
             pointshape = c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                            24,24,24,24,24,25,25,25,25,25,22,22),
             col.ind = all_ratio_pheno_delta$Age, 
             fill.ind = all_ratio_pheno_delta$Age,
             #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             alpha.ind = 0.4 , 
             #habillage = groups,# color by groups
             #palette = c(), 
             #legend.title = "Age" ) #+ theme(legend.position = "none")
             title = "Ratio of Distances x FC between Cotwins") +
  scale_fill_gradientn(colours = heat.colors(5))+
  theme_classic()+
  theme(legend.position = "none")

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 30)  +
  theme(axis.text.x = element_text(angle=90))

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 30)  +
  theme(axis.text.x = element_text(angle=90))
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 30)+
  theme(axis.text.x = element_text(angle=90))


##### t-SNE#######
library(tidyverse)

library(Rtsne)

library(ggplot2)

#all dist pheno fc

scaled_all_dist_pheno_fc<-scale(all_dist_pheno_fc[2:352])
scaled_all_dist_pheno_fc<-cbind(all_dist_pheno_fc[c(1,353:358)],scaled_all_dist_pheno_fc)

tsne <- Rtsne(scaled_all_dist_pheno_fc[,8:358],check_duplicates = FALSE, perplexity=5, max_iter = 20000, theta = 0.0) 

#plot t-sne results (you put it in a plot)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_fc$Co_twins_ID, 
                                   fill=scaled_all_dist_pheno_fc$Co_twins_ID), 
                               alpha = 0.3, size = 10, shape=21) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_fc$Co_twins_ID)) +
  labs(title = "t-SNE of distances in pop", x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

#all dist pheno delta: both sexes x 2, just male, just female

scaled_all_dist_pheno_delta<-scale(all_dist_pheno_delta[2:352])
scaled_all_dist_pheno_delta<-cbind(all_dist_pheno_delta[c(1,353:358)],scaled_all_dist_pheno_delta)

tsne <- Rtsne(scaled_all_dist_pheno_delta[,8:358],check_duplicates = FALSE, perplexity=5, max_iter = 20000, theta = 0.0) 

#plot t-sne results (you put it in a plot)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2]) # 1 and 2 are the dimensions of t-sne

#color and modify t-sne plot
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_delta$Sex.x, 
                                   fill=scaled_all_dist_pheno_delta$Sex.x), 
                               alpha = 0.3, size = 10, shape=21) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_delta$Co_twins_ID)) +
  labs(title = "t-SNE of distances in pop", x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_delta$Sex.x, 
                                   fill=scaled_all_dist_pheno_delta$Sex.x), 
                               alpha = 0.3, size = 10, shape=21) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_delta$Co_twins_ID)) +
  labs(title = "t-SNE of distances in pop", x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=scaled_all_dist_pheno_delta$Sex.x, 
                                   fill=scaled_all_dist_pheno_delta$Sex.x), 
                               alpha = 0.3, size = 10, shape=21) + 
  theme_classic() + # alpha is the opacity
  geom_text(aes(x=x, y=y, label=scaled_all_dist_pheno_delta$Co_twins_ID)) +
  labs(title = "t-SNE of distances in pop", x="tSNE dim1", y="tSNE dim2", color="Co_twins_ID", fill="Co_twins_ID") + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 


#import of tsnes from server and plots of ratio matrices 

#tab_all_ratio_pheno_absdelta<-read.table(file="/Volumes/projects_secondary/mnp/warren/test/tsne_all_ratio_delta.tab", sep ="\t", header=F)
#tab_all_ratio_pheno_fc<-read.table(file="/Volumes/projects_secondary/mnp/warren/test/tsne_all_ratio_fc.tab", sep ="\t", header=F)
tab_all_ratio_pheno_delta_nops<-read.table(file="/Volumes/projects_secondary/mnp/warren/test/tsne_all_ratio_delta_nops.tab",sep="\t",header=F)
tab_all_ratio_pheno_fc_nops<-read.table(file="/Volumes/projects_secondary/mnp/warren/test/tsne_all_ratio_fc_nops.tab",sep="\t",header=F)

#all_ratio_pheno_delta_tsne<-cbind(all_dist_pheno_fc[c(1,406:411)],tab_all_ratio_pheno_absdelta)
#all_ratio_pheno_fc_tsne<-cbind(all_dist_pheno_fc[c(1,407:413)],tab_all_ratio_pheno_fc)
all_ratio_pheno_delta_nops_tsne<-cbind(all_dist_pheno_fc[c(1,353:358)],tab_all_ratio_pheno_delta_nops)
all_ratio_pheno_fc_nops_tsne<-cbind(all_dist_pheno_fc[1:83,c(1,353:358)],tab_all_ratio_pheno_fc_nops)

all_ratio_pheno_fc_nops_tsne_male<-all_ratio_pheno_fc_nops_tsne%>%filter(.,Sex=="M")
all_ratio_pheno_fc_nops_tsne_female<-all_ratio_pheno_fc_nops_tsne%>%filter(.,Sex=="F")
all_ratio_pheno_delta_nops_tsne_male<-all_ratio_pheno_delta_nops_tsne%>%filter(.,Sex.x=='M')
all_ratio_pheno_delta_nops_tsne_female<-all_ratio_pheno_delta_nops_tsne%>%filter(.,Sex.x=='F')

#ratio matrix tsne unranked/uncorrelated
ggplot(all_ratio_pheno_delta_nops_tsne[,8:9]) + 
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

#delta ratio matrix tsne: only males, colored by age
ggplot(all_ratio_pheno_delta_nops_tsne_male[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_delta_nops_tsne_male[,8], y=all_ratio_pheno_delta_nops_tsne_male[,9], color=all_ratio_pheno_delta_nops_tsne_male[,8], fill=all_ratio_pheno_delta_nops_tsne_male[,7]), alpha = 0.3, size = 10, shape=21) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_delta_nops_tsne_male[,8], y=all_ratio_pheno_delta_nops_tsne_male[,9], label=all_ratio_pheno_delta_nops_tsne_male[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_delta_nops_tsne_male[,7], fill=all_ratio_pheno_delta_nops_tsne_male[,7]) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

#delta ratio matrix tsne: only females, colored by age
ggplot(all_ratio_pheno_delta_nops_tsne_female[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_delta_nops_tsne_female[,8], y=all_ratio_pheno_delta_nops_tsne_female[,9], color=all_ratio_pheno_delta_nops_tsne_female[,8], fill=all_ratio_pheno_delta_nops_tsne_female[,7]), alpha = 0.3, size = 10, shape=21) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_delta_nops_tsne_female[,8], y=all_ratio_pheno_delta_nops_tsne_female[,9], label=all_ratio_pheno_delta_nops_tsne_female[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_delta_nops_tsne_female[,7], fill=all_ratio_pheno_delta_nops_tsne_female[,7]) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

#all ratio pheno fc coloring just the controls
ggplot(all_ratio_pheno_fc_nops_tsne[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_fc_nops_tsne[,8], 
                 y=all_ratio_pheno_fc_nops_tsne[,9],
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
             alpha = 0.8, size = 20, shape= c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              21,21,21,21,21,21,21,21,21,21,21,21,21,
                                              24,24,24,24,24,
                                              22,22,22,22,22)) +  
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_fc_nops_tsne[,8], y=all_ratio_pheno_fc_nops_tsne[,9], label=all_ratio_pheno_fc_nops_tsne[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_fc_nops_tsne[,7], fill=all_ratio_pheno_fc_nops_tsne[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22),
        legend.position = "none",
        legend.title=element_blank(), 
        legend.text=element_text(size=20))

#fc ratio matrix tsne: both sexes, colored by age
ggplot(all_ratio_pheno_fc_nops_tsne[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_fc_nops_tsne[,8], y=all_ratio_pheno_fc_nops_tsne[,9], color=all_ratio_pheno_fc_nops_tsne[,7], fill=all_ratio_pheno_fc_nops_tsne[,7]), alpha = 0.3, size = 10, shape=21) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_fc_nops_tsne[,8], y=all_ratio_pheno_fc_nops_tsne[,9], label=all_ratio_pheno_fc_nops_tsne[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_fc_nops_tsne[,7], fill=all_ratio_pheno_fc_nops_tsne[,7]) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

#all ratio pheno fc by sex
ggplot(all_ratio_pheno_fc_nops_tsne[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_fc_nops_tsne[,8], y=all_ratio_pheno_fc_nops_tsne[,9], fill=all_ratio_pheno_fc_nops_tsne[,2]), alpha = 0.3, size = 10, shape=21) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_fc_nops_tsne[,8], y=all_ratio_pheno_fc_nops_tsne[,9], label=all_ratio_pheno_fc_nops_tsne[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_fc_nops_tsne[,7], fill=all_ratio_pheno_fc_nops_tsne[,7]) + # labs = labels
  #scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))

#fc ratio matrix tsne: only males, colored by age
ggplot(all_ratio_pheno_fc_nops_tsne_male[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_fc_nops_tsne_male[,8], y=all_ratio_pheno_fc_nops_tsne_male[,9], color=all_ratio_pheno_fc_nops_tsne_male[,7], fill=all_ratio_pheno_fc_nops_tsne_male[,7]), alpha = 0.3, size = 10, shape=21) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_fc_nops_tsne_male[,8], y=all_ratio_pheno_fc_nops_tsne_male[,9], label=all_ratio_pheno_fc_nops_tsne_male[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_fc_nops_tsne_male[,2], fill=all_ratio_pheno_fc_nops_tsne_male[,7]) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

ggplot(all_ratio_pheno_fc_nops_tsne_female[,8:9]) + 
  geom_point(aes(x=all_ratio_pheno_fc_nops_tsne_female[,8], y=all_ratio_pheno_fc_nops_tsne_female[,9], color=all_ratio_pheno_fc_nops_tsne_female[,7], fill=all_ratio_pheno_fc_nops_tsne_female[,7]), alpha = 0.3, size = 10, shape=21) + 
  theme_classic() +# alpha is the opacity 
  geom_text(aes(x=all_ratio_pheno_fc_nops_tsne_female[,8], y=all_ratio_pheno_fc_nops_tsne_female[,9], label=all_ratio_pheno_fc_nops_tsne_female[,1])) +
  labs(title = "t-SNE of ratio distances", x="tSNE dim1", y="tSNE dim2", color=all_ratio_pheno_fc_nops_tsne_female[,2], fill=all_ratio_pheno_fc_nops_tsne_female[,7]) + # labs = labels
  scale_colour_gradientn(aesthetics = c("colour", "fill"), colors =rev(heat.colors(10))) + # this is to use a color scale rather then discrete colors ## heat.colors = palette
  theme(plot.title = element_text(hjust = 0.5, size = 25), ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_blank(), 
        legend.text=element_text(size=20))#, 
#legend.position = "none") 

#####Density#####

# skip until next hashtag comment
obama<-read.csv(file = "Obama_cf.csv", header = T, stringsAsFactors = FALSE)
andrew<-read.csv(file = "Andrew_cf.csv", header = T, stringsAsFactors = FALSE)

andrew<-slice(andrew, 1)
obama<-slice(obama, 2)
andrew_obama<-rbind(andrew, obama)
andrew_obama[2,c(2,3)]<-c("Barack_2","Barack_1")
andrew_obama = andrew_obama %>% rename(.,status = X)
andrew_obama = andrew_obama %>% rename(.,Names = target)
andrew_obama = left_join(andrew_obama, pheno_data,by = "Names")
andrew_obama = andrew_obama %>% 
  rename(., target = Names) %>%
  rename(., Names = source)
andrew_obama = left_join(andrew_obama, pheno_data,by = "Names")
andrew_obama = andrew_obama %>% 
  rename(., source = Names)
andrew_obama$status<-c("negative_ctrl","negative_ctrl")
similarity_scores_pheno_aldi<-rbind(andrew_obama,similarity_scores_pheno_aldi)

rd_0 <- read.csv(file = 'pos_ctrl_all_cf.csv', header = TRUE, stringsAsFactors = FALSE)

status <- c("positive_ctrl","positive_ctrl","positive_ctrl",'positive_ctrl',"positive_ctrl")

rd_0 = rd_0 %>% rename(.,status = X)
rd_0 = rd_0 %>% rename(.,Names = target)
rd_0 = left_join(rd_0, cf_phenotype_data,by = "Names")
rd_0 = rd_0 %>% 
  rename(., target = Names) %>%
  rename(., Names = source)
rd_0 = left_join(rd_0, cf_phenotype_data,by = "Names")
rd_0 = rd_0 %>% 
  rename(., source = Names)

similarity_scores_3<-mutate(similarity_scores,target=as.character(target))
similarity_scores_3<-mutate(similarity_scores_3,target=sapply(strsplit(similarity_scores_3$target, split='.', fixed=TRUE),function(x) (x[1])))
similarity_scores_3<-mutate(similarity_scores_3,source=as.character(source))
similarity_scores_3<-mutate(similarity_scores_3,source=sapply(strsplit(similarity_scores_3$source, split='.', fixed=TRUE),function(x) (x[1])))
similarity_scores_3$target <- ifelse(similarity_scores_3$target %in% c('Prijatel', 'Martin_emily'), c('Prijatel_Karen',"Martin_Emily"), similarity_scores_3$target)
similarity_scores_3$source <- ifelse(similarity_scores_3$source %in% c('Prijatel', 'Martin_emily'), c('Prijatel_Karen',"Martin_Emily"), similarity_scores_3$source)

similarity_scores_3 = similarity_scores_3[,-5]

similarity_scores_3 = similarity_scores_3 %>% rename(.,Names = target)

similarity_scores_pheno = left_join(similarity_scores_3, cf_phenotype_data,by = "Names")

similarity_scores_pheno = similarity_scores_pheno %>% 
  rename(., target = Names) %>%
  rename(., Names = source)

similarity_scores_pheno = left_join(similarity_scores_pheno, cf_phenotype_data,by = "Names")

similarity_scores_pheno = similarity_scores_pheno %>% 
  rename(., source = Names)

similarity_scores_pheno<-similarity_scores_pheno[!(similarity_scores_pheno$target=="Parks_Katie copy"),]
similarity_scores_pheno = left_join(similarity_scores_pheno, phenotype_data,by = "Names")
similarity_scores_pheno<-distinct(similarity_scores_pheno)

similarity_scores_pheno<-similarity_scores_pheno[-c(6482:6485),]
similarity_scores_pheno_aldi <- bind_rows(similarity_scores_pheno,rd_0)
library("plyr")
similarity_scores_pheno_aldi <- transform(similarity_scores_pheno_aldi,
                                          status=revalue(status,c("positive_ctrl"="pos_ctrl")))

#begin here for density plots

similarity_scores_pheno_aldi_sex<-similarity_scores_pheno_aldi[!(similarity_scores_pheno_aldi$Sex.x=="M" & similarity_scores_pheno_aldi$Sex.y=="F"),]
similarity_scores_pheno_aldi_sex<-similarity_scores_pheno_aldi_sex[!(similarity_scores_pheno_aldi_sex$Sex.x=="F" & similarity_scores_pheno_aldi_sex$Sex.y=="M"),]
similarity_scores_pheno_aldi_sex_female<-similarity_scores_pheno_aldi_sex[!(similarity_scores_pheno_aldi_sex$Sex.x=="M"),]
similarity_scores_pheno_aldi_sex_male<-similarity_scores_pheno_aldi_sex[!(similarity_scores_pheno_aldi_sex$Sex.x=="F"),]

similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex[!(similarity_scores_pheno_aldi_sex$Ethnicity.x=="caucasian" & similarity_scores_pheno_aldi_sex$Ethnicity.y=="hispanic"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="hispanic" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="caucasian"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="caucasian" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="black"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="black" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="caucasian"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="caucasian" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="asian"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="caucasian"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="black" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="asian"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="black" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="hispanic"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="black"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="hispanic" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="black"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="hispanic" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="asian"),]
similarity_scores_pheno_aldi_sex_ethnicity<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.y=="hispanic"),]

similarity_scores_pheno_aldi_sex_ethnicity_black<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="caucasian" | similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="hispanic" | similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="asian"),]
similarity_scores_pheno_aldi_sex_ethnicity_white<-similarity_scores_pheno_aldi_sex_ethnicity[!(similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="black" | similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="hispanic" | similarity_scores_pheno_aldi_sex_ethnicity$Ethnicity.x=="asian"),]


similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi[!(similarity_scores_pheno_aldi$Ethnicity.x=="caucasian" & similarity_scores_pheno_aldi$Ethnicity.y=="hispanic"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="hispanic" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="caucasian"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="caucasian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="black"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="black" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="caucasian"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="caucasian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="asian"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="caucasian"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="black" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="asian"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="black" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="hispanic"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="black"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="hispanic" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="black"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="hispanic" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="asian"),]
similarity_scores_pheno_aldi_ethnicity<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="hispanic"),]

similarity_scores_pheno_aldi_ethnicity_black<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="hispanic"),]
similarity_scores_pheno_aldi_ethnicity_white<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Ethnicity.x=="asian" & similarity_scores_pheno_aldi_ethnicity$Ethnicity.y=="hispanic"),]

similarity_scores_pheno_aldi_ethnicity_white_sex<-similarity_scores_pheno_aldi_ethnicity_white[!(similarity_scores_pheno_aldi_ethnicity_white$Sex.x=="M" & similarity_scores_pheno_aldi_ethnicity_white$Sex.y=="F"),]
similarity_scores_pheno_aldi_ethnicity_white_sex<-similarity_scores_pheno_aldi_ethnicity_white_sex[!(similarity_scores_pheno_aldi_ethnicity_white_sex$Sex.x=="F" & similarity_scores_pheno_aldi_ethnicity_white_sex$Sex.y=="M"),]
similarity_scores_pheno_aldi_ethnicity_white_sex_male<-similarity_scores_pheno_aldi_ethnicity_white_sex[!(similarity_scores_pheno_aldi_ethnicity_white_sex$Sex.x=="F"),]
similarity_scores_pheno_aldi_ethnicity_white_sex_female<-similarity_scores_pheno_aldi_ethnicity_white_sex[!(similarity_scores_pheno_aldi_ethnicity_white_sex$Sex.x=="M"),]

similarity_scores_pheno_aldi_ethnicity_black_sex<-similarity_scores_pheno_aldi_ethnicity_black[!(similarity_scores_pheno_aldi_ethnicity_white$Sex.x=="M" & similarity_scores_pheno_aldi_ethnicity_white$Sex.y=="F"),]
similarity_scores_pheno_aldi_ethnicity_black_sex<-similarity_scores_pheno_aldi_ethnicity_black_sex[!(similarity_scores_pheno_aldi_ethnicity_white$Sex.x=="F" & similarity_scores_pheno_aldi_ethnicity_white$Sex.y=="M"),]
similarity_scores_pheno_aldi_ethnicity_black_sex_male<-similarity_scores_pheno_aldi_ethnicity_black_sex[!(similarity_scores_pheno_aldi_ethnicity_white$Sex.x=="F"),]
similarity_scores_pheno_aldi_ethnicity_black_sex_female<-similarity_scores_pheno_aldi_ethnicity_black_sex[!(similarity_scores_pheno_aldi_ethnicity_white$Sex.x=="M"),]

similarity_scores_pheno_aldi_ethnicity_sex<-similarity_scores_pheno_aldi_ethnicity[!(similarity_scores_pheno_aldi_ethnicity$Sex.x=="M" & similarity_scores_pheno_aldi_ethnicity$Sex.y=="F"),]
similarity_scores_pheno_aldi_ethnicity_sex<-similarity_scores_pheno_aldi_ethnicity_sex[!(similarity_scores_pheno_aldi_ethnicity_sex$Sex.x=="F" & similarity_scores_pheno_aldi_ethnicity_sex$Sex.y=="M"),]
similarity_scores_pheno_aldi_ethnicity_sex_male<-similarity_scores_pheno_aldi_ethnicity_sex[!(similarity_scores_pheno_aldi_ethnicity_sex$Sex.x=="F"),]
similarity_scores_pheno_aldi_ethnicity_sex_female<-similarity_scores_pheno_aldi_ethnicity_sex[!(similarity_scores_pheno_aldi_ethnicity_sex$Sex.x=="M"),]

similarity_scores_pheno_aldi_ethnicity<-na.omit(similarity_scores_pheno_aldi_ethnicity)
similarity_scores_pheno_aldi_ethnicity_black<-na.omit(similarity_scores_pheno_aldi_ethnicity_black)
similarity_scores_pheno_aldi_ethnicity_white<-na.omit(similarity_scores_pheno_aldi_ethnicity_white)
similarity_scores_pheno_aldi_ethnicity_sex<-na.omit(similarity_scores_pheno_aldi_ethnicity_sex)
similarity_scores_pheno_aldi_sex_ethnicity<-na.omit(similarity_scores_pheno_aldi_sex_ethnicity)

library(facetscales)
library(ggplot2)
library(cowplot)
library(ggpubr)

theme <- theme(plot.title = element_text(hjust = 1, size = 20), 
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),  
               axis.text.x= element_text(size = 10), 
               axis.text.y= element_blank(), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               legend.position = "none")

scales_y <- list(
  `co-twins` = scale_y_continuous(limits = c(0, 0.7), breaks = NULL),
  `negative_ctrl` = scale_y_continuous(limits = c(0, 0.1), breaks = NULL),
  `positive_ctrl` = scale_y_continuous(limits = c(0, 1500), breaks = NULL)
)



levels(similarity_scores_pheno$Compared)
levels(similarity_scores_pheno$Compared) <- c("co-twins" ,   "negative_ctrl" ,  "pos_ctrl" ,  "neg_ctrl_*", "positive_ctrl")
similarity_scores_pheno$Compared = factor(similarity_scores_pheno$status, levels=c('negative_ctrl','co-twins','positive_ctrl','neg_ctrl_*', 'pos_ctrl'))
levels(similarity_scores_pheno$Compared)


ggplot(subset(similarity_scores_pheno_aldi, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
#y = "Density"


ggplot(subset(similarity_scores_pheno_aldi_sex, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggplot(subset(similarity_scores_pheno_aldi_sex_male, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggplot(subset(similarity_scores_pheno_aldi_sex_female, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggplot(subset(similarity_scores_pheno_aldi_ethnicity, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggplot(subset(similarity_scores_pheno_aldi_ethnicity_sex, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggplot(subset(similarity_scores_pheno_aldi_ethnicity_white, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

ggplot(subset(similarity_scores_pheno_aldi_ethnicity_black, status == 'negative_ctrl' |status == 'co-twins'| status == 'pos_ctrl'), 
       aes(x=Similarity, color=status)) +
  stat_density(aes(x=Similarity, y=..scaled..,color=status), position=position_dodge(width = .1), geom="line") +
  scale_color_manual(values = c( "#615292", "#EA5454","#AEC960")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  scale_fill_discrete(name= "Similarity Scores",labels = c("cowtins",'negative control','positive control'))+
  labs(title = "AWS-Rekognition controls", y = "Density",color = "Similarity Scores")+
  expand_limits(x=102)+
  theme_classic()+
  theme(#legend.position=c(0.1,-0.1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 30), 
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),  
    axis.text.x= element_text(size = 15, angle = 0), 
    axis.text.y= element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())

#######UMAP#######

library(umap)

scaled_adpd<-scale(all_dist_pheno_delta)
umap_scaled_adpd<-umap(scaled_adpd)
View(umap_scaled_adpd)

library(M3C)

pheno_data_paired = pheno_data[,c(3:9)]
Co_twins_ID <- c('1','1','1','2','3','4','5','6','7','8','8','8','8','8','8','9','10','11','12','13','14','15','16','17','18','18','18','18',
                 '18','18','19','20','21','22','23','23','23','24','24','24','24','24',
                 '24','25','26','27','28','29','30','31','32','33','34','35','36','37','38',
                 '39','40','41','42','43','44','45','46','47','48','49','50','51','52','53',
                 '54',
                 'pos1','pos2','pos3','pos4','pos5','neg1','neg2','neg3','neg4','neg5','55','56')

rtc_df<-rownames_to_column(all_dist_pheno_fc)
rowname<-rtc_df[,1]

all_dist_pheno_fc <- data.frame(all_dist_pheno_fc)
all_dist_pheno_fc <- cbind(Co_twins_ID, all_dist_pheno_fc)
all_dist_pheno_fc <- left_join(all_dist_pheno_fc, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_fc <- distinct(all_dist_pheno_fc,ChinBottom_LeftPupil,.keep_all= TRUE)
all_dist_pheno_fc <- cbind(pairwise_rownames,all_dist_pheno_fc)
all_dist_pheno_fc <- column_to_rownames(all_dist_pheno_fc,var="pairwise_rownames")

all_dist_pheno_delta<-data.frame(all_dist_pheno_delta)
all_dist_pheno_delta<- cbind(Co_twins_ID, all_dist_pheno_delta)
all_dist_pheno_delta<-left_join(all_dist_pheno_delta, pheno_data_paired,by = "Co_twins_ID")
all_dist_pheno_delta <- distinct(all_dist_pheno_delta,ChinBottom_LeftPupil,.keep_all= TRUE)
all_dist_pheno_delta <- cbind(pairwise_rownames,all_dist_pheno_delta)
all_dist_pheno_delta <- column_to_rownames(all_dist_pheno_delta,var="pairwise_rownames")

all_ratio_pheno_delta<-data.frame(all_ratio_pheno_delta)
all_ratio_pheno_delta<- cbind(Co_twins_ID, all_ratio_pheno_delta)
all_ratio_pheno_delta<-left_join(all_ratio_pheno_delta, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_delta <- distinct(all_ratio_pheno_delta,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)
all_ratio_pheno_delta <- cbind(pairwise_rownames,all_ratio_pheno_delta)
all_ratio_pheno_delta <- column_to_rownames(all_ratio_pheno_delta,var="pairwise_rownames")

all_ratio_pheno_fc<-data.frame(all_ratio_pheno_fc)
all_ratio_pheno_fc<- cbind(Co_twins_ID, all_ratio_pheno_fc)
all_ratio_pheno_fc<-left_join(all_ratio_pheno_fc, pheno_data_paired,by = "Co_twins_ID")
all_ratio_pheno_fc <- distinct(all_ratio_pheno_fc,ChinBottom_LeftPupil__ChinBottom_RightPupil,.keep_all= TRUE)
all_ratio_pheno_fc <- cbind(pairwise_rownames,all_ratio_pheno_fc)
all_ratio_pheno_fc <- column_to_rownames(all_ratio_pheno_fc,var="pairwise_rownames")

scaled_adpfc<-scale(all_dist_pheno_fc[,2:352])
scaled_adpfc<-cbind(all_dist_pheno_fc[,c(1,353:358)],scaled_adpfc)
umap(scaled_adpfc[,8:358],colvec=c('skyblue'))
umap(t(scaled_adpfc[,8:358]),labels=as.factor(scaled_adpfc[,2]), text = as.character(scaled_adpfc[,1]),)

scaled_adpd<-scale(all_dist_pheno_delta[,2:352])
scaled_adpd<-cbind(all_dist_pheno_delta[,c(1,353:358)],scaled_adpd)
umap(t(scaled_adpd[,8:358]),labels=as.factor(scaled_adpd[,2]), text = as.character(scaled_adpd[,1]))

umap(t(all_ratio_pheno_delta[,2:61426]),labels=as.factor(all_ratio_pheno_delta[,61427]), text = as.character(all_ratio_pheno_delta[,1]))
umap(t(all_ratio_pheno_fc[,2:61426]),labels=as.factor(all_ratio_pheno_fc[,61427]), text = as.character(all_ratio_pheno_delta[,1]))

scaled_arpd<-scale(all_ratio_pheno_delta[,2:61426])
scaled_arpd<-cbind(all_dist_pheno_delta[,c(1,353:358)],scaled_arpd)

scaled_arpfc<-scale(all_ratio_pheno_fc[,2:61426])
scaled_arpfc<-cbind(all_dist_pheno_delta[,c(1,353:358)],scaled_arpfc)

umap(t(scaled_arpd[,8:61432]),labels=as.factor(scaled_arpd[,2]), text = as.character(scaled_arpd[,1]))
umap(t(scaled_arpfc[,8:61432]),labels=as.factor(scaled_arpfc[,2]), text = as.character(scaled_arpfc[,1]))


library(uwot)

t <- theme(plot.title = element_text(hjust = 0.5, size = 25), 
           axis.title.x = element_text(size = 20),
           axis.title.y = element_text(size = 20),  
           axis.text.x= element_text(size = 20), 
           axis.text.y= element_text(size = 20),
           legend.position = c(0.15,0.80),
           legend.title = element_text(size = 15), 
           legend.text = element_text(size = 15))

umap_arpd <- umap(all_ratio_pheno_delta[,c(2:61426)])
umap_arpd<-as.data.frame(umap_arpd)
pdf("UMAP_arpd.pdf", w=10, h=8)
ggplot(umap_arpd) + 
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
  geom_text(aes(x=V1, y=V2, label=all_ratio_pheno_delta$Co_twins_ID), size = 3) +
  labs(title = "UMAP - individuals", x="UMAP x", y="UMAP y") + theme_classic() + t+
  scale_color_manual(values=c( "#EA5454", "#AEC960", "#615292"), guide=FALSE) +
  scale_fill_manual(values=c( "#EA5454", "#AEC960", "white"), guide=FALSE)
dev.off()

pdf("UMAP_arpd_noName.pdf", w=10, h=8)
ggplot(umap_arpd) + 
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

pdf("UMAP_arpd_noName_TWINS.pdf", w=10, h=8)
ggplot(umap_arpd) + 
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









