setwd("/Users/Warren.Sink/Desktop/R projects/facial-polyphenism/CSVs")

library(tidyverse)

df <- read.csv(file = 'schoeller_coord_non_nose_norm.csv', header = TRUE, stringsAsFactors = FALSE)

colnames(df)<-df[1,]
df2<-df[-1,]

#df2[46,1]<-"RightEyebrowLeft_y"
rownames(df2)<-df2[,1]
df2<-df2[,-1]
schoeller_coord_3n<-df2

#Dimensions = 1275 * 2100
dim(schoeller_coord_3n)

rownames(schoeller_coord_3n)[[2]]

for (i in (1:ncol(schoeller_coord_3n))){
  for (j in (1:nrow(schoeller_coord_3n))){
    schoeller_coord_3n[j,i]<-str_remove(string = schoeller_coord_3n[j,i], pattern = "\\[|\\]")
    schoeller_coord_3n[j,i]<-str_remove(string = schoeller_coord_3n[j,i], pattern = "\\]|\\[")
  }
}


#needs work
rnames_df2 <- row.names(df2)
schoeller_coord_3n <- lapply(schoeller_coord_3n, function(x) as.numeric(as.character(x)))
schoeller_coord_3n_df <- as.data.frame(schoeller_coord_3n)
rownames(schoeller_coord_3n) <- rnames_df2

#for (i in (1:ncol(schoeller_coord_3n_df))){
#  for (j in (1:nrow(schoeller_coord_3n_df))){
#    if (j %% 2 > 0){
#      schoeller_coord_3n_df[j,i]<-(1-schoeller_coord_3n_df[j,i])
#    } else {
#      schoeller_coord_3n_df[j,i]<-(1-schoeller_coord_3n_df[j,i])
#    }
#  }
#}

for (i in (1:ncol(schoeller_coord_3n_df))){
  for (j in (1:nrow(schoeller_coord_3n_df))){
    if (j %% 2 > 0){
      schoeller_coord_3n_df[j,i]<-1275*(1-schoeller_coord_3n_df[j,i])
    } else {
      schoeller_coord_3n_df[j,i]<-2100*(1-schoeller_coord_3n_df[j,i])
    }
  }
}

rnames<-row.names(df2)
rownames(schoeller_coord_3n_df)<-rnames

schoeller_coord_y <- seq(0, nrow(schoeller_coord_3n_df), 2)
schoeller_coord_x <- seq(1, nrow(schoeller_coord_3n_df), 2)

schoeller_coord_3n_df_y <- schoeller_coord_3n_df[schoeller_coord_y,]
schoeller_coord_3n_df_x <- schoeller_coord_3n_df[schoeller_coord_x,]


write.csv(schoeller_coord_3n_df, file = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/schoeller_coord_3n_df.csv")
write.csv(schoeller_coord_3n_df_y, file = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/schoeller_coord_3n_df_y.csv")
write.csv(schoeller_coord_3n_df_x, file = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/schoeller_coord_3n_df_x.csv")

schoeller_coord_3n_csv <- read_csv("Schoeller_IMGs/schoeller_coord_3n_df.csv")

