setwd("/Users/Warren.Sink/github/facial-polyphenism/CSVs")

library(tidyverse)

df <- read.csv(file = 'msu_coord_nonnosenorm.csv', header = TRUE, stringsAsFactors = FALSE)

colnames(df)<-df[1,]
df2<-df[-1,]

rownames(df2)<-df2[,1]
df2<-df2[,-1]
msu_coord_3n<-df2

for (i in (1:ncol(msu_coord_3n))){
  for (j in (1:nrow(msu_coord_3n))){
    msu_coord_3n[j,i]<-str_remove(string = msu_coord_3n[j,i], pattern = "\\[|\\]")
    msu_coord_3n[j,i]<-str_remove(string = msu_coord_3n[j,i], pattern = "\\]|\\[")
    
  }
}

t_msu_coord_3n <- as.data.frame(t(msu_coord_3n))
rnames_t_msu_coord_3n <- colnames(msu_coord_3n)

t_msu_coord_3n <- splitstackshape::cSplit(t_msu_coord_3n, names(t_msu_coord_3n)) 
t_msu_coord_3n <- as.data.frame(t_msu_coord_3n)
rownames(t_msu_coord_3n) <- rnames_t_msu_coord_3n

msu_coord_3n <- t_msu_coord_3n %>% t() 

testdf <- data.frame()
colnam_msu_coord_3n <- colnames(msu_coord_3n) 
rownam_msu_coord_3n <- rownames(msu_coord_3n)
for (i in (1:ncol(msu_coord_3n))){
  for (j in (1:nrow(msu_coord_3n))){
    if (j %% 2 > 0){
      testdf[j,i]<-jpeg_dim[i,2]*(1-msu_coord_3n[j,i])
    } else {
      testdf[j,i]<-jpeg_dim[i,3]*(1-msu_coord_3n[j,i])
    }
  }
}
colnames(testdf) <- colnam_msu_coord_3n
rownames(testdf) <- rownam_msu_coord_3n

msu_coord_y <- seq(0, nrow(msu_coord_3n), 2)
msu_coord_x <- seq(1, nrow(msu_coord_3n), 2)

msu_coord_3n <- testdf
msu_coord_3n_y <- msu_coord_3n[msu_coord_y,]
msu_coord_3n_x <- msu_coord_3n[msu_coord_x,]

msu_coord_3n <- as.data.frame(msu_coord_3n)
msu_coord_3n_y <- as.data.frame(msu_coord_3n_y)
msu_coord_3n_x <- as.data.frame(msu_coord_3n_x)

write.csv(msu_coord_3n, file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/msu_extract_nonnosenorm.csv")
write.csv(msu_coord_3n_y, file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/msu_extract_nonnosenorm_y.csv")
write.csv(msu_coord_3n_x, file = "/Users/Warren.Sink/github/facial-polyphenism/CSVs/msu_extract_nonnosenorm_x.csv")








