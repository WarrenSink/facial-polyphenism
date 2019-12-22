  setwd("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute")
  
library(tidyverse)
  
bb <- read.csv(file = 'schoeller_coord_bb.csv', header = TRUE, stringsAsFactors = FALSE)
  
bb <- bb %>%  filter(X %in% c("ID","BB_Height","BB_Left","BB_Top","BB_Width"))
  
for (i in (1:ncol(bb))){
  for (j in (1:nrow(bb))){
    bb[j,i]<-str_remove(string = bb[j,i], pattern = "\\[|\\]")
    bb[j,i]<-str_remove(string = bb[j,i], pattern = "\\]|\\[")
  }
}
  
colnames(bb)<-bb[1,]
bb<-bb[-1,]
  
rownames(bb)<-bb[,1]
bb<-bb[,-1]
  
rnames_bb <- rownames(bb)

bb <- lapply(bb, function(x) as.numeric(as.character(x)))
bb <- as.data.frame(bb)
rownames(bb) <- rnames_bb

for (i in (1:ncol(bb))){
  for (j in (1:nrow(bb))){
    if (j %% 2 > 0){
      bb[j,i]<-2100*(bb[j,i])
    } else {
      bb[j,i]<-1275*(bb[j,i])
    }
  }
}

bb_t <- t(bb)

bb_corners = as.data.frame(bb_t) %>%
  rename("BB_Bottom" = "BB_Top") %>%
  mutate(BB_Right = BB_Left + BB_Width) %>%
  mutate(BB_Top = BB_Bottom + BB_Height) %>%
  select(BB_Left, BB_Right, BB_Bottom,  BB_Top) 

bb_x <- bb_corners %>%
  select(BB_Left, BB_Right)
bb_y <- bb_corners %>%
  select(BB_Bottom, BB_Top)

bb_x <- bb_x %>%
  mutate(BB_TopLeft = BB_Left) %>%
  mutate(BB_RightxTop = BB_Right) %>% 
  mutate(BB_return = BB_Left) %>%
  .[,c(1,3,4,2,5)]

bb_y <- bb_y %>%
  mutate(BB_TopRight = BB_Top) %>%
  mutate(BB_BottomRight = BB_Bottom) %>%
  mutate(BB_return = BB_Bottom)

bb_y = t(bb_y)
bb_x = t(bb_x)

write.csv(bb_corners, file = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/schoeller_bb_corners.csv", row.names=FALSE)
write.csv(bb_y, file = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/schoeller_bb_corners_y.csv", row.names=FALSE)
write.csv(bb_x, file = "/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/Schoeller_IMGs/schoeller_bb_corners_x.csv", row.names=FALSE)




