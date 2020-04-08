library(magick)
library(tidyverse)

setwd("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/msu-images-extract")
img_file_list <- list.files(pattern = "\\.jpg$")

jpeg_dim <- data.frame()

for (i in 1:length(img_file_list)){
  temp_data <- image_read(img_file_list[i])
  img <- image_info(temp_data)
  jpeg_dim <- rbind(jpeg_dim, img)
}








