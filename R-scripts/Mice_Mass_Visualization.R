
library(tidyverse)
library(readxl)

trim_weights <- read_excel(path = "Downloads/FVB.Trim28_2904-3258.xlsx", col_names = FALSE)

litter <- c("1","1",'1','1','1','1','1','1','1','1','1','1','1','1',
            '2','2','2','2','2','2','2','2','2',
            '3','3','3','3','3','3','3','3','3','3','3','3',
            '4','4','4','4','4','4','4','4','4',
            '5','5','5','5','5','5','5','5',
            '6','6','6','6','6','6','6','6')

trim_weights[,2] <- litter



pdf("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/FVB.Trim28_2904-3258_all.pdf", w = 12, h = 10)
trim_weights %>%
  drop_na() %>%
  ggplot(aes(x = ...10, y = ...9, color = ...2)) +
  geom_point(aes(shape = ...4), size = 5) +
  labs(x = "Lean Mass (g)", y = "Fat mass (g)", title = "FVB.Trim28(D9) Male/Female Mice from 2904-3258", shape="Genotype", colour="Litter") +
  scale_color_discrete(labels = c("2904-2917", "3033-3041","3042-3052","3128-3136","3140-3147","3251-3258")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        #legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/FVB.Trim28_2904-3258_males.pdf", w = 12, h = 10)
trim_weights %>%
  filter(...5 == "m") %>%
  drop_na() %>%
  ggplot(aes(x = ...10, y = ...9, color = ...2)) +
  geom_point(aes(shape = ...4), size = 5) +
  labs(x = "Lean Mass (g)", y = "Fat mass (g)", title = "FVB.Trim28(D9) Male Mice from 2904-3258", shape="Genotype", colour="Litter") +
  scale_color_discrete(labels = c("2904-2917", "3033-3041","3042-3052","3128-3136","3140-3147","3251-3258")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        #legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

pdf("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/FVB.Trim28_2904-3258_females.pdf", w = 12, h = 10)
trim_weights %>%
  filter(...5 == "f") %>%
  drop_na() %>%
  ggplot(aes(x = ...10, y = ...9, color = ...2)) +
  geom_point(aes(shape = ...4), size = 5) +
  labs(x = "Lean Mass (g)", y = "Fat mass (g)", title = "FVB.Trim28(D9) Female Mice from 2904-3258", shape="Genotype", colour="Litter") +
  scale_color_discrete(labels = c("2904-2917", "3033-3041","3042-3052","3128-3136","3140-3147","3251-3258")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        #legend.position = "none", ## all these are just sixe of texts, position of title
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x= element_text(size = 22),
        axis.text.y= element_text(size = 22), 
        legend.title=element_text(size=20), 
        legend.text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

trim_weights_male <- trim_weights %>%
  drop_na() %>%
  filter(...5 == "m") 

trim_weights_female <- trim_weights %>%
  drop_na() %>%
  filter(...5 == "f")

pdf("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/FVB.Trim28_2904-3258_all_barplot.pdf", w = 12, h = 10)
test <- trim_weights %>% 
  drop_na() %>%
  group_by(...2, ...4) %>%
  summarise(mean_fat = mean(...9)) %>%
  rename( "litter" = "...2") %>%
  rename("Genotype" = "...4") 
  #filter(litter != "4") %>%
  ggplot(aes(x = litter, y = mean_fat, fill = Genotype)) +
  geom_bar(stat="identity", position = position_dodge2(preserve = ("single")), alpha=0.5) +
  geom_point(data = trim_weights, aes(y = ...9, fill =...4, color = ...4), show.legend = FALSE) +
  labs(x = "Litter", y = "Fat mass (g)", title = "FVB.Trim28(D9) Male/Female Mice from 2904-3258") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        #axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white')) 
dev.off()

pdf("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/FVB.Trim28_2904-3258_males_barplot.pdf", w = 12, h = 10)
trim_weights %>% 
  drop_na() %>%
  filter(...5 == "m") %>%
  group_by(...2, ...4) %>%
  summarise(mean_fat = mean(...9)) %>%
  rename( "litter" = "...2") %>%
  rename("Genotype" = "...4") %>%
  #filter(litter != "4") %>%
  ggplot(aes(x = litter, y = mean_fat, fill = Genotype)) +
  geom_bar(stat="identity", position = position_dodge2(preserve = ("single")), alpha=0.5) +
  geom_point(data = trim_weights_male, aes(x = ...2, y = ...9, fill =...4, color = ...4), show.legend = FALSE) +
  labs(x = "Litter", y = "Fat mass (g)", title = "FVB.Trim28(D9) Male Mice from 2904-3258") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        #axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'))
dev.off()

pdf("/Users/Warren.Sink/Desktop/OneDrive-Van-Andel-Institute/FVB.Trim28_2904-3258_females_barplot.pdf", w = 12, h = 10)
trim_weights %>% 
  drop_na() %>%
  filter(...5 == "f") %>%
  group_by(...2, ...4) %>%
  summarise(mean_fat = mean(...9)) %>%
  rename( "litter" = "...2") %>%
  rename("Genotype" = "...4") %>%
  #filter(litter != "4") %>%
  ggplot(aes(x = litter, y = mean_fat, fill = Genotype)) +
  geom_bar(stat="identity", position = position_dodge2(preserve = ("single")), alpha=0.5) +
  geom_point(data = trim_weights_female, aes(x = ...2, y = ...9, fill =...4, color = ...4), show.legend = FALSE) +
  labs(x = "Litter", y = "Fat mass (g)", title = "FVB.Trim28(D9) Female Mice from 2904-3258") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        #axis.text.y= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'))
dev.off()


ggplot(test, aes(x = litter, y = mean_fat, fill = Genotype)) +
  geom_bar(stat="identity", position = position_dodge()) +
  labs(x = "Litter", y = "Fat mass (g)", title = "FVB.Trim28(D9) Male/Female Mice from 2904-3258") +
  theme(plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),  
        axis.text.x= element_text(size = 15, angle = 0),
        axis.text.y= element_blank())
  
  
trim_weights %>% 
  group_by(...2) %>%
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





