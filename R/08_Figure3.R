library(tidyverse)
library(patchwork)
library(janitor)

source("./R/palettes.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]
inoc_colors <- c("#6ae802","#1b9e2c",
                 "#ede502","#c7b602",
                 "#e60712","#f2525a",
                 "gray30")


# load plant time series data
plants <- read_csv("./Data/MasterMeasurement_Sep_Sheet1.csv") %>% clean_names()

# load ordination plots
grandfir_bact_ord_plot <- readRDS("./Output/grandfir_bact_ord_plot.RDS") + theme(axis.title = element_blank()) + labs(y="Bacteria") + theme(axis.title.y = element_text(face='bold',size=14))
snowbrush_bact_ord_plot <- readRDS("./Output/snowbrush_bact_ord_plot.RDS") + theme(axis.title = element_blank()) 
grandfir_fungi_ord_plot <- readRDS("./Output/grandfir_fungi_ord_plot.RDS")  + theme(axis.title = element_blank()) + labs(y="Fungi") + theme(axis.title.y = element_text(face='bold',size=14))
snowbrush_fungi_ord_plot <- readRDS("./Output/snowbrush_fungi_ord_plot.RDS") + theme(axis.title = element_blank())
snowbrush_18S_ord_plot <- readRDS("./Output/snowbrush_18S_ord_plot.RDS") + theme(axis.title = element_blank())


# clean plant data
names(plants)
plants <- 
plants %>% 
  dplyr::select(all_of(c("unique_id","block","unique_id_2","block_2","seed_species","drought_trt","soil_inoculum","burn_frequency")),
                contains("leaf_num")) %>% 
  pivot_longer(contains("leaf_num"),names_to = "month",values_to = "leaf_number") %>% 
  mutate(month = month %>% str_remove("tru.*_leaf_num_") %>% str_to_sentence() %>% str_replace("June","Jun")) %>% 
  mutate(month = factor(month, levels = c("Nov","Dec","Jan","Feb","Mar","May","Jun","Aug"))) %>% 
  # mutate(plantgrowth = scale(leaf_number)[,1]) %>% 
  mutate(host = case_when(seed_species == "GrandFir" ~ "Abies grandis",
                          seed_species == "Snowbrush" ~ "Ceanothus velutinus"))

scaled_plantgrowth <- 
plants %>% 
  group_by(host) %>% 
  summarize(plantgrowth = scale(leaf_number)) %>% 
  pluck("plantgrowth")
plants$scaled_plantgrowth <- c(scaled_plantgrowth) %>% as.numeric()

plant_means <- 
plants %>% 
  group_by(host,drought_trt,month) %>% 
  summarize(plantgrowth = mean(scale(leaf_number)[,1]))

plants$inoc_source <-   
ifelse(grepl("_W",plants$soil_inoculum),plants$soil_inoculum,"_Sterile") %>% 
  str_split("_") %>% 
  map_chr(2) %>% 
  str_remove("W")


# plot plant growth over time
plants %>% 
  dplyr::filter(!is.na(host)) %>% 
  ggplot(aes(x=month,y=scaled_plantgrowth,color=drought_trt)) +
  geom_jitter(size=3,width = .2,alpha=.25) +
  geom_boxplot(alpha=.2,size=1.5,fill='black') +
  facet_wrap(~host,ncol = 2) +
  scale_color_manual(values = drought_colors,labels=c("Low","High")) +
  labs(color="Moisture",x="\nTimepoint",y="Scaled plant growth") +
  theme_bw() +
  theme(strip.text = element_text(face='bold.italic',size=18),
        strip.background = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5))
ggsave("./Output/figs/plantgrowth_time_series_host.png",dpi=300,height = 6,width = 10)


plants %>% 
  dplyr::filter(!is.na(host)) %>% 
  ggplot(aes(x=month,y=scaled_plantgrowth,color=drought_trt)) +
  geom_jitter(size=3,width = .2,alpha=.25) +
  geom_boxplot(alpha=.2,size=1.5,fill='black') +
  facet_wrap(~inoc_source) +
  scale_color_manual(values = drought_colors,labels=c("Low","High")) +
  labs(color="Moisture",x="\nTimepoint",y="Scaled plant growth") +
  theme_bw() +
  theme(strip.text = element_text(face='bold.italic',size=18),
        strip.background = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5))
ggsave("./Output/figs/plantgrowth_time_series_inoc.png",dpi=300,height = 6,width = 10)




ordplots <- 
(grandfir_bact_ord_plot + snowbrush_bact_ord_plot) / (grandfir_fungi_ord_plot + snowbrush_fungi_ord_plot) +
  plot_layout(guides = 'collect') # + plot_annotation(tag_levels = "A")
