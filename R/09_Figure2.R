# Build network complexity vs plant health figure for manuscript

# setup ####
# detach("package:GGally", unload = TRUE);detach("package:ggstats", unload = TRUE)
# detach("package:ggplot2",unload = TRUE)
library(ggplot2,lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.3/")
library(dplyr); library(stringr); library(purrr); library(tidyr);library(readr)
packageVersion('ggplot2') # should be 3.4.?

library(janitor)
library(patchwork)
library(phyloseq)
source("./R/palettes.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]
inoc_colors <- c("#4cbfe6","#2443f0",
                 "#f0c424","#db7e04",
                 "#cc6866","#9e0402",
                 "gray20")
'%ni%' <- Negate('%in%')
# data ####

# network data

# grouped
nets <- readRDS("./Output/network_attributes_grouped.RDS")
fung_nets <- nets %>% dplyr::filter(domain == "fungi")
bact_nets <- nets %>% dplyr::filter(domain == "bacteria")

# for each plant
all_plants_fung_nets <- readRDS("./Output/all_plants_fungi_network_stats.RDS")
all_plants_bact_nets <- readRDS("./Output/all_plants_bacteria_network_stats.RDS")

# plant growth data
plants <- 
  readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS") %>% 
  microbiome::meta() %>% 
  dplyr::select(unique_id,sample_name,block,leaf_number,wilting_scale,species,
                final_root_dm,shoot_dm,height,
                burn_frequency = fire_freq,
                soil_inoculum = inoculum,
                drought_trt = drought) %>% 
  mutate(host = case_when(species == "GrandFir" ~ "Abies grandis",
                          species == "Snowbrush" ~ "Ceanothus velutinus"),
         total_dry_mass = final_root_dm + shoot_dm)
    
# confirm that scaled values are well correlated with raw values
sb <- plants %>% 
  filter(host == "Ceanothus velutinus")
gf <- plants %>% 
  filter(host == "Abies grandis")

sb %>% 
  mutate(scaled_leaf_number = scale(leaf_number),
         scaled_drymass = scale(total_dry_mass),
         scaled_height = scale(height)) %>% 
  mutate(composite_score = (scaled_leaf_number + scaled_drymass + scaled_height) / 3 ) %>% 
  pivot_longer(starts_with("scaled_"),names_to = "variable") %>% 
  ggplot2::ggplot(aes(x=composite_score,y=value,color=variable)) +
  geom_smooth()


x <- sb %>% 
  mutate(scaled_leaf_number = scale(leaf_number) %>% as.numeric,
         scaled_drymass = scale(total_dry_mass) %>% as.numeric,
         scaled_height = scale(height) %>% as.numeric)
y <- x[,grep("leaf_number|dry_mass|height|drymass",names(x))] %>% janitor::clean_names()

GGally::ggpairs(y)


# plants <- read_csv("./Data/MasterMeasurement_Sep_Sheet1.csv") %>% clean_names()
# plants <- 
#   plants %>% 
#   dplyr::select(all_of(c("unique_id","block","unique_id_2","block_2","seed_species","drought_trt","soil_inoculum","burn_frequency")),
#                 contains("leaf_num")) %>% 
#   pivot_longer(contains("leaf_num"),names_to = "month",values_to = "leaf_number") %>% 
#   mutate(month = month %>% str_remove("tru.*_leaf_num_") %>% str_to_sentence() %>% str_replace("June","Jun")) %>% 
#   mutate(month = factor(month, levels = c("Nov","Dec","Jan","Feb","Mar","May","Jun","Aug"))) %>% 
#   # mutate(plantgrowth = scale(leaf_number)[,1]) %>% 
#   mutate(host = case_when(seed_species == "GrandFir" ~ "Abies grandis",
#                           seed_species == "Snowbrush" ~ "Ceanothus velutinus"))

sb <- 
sb %>% 
  mutate(scaled_leaf_number = scale(leaf_number),
         scaled_drymass = scale(total_dry_mass),
         scaled_height = scale(height)) %>% 
  mutate(composite_score = (scaled_leaf_number + scaled_drymass + scaled_height) / 3 )

gf <- 
  gf %>% 
  mutate(scaled_leaf_number = scale(leaf_number),
         scaled_drymass = scale(total_dry_mass),
         scaled_height = scale(height)) %>% 
  mutate(composite_score = (scaled_leaf_number + scaled_drymass + scaled_height) / 3 )


plants <- rbind(gf,sb)
plants$composite_score

# scaled_plantgrowth <-
#   plants %>%
#   group_by(host) %>%
#   reframe(plantgrowth = scale(leaf_number)) %>%
#   pluck("plantgrowth")
# plants$scaled_plantgrowth <- c(scaled_plantgrowth) %>% as.numeric()


# ps
# plant_means <- 
#   plants %>% 
#   group_by(host,drought_trt,month) %>% 
#   summarize(plantgrowth = mean(scale(leaf_number)[,1]))

plants$inoc_source <-   
  ifelse(grepl("_W",plants$soil_inoculum),plants$soil_inoculum,"_Sterile") %>% 
  str_split("_") %>% 
  map_chr(2) %>% 
  str_remove("W")
# plants <- plants %>% dplyr::filter(month == "Aug" & !is.na(unique_id))

plants <- 
plants %>% 
  dplyr::mutate(fung_sample_name = unique_id %>% str_replace("Q","F-R-"),
                bact_sample_name = unique_id %>% str_replace("Q","R-"))



# read in physeq objects
inoc_16s <- readRDS("./Output/phyloseq_objects/16S_inoculum_samples_clean_phyloseq_object.RDS")
inoc_its <- readRDS("./Output/phyloseq_objects/ITS_inoculum_samples_clean_phyloseq_object.RDS")
diversity <- estimate_richness(inoc_16s,measures=c("Observed","Shannon")) %>% 
  mutate(inoc_id=row.names(.),amplicon='bacteria') %>% 
  full_join(
    estimate_richness(inoc_its,measures=c("Observed","Shannon")) %>% 
      mutate(inoc_id=row.names(.),amplicon='fungi')
    
  ) %>% 
  mutate(inoc_source = str_sub(inoc_id,start=-1) %>% as.character) %>% 
  group_by(inoc_source) %>% 
  summarize(Shannon = mean(Shannon,na.rm=TRUE),
            Observed = mean(Observed,na.rm=TRUE))

plants %>% full_join(diversity,by="inoc_source") %>% 
  pivot_longer(all_of(c("Shannon","Observed")),
               names_to = "measure",
               values_to = "value") %>% 
  ggplot(aes(y=composite_score,x=value,color=inoc_source)) +
    geom_point() +
    facet_wrap(~measure*host,scales='free')

plants %>% 
  full_join(diversity,by="inoc_source") %>%
  glm(data= .,
      formula = composite_score ~ Shannon * host * drought_trt) %>% 
  summary
# No, alpha diversity of inoculum has no relationship with plant growth

# Try with diversity of root communities
root_16s <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS")
root_its <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")

diversity <- 
estimate_richness(root_16s,measures=c("Observed","Shannon")) %>% 
  mutate(sample_name=row.names(.),amplicon='bacteria') %>% 
  full_join(
    estimate_richness(root_its,measures=c("Observed","Shannon")) %>% 
      mutate(sample_name=row.names(.),amplicon='fungi')
    
  ) %>% 
  mutate(sample_name = sample_name %>% str_replace_all("\\.","-"))
diversity$sample_name2 <- diversity$sample_name
diversity <- 
plants %>% 
  pivot_longer(all_of(c("bact_sample_name","fung_sample_name")),
               values_to = "sample_name2") %>% 
  full_join(diversity,by="sample_name2")

lmerTest::lmer(data=diversity,
    formula = composite_score ~ Shannon * amplicon * host * drought_trt + (1|block)) %>% 
  summary

diversity %>% 
  ggplot(aes(x=Shannon,y=composite_score,color=drought_trt)) +
  geom_point() + geom_smooth(method='lm') +
  facet_wrap(~amplicon*host,scales = 'free_x') +
  scale_color_manual(values=drought_colors)

diversity %>% 
  glm(data=.,
      formula=composite_score ~ inoc_source) %>% 
  summary()

diversity %>% 
  mutate(moisture=case_when(drought_trt == "D" ~ "Low",
                            drought_trt == "ND" ~"High") %>% 
           factor(levels=c("Low","High"))) %>% 
  ggplot(aes(x=inoc_source,y=(composite_score),fill=moisture)) +
  geom_boxplot() +
  facet_wrap(~host) +
  labs(x="Inoculum source",y="Scaled\nplant growth") +
  scale_fill_manual(values=drought_colors) +
  theme_bw() +
  theme(axis.title.x = element_text(face='bold',size=18),
        axis.text.x = element_text(angle=60,face='bold',size=14,hjust=1,vjust=1),
        axis.text.y = element_text(face='bold',size=14),
        axis.title.y = element_text(face='bold',size=16),
        plot.title = element_text(face='bold.italic',size=18,hjust=.5),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face='bold.italic',size=18)
  )


# pool network info with plant scores

# need a function/loop to extract all this from the lists
selected_attr <- c("n_vertices","n_edges","global_effic","clustering_coeficient","max_degree")

fung_newlist <- list()
for(i in selected_attr){
  fung_newlist[[i]] <- map_dbl(all_plants_fung_nets,function(x){pluck(x,i)})
}
fung_newlist <- fung_newlist %>% as.data.frame()
fung_newlist$fung_sample_name <- row.names(fung_newlist)
fung_newlist$amplicon <- "fungi"

bact_newlist <- list()
for(i in selected_attr){
  bact_newlist[[i]] <- map_dbl(all_plants_bact_nets,function(x){pluck(x,i)})
}
bact_newlist <- bact_newlist %>% as.data.frame()
bact_newlist$bact_sample_name <- row.names(bact_newlist)
bact_newlist$amplicon <- "bacteria"

# Combine network stats with plant responses and treatments ####

fung <- plants %>% full_join(fung_newlist)
bact <- plants %>% full_join(bact_newlist)

# Come up with 'composite network complexity' metric (mean scaled network stat score)

# Invert appropriate metrics
fung$clustering_coeficient <- 1/fung$clustering_coeficient
fung$global_effic <- 1/fung$global_effic

fung$mean_scaled_network_score <- 
fung %>% 
  dplyr::select(all_of(selected_attr)) %>% 
  mutate(across(everything(),scale)) %>% 
  rowMeans(na.rm=TRUE)

bact$mean_scaled_network_score <- 
  bact %>% 
  dplyr::select(all_of(selected_attr)) %>% 
  mutate(across(everything(),scale)) %>% 
  rowMeans(na.rm=TRUE)


# Build full data frame
dat <- full_join(fung, bact)
dat <- dat %>% 
       mutate(moisture = case_when(drought_trt == "D" ~ "Low", drought_trt == "ND" ~ "High"),
              amplicon = factor(amplicon,levels=c("fungi","bacteria")))



# PLOT ####
names(dat)
dat %>% 
  dplyr::filter(!is.na(amplicon)) %>% 
  ggplot(aes(x=mean_scaled_network_score,
             y=scaled_plantgrowth,
             color=moisture)) +
  geom_point() + 
  geom_smooth(method='lm') +
  facet_wrap(~amplicon*host,scales='free')

# plotting function
max_y <- max(dat$scaled_plantgrowth,na.rm=TRUE)
min_y <- min(dat$scaled_plantgrowth,na.rm=TRUE)
plot_fun <- 
  function(ampvar,hostvar){
  dat %>% 
    dplyr::filter(amplicon == ampvar & host == hostvar) %>% 
    ggplot(aes(x=mean_scaled_network_score,
               y=scaled_plantgrowth,
               color=moisture)) +
    geom_smooth(method='lm',fill='gray80',linetype=0) +
    geom_point(size=3) +
    geom_smooth(method='lm',se=FALSE) +
    # facet_wrap(~host,scales='free_x') +
    labs(color="Moisture",
         x = "Scaled composite network complexity",
         y = "Scaled\nplant growth") +
    scale_color_manual(values=drought_colors) +
    coord_cartesian(ylim = c(min_y,max_y)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text = element_text(face='bold.italic'),
          legend.title = element_text(face='bold',size=16),
          legend.text = element_text(face='bold',size=14),
          axis.title.y = element_text(face='bold',size=18),
          axis.title.x = element_text(face='bold',size=18),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=14,face='bold'),
          plot.title = element_text(face='bold.italic',size=18,hjust=.5))
}

# run base plots for all 4 groups
fungal_network_plot_gf <- 
  plot_fun(ampvar = "fungi",hostvar = "Abies grandis") + 
  ggtitle("Abies grandis") + 
  labs(y="Fungi") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face='bold',size=14))
fungal_network_plot_sb <- 
  plot_fun(ampvar = "fungi",hostvar = "Ceanothus velutinus") + 
  ggtitle("Ceanothus velutinus") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face='bold',size=14))
bacterial_network_plot_gf <- 
  plot_fun(ampvar = "bacteria",hostvar = "Abies grandis") +
  labs(y="Bacteria") +
  theme(axis.text.x = element_text(face='bold',size=14))
bacterial_network_plot_sb <- 
  plot_fun(ampvar = "bacteria",hostvar = "Ceanothus velutinus") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face='bold',size=14))

# patch together...
patch <- (fungal_network_plot_gf | fungal_network_plot_sb) /
  (bacterial_network_plot_gf | bacterial_network_plot_sb) +
  plot_layout(guides='collect') & xlab(NULL)
patch
patch2 <- 
wrap_elements(panel = patch) +
  labs(tag = "Scaled composite network complexity") +
  theme(
    plot.tag = element_text(size = 18,face='bold'),
    plot.tag.position = "bottom"
  )
patch2 


# build 3-panel plot for manuscript
gf_fungal_mutualist_plot <- readRDS("./Output/figs/ITS_gf_mutualist_plot.RDS") + 
  labs(x="\nProportion of putative\nfungal mutualists",
       y="Scaled\nplant growth",
       title = "Abies grandis") + 
  theme(legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        legend.position = 'none', 
        plot.title = element_text(face='bold.italic',size=18), 
        axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14)) + 
  geom_point(size=3) + 
  geom_smooth(method = 'lm',color='black')

gf_bacterial_mutualist_plot <- readRDS("./Output/figs/16S_gf_mutualist_plot.RDS") + 
  labs(x='\nProportion of putative\nbacterial mutualists',
       y="Scaled\nplant growth",
       title="Abies grandis",
       color="Moisture") + 
  theme(legend.position = 'none', 
        plot.title = element_text(face='bold.italic',size=18), 
        axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face = 'bold',size=14)) + 
  geom_point(size=3) + 
  geom_smooth(method='lm',color='black')

sb_bacterial_pathogen_plot <- readRDS("./Output/figs/16S_sb_pathogen_plot.RDS") + 
  labs(x='\nProportion of putative\nbacterial pathogens',
       color="Moisture",
       title="Ceanothus velutinus") + 
  theme(plot.title = element_text(face='bold.italic',size=18), 
        axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face = 'bold',size=14)) + 
  geom_point(size=3) + 
  geom_smooth(method='lm')

pw <- 
  (gf_fungal_mutualist_plot + gf_bacterial_mutualist_plot + sb_bacterial_pathogen_plot) &
  ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
guild_plots <- 
  wrap_elements(pw) +
  labs(tag = "Scaled\nplant growth\n") +
  theme(
    plot.tag = element_text(angle = 90,face='bold',size=18),
    plot.tag.position = "left"
  ) +
  plot_layout(guides = 'collect',axis_title="collect")
# This ensures that italics remain in saved versions of the plots
showtext::showtext_auto()
showtext::showtext_opts(dpi = 140)
guild_plots

saveRDS(guild_plots,"./Output/figs/ALL_guild_plots.RDS")
ggsave("./Output/figs/ALL_guild_plots.tiff",device = 'tiff',width = 24,height = 8,dpi=400)
ggsave("./Output/figs/ALL_guild_plots.png",device = 'png',width = 24,height = 8,dpi=400)
ggsave("./Output/figs/ALL_guild_plots.pdf",device = 'pdf',width = 24,height = 8,dpi=600)


patch2 / guild_plots + plot_layout(guides='collect')

patch2
ggsave(plot = patch2,"./Output/figs/ALL_network_by_plant_growth_plots.pdf",device = 'pdf',
       width = 8,height = 6,dpi=600)
# maybe group by inoc source, like before and do a boxplot after all !?
new <- dat %>%
  dplyr::filter(!is.na(amplicon)) %>%
  group_by(inoc_source, moisture,host,amplicon) %>%
  summarize(new = median(mean_scaled_network_score,na.rm=TRUE))

gf_hi_fung <- 
new %>% 
  dplyr::filter(host == "Abies grandis") %>% 
  dplyr::filter(moisture == "High") %>% 
  dplyr::filter(amplicon == "fungi") %>% 
  arrange(new) %>% 
  pluck('inoc_source')
  
dat %>% 
  dplyr::filter(host == "Abies grandis") %>% 
  dplyr::filter(moisture == "High") %>% 
  dplyr::filter(amplicon == "fungi") %>% 
  mutate(inoc_source = factor(inoc_source,levels = new %>% 
                                dplyr::filter(host == "Abies grandis") %>% 
                                dplyr::filter(moisture == "High") %>% 
                                dplyr::filter(amplicon == "fungi") %>% 
                                arrange(new) %>% 
                                pluck('inoc_source'))) %>% 
  ggplot(aes(x=inoc_source,y=scaled_plantgrowth)) +
  geom_boxplot() +
  theme_bw()

rm(a)
rm(m)
rm(h)



plot_fun2 <- function(h,m,a){
  
  filt <- dat %>% 
    dplyr::filter(host == h) %>% 
    dplyr::filter(moisture == m) %>% 
    dplyr::filter(amplicon == a)
  filt$amplicon <- filt$amplicon %>% as.character()
  ord <- filt %>% 
    group_by(inoc_source) %>% 
    summarize(med = median(scaled_plantgrowth,na.rm=TRUE)) %>% 
    arrange(med) %>% 
    pluck("inoc_source")
  meds <- filt %>% 
    group_by(inoc_source) %>% 
    summarize(med = median(scaled_plantgrowth,na.rm=TRUE)) %>% 
    arrange(med) %>% 
    pluck("med")


    
  filt %>% 
    mutate(inoc_source = factor(inoc_source,levels = ord)) %>% 
    ggplot(aes(x=inoc_source,y=scaled_plantgrowth)) +
    geom_boxplot(alpha=.5,color=ifelse(m=="Low",drought_colors[1],drought_colors[2])) +
    geom_path(data=data.frame(ord,meds),aes(group=1,x=ord,y=meds),linewidth=1,color=ifelse(m=="Low",drought_colors[1],drought_colors[2])) +
    labs(title=h,y="Scaled\nplant growth") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(face='bold',size=14),
          axis.title.y = element_text(face='bold',size=16),
          plot.title = element_text(face='bold.italic',size=18,hjust=.5),
          panel.grid = element_blank()
          )
}


p1 <- plot_fun2(h="Abies grandis",m="Low",a="fungi") + theme(axis.title.y = element_blank(),plot.title = element_blank())
p2 <- plot_fun2(h="Abies grandis",m="High",a="fungi") + theme(axis.title.y = element_blank())
p3 <- plot_fun2(h="Ceanothus velutinus",m="Low",a="fungi") + theme(axis.title.y = element_blank(),plot.title = element_blank())
p4 <- plot_fun2(h="Ceanothus velutinus",m="High",a="fungi") + theme(axis.title.y = element_blank())


pw2 <- (p2 | p4) / (p1 | p3)

guild_plots2 <- 
  wrap_elements(pw2) +
  labs(tag = "Scaled plant growth") +
  theme(
    plot.tag = element_text(angle = 90,face='bold',size=18),
    plot.tag.position = "left"
  ) +
  plot_layout(guides = 'collect',axis_title="collect")
guild_plots2

ggsave("./Output/figs/network_plots_draft_boxplots.pdf",width = 8,height = 6,dpi=600)
# ggsave(guild_plots2,"./Output/figs/network_plots_draft_boxplots.pdf",width = 8,height = 6,dpi=600,device = 'pdf')


dat %>%
  dplyr::filter(!is.na(amplicon)) %>%
  full_join(new) %>%
  ggplot(aes(x=new,y=scaled_plantgrowth,color=moisture)) +
  geom_smooth(method='lm') +
  facet_wrap(~amplicon*host,scales='free')



# step <- 
# glm(data=dat,
#     formula = scaled_plantgrowth ~ mean_scaled_network_score * host * amplicon * moisture) %>% 
#   MASS::stepAIC()
# 
# glm(data=dat,
#     formula = step$formula) %>% summary
# 
# dat %>% 
#   dplyr::filter(amplicon == "fungi" & host == "Ceanothus velutinus" & moisture == "Low") %>% 
#   glm(data=.,formula = scaled_plantgrowth ~ mean_scaled_network_score) %>% 
#   summary

# fungal_network_plot_gf <- 
#   dat %>% 
#   dplyr::filter(amplicon == "fungi" & host == "Abies grandis") %>% 
#   ggplot(aes(x=mean_scaled_network_score,
#              y=scaled_plantgrowth,
#              color=moisture)) +
#   geom_smooth(method='lm',fill='gray80',linetype=0) +
#   geom_point(size=3) +
#   geom_smooth(method='lm',se=FALSE) +
#   facet_wrap(~host,scales='free_x') +
#   labs(color="Moisture",
#        x = "Scaled composite network complexity",
#        y = "Scaled\nplant growth") +
#   scale_color_manual(values=drought_colors) +
#   theme(strip.background = element_blank(),
#         strip.text = element_text(face='bold.italic'),
#         legend.title = element_text(face='bold',size=16),
#         legend.text = element_text(face='bold',size=14),
#         axis.title.y = element_text(face='bold',size=18),
#         axis.title.x = element_text(face='bold',size=18),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(size=14,face='bold'),
#         plot.title = element_text(face='bold.italic',size=18,hjust=.5))
# fungal_network_plot_gf



# build summarized table
grouped_summary <- 
dat %>% 
  group_by(amplicon, inoc_source, host, moisture) %>% 
  dplyr::select(amplicon, inoc_source, host, moisture, all_of(selected_attr),mean_scaled_network_score,scaled_plantgrowth) %>% 
  dplyr::filter(!is.na(amplicon)) %>% 
  summarize_at(c(selected_attr,"mean_scaled_network_score","scaled_plantgrowth"),
               list(mean,sd),na.rm=TRUE)
#rename sanely
names(grouped_summary) <- 
names(grouped_summary) %>% 
  str_replace("fn1","mean") %>% 
  str_replace("fn2","sd")

grouped_summary %>% names

grouped_summary %>% ggplot(aes(x=scaled_plantgrowth_mean)) + geom_histogram()
grouped_summary %>% 
  ggplot(aes(x=mean_scaled_network_score_mean,
             y=scaled_plantgrowth_mean,
             ymin=scaled_plantgrowth_mean - scaled_plantgrowth_sd,
             ymax=scaled_plantgrowth_mean + scaled_plantgrowth_sd,
             color=moisture)) +
  geom_errorbar() +
  geom_smooth(method='lm') +
  facet_wrap(~host*amplicon,scales = 'free')

