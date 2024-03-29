# -----------------------------------------------------------------------------#
# Pulling and processinf bacterial traits from BacDive database 
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.42.0
#                     patchwork v 1.1.2
#                     BacDive v 0.8.0
#                     lmerTest v 3.1.3
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
#  install.packages("BacDive", repos="http://R-Forge.R-project.org")

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")
library(lmerTest); packageVersion("lmerTest")
source("./R/palettes.R")
source("./R/scale01.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]


clean_model_df <- function(x){
  broom.mixed::tidy(x) %>% 
    mutate(term=term %>% str_remove("indicator")) %>% 
    mutate(across(where(is.numeric),function(z){round(z,4)}))
}


# Random seed
set.seed(666)

readRenviron("./.Renviron")

# Data
bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")

# initialize bacdive connection (add your userid and password)
# BacDive URL: https://bacdive.dsmz.de/
bacdive <- BacDive::open_bacdive(username = Sys.getenv("BACDIVE_USER"),
                                 password = Sys.getenv("BACDIVE_PW"))


# get list of unique genera
genus_list <- bact@tax_table[,6] %>% 
  table()
genus_list <- genus_list %>% 
  as.data.frame() %>% 
  pluck(".") %>% 
  levels() 
names(genus_list) <- genus_list

# BacDive API ####

# run in for-loop to build giant database of all traits for our taxa
bact_trait_db <- list()

for(i in genus_list){
  x <- BacDive::request(object = bacdive,
                        query = genus_list[i],
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  # if no results found, return NA
  if(length(y) == 0){
    bact_trait_db[[i]] <- NA
  } else {
    # otherwise, get info on that id
    z <- BacDive::fetch(bacdive,y)
    bact_trait_db[[i]] <- z$results  
  }
}

# Find any that didn't work... and re-run (in case of database hang-ups)
for(i in genus_list[which(!genus_list %in% names(bact_trait_db))]){
  x <- BacDive::request(object = bacdive,
                        query = genus_list[i],
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  # if no results found, return NA
  if(length(y) == 0){
    bact_trait_db[[i]] <- NA
  } else {
    # otherwise, get info on that id
    z <- BacDive::fetch(bacdive,y)
    bact_trait_db[[i]] <- z$results  
  }
}

# Find any that didn't work... and re-run (in case of database hang-ups)
for(i in genus_list[which(!genus_list %in% names(bact_trait_db))]){
  x <- BacDive::request(object = bacdive,
                        query = genus_list[i],
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  # if no results found, return NA
  if(length(y) == 0){
    bact_trait_db[[i]] <- NA
  } else {
    # otherwise, get info on that id
    z <- BacDive::fetch(bacdive,y)
    bact_trait_db[[i]] <- z$results  
  }
}

# save output from previous slow steps
saveRDS(bact_trait_db,"./Output/16S_Bacterial_Trait_Database.RDS")
# reload point, for convenience
bact_trait_db <- readRDS("./Output/16S_Bacterial_Trait_Database.RDS")
        
# IDENTIFY PATHOGENIC GENERA ####
# This database is a highly nested and confusing pile of information
# access looks something like
bact_trait_db$Rhizobium$`132155`$`Physiology and metabolism`
# DB_OBJECT   Genus   AccessionID   Trait Aspect

# 
# # Best way night be a nested for-loop with lots of tests for missing data :(
# 
# genus_physiology <- list()
# for(genus in names(bact_trait_db)){
#   
#   all_strains <- bact_trait_db[[genus]]
#   
#   
#   
#   for(strain in names(all_strains)){
#     
#     strain_x <- all_strains[[strain]]
#     physiology_x <- strain_x[["Physiology and metabolism"]]
#     
#     
#     
#     if(length(physiology_x) < 1){
#       genus_physiology[[genus]] <- NA
#     }
#   }  
# }
# 
# 
# genus_name <- bact_trait_db %>% names
# 
# for(i in genus_name[310]){
#   bact_trait_db %>%
#     pluck(genus_name[310]) %>% 
#     # pluck(i) %>% 
#     names %>% 
#     # pluck(1) %>% 
#     # pluck("General") %>%
#     print()
# }
# 
# 3
# get names of all the traits
genus_trait_names <- 
map(bact_trait_db, function(x){
  x %>% pluck(1) %>% names
})


# Pluck keyword lists from each taxon
genus_general_keywords <- 
map(bact_trait_db, function(x){
  x %>% pluck(1) %>% pluck("General") %>% pluck("keywords")
}) %>% map(unlist)


# Search for "pathogen" in any field
pathogen_list <- c()
for(i in genus_list){
  pathogen_list[i] <- 
  bact_trait_db[[i]] %>% 
    unlist() %>% 
    grepl(pattern="pathogen", ignore.case = TRUE) %>% 
    any()
}
pathogen_list <- pathogen_list %>% which %>% names
# export
saveRDS(pathogen_list,"./Output/list_of_pathogenic_bacterial_genera.RDS")

# ADD TO PHYSEQ ####
ps <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")
tax_df <- ps@tax_table[,c("Genus","Species")] %>% as.data.frame()
tax_df2 <- data.frame(Genus = pathogen_list,Guild = "Pathogen")
# Add it to "Species" slot of tax_table
ps@tax_table[,"Species"] <- left_join(tax_df,tax_df2) %>% pluck("Guild")
colnames(ps@tax_table)[7] <- "Guild"

saveRDS(ps,"Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS")

# CALCULATE BACTERIAL GUILD PROPORTIONS ####

# find proportions of "pathogenic" genera
pathogen_proportions <- 
  bact %>% 
  tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(bact@tax_table[,6] %in% pathogen_list) %>% 
  sample_sums()

# build data frame for modeling
pathogen_df <- 
  microbiome::meta(bact) %>% 
  mutate(proportion_pathogen = pathogen_proportions)

# find proportions of "Frankia" & "Rhizobium" for 'selected mutualists'
# bact@tax_table[,6] %>% unname %>% unique %>% sort
mutualist_list <- c("Frankia","Rhizobium")
b <- bact %>% 
  tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)})
mutualist_proportions <- 
  b %>% 
  subset_taxa(b@tax_table[,6] %in% mutualist_list) %>% 
  sample_sums()
# build data frame for modeling
mutualist_df <- 
  microbiome::meta(bact) %>% 
  mutate(proportion_mutualist = mutualist_proportions)

guild_df <- full_join(pathogen_df,mutualist_df)

# REGRESSION PLOTS ####
grandfir_pathogen_plot <-
  pathogen_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_pathogen,y=value,color=drought)) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of pathogenic bacterial genera",y="Scaled/Centered Value",color="Drought") +
  scale_color_manual(values = drought_colors)
saveRDS(grandfir_pathogen_plot,"./Output/figs/16S_Pathogen_Plot_grandfir.RDS")


grandfir_pathogen_plot2 <-
  pathogen_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_pathogen,y=value,color=ordered(fire_freq,levels=c("0","1","3")))) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of pathogenic bacterial genera",y="Scaled/Centered Value",color="Fire frequency") +
  scale_color_manual(values = fire_colors)
saveRDS(grandfir_pathogen_plot2,"./Output/figs/16S_Pathogen_Plot_grandfir_by_fire.RDS")

snowbrush_pathogen_plot <-
  pathogen_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_pathogen,y=value,color=drought)) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of pathogenic bacterial genera",y="Scaled/Centered Value",color="Drought") +
  scale_color_manual(values = pal.discrete[c(2,7)])
saveRDS(snowbrush_pathogen_plot,"./Output/figs/16S_Pathogen_Plot_snowbrush.RDS")

snowbrush_pathogen_plot2 <-
  pathogen_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_pathogen,y=value,color=ordered(fire_freq,levels=c("0","1","3")))) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of pathogenic bacterial genera",y="Scaled/Centered Value",color="Drought") +
  scale_color_manual(values = fire_colors)
saveRDS(snowbrush_pathogen_plot2,"./Output/figs/16S_Pathogen_Plot_snowbrush_by_fire.RDS")


# MODELS ####
grandfir_pathogen_glm <- 
  pathogen_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  dplyr::select(all_of(c("leaf_number","shoot_dm","final_root_dm","proportion_pathogen","drought","block"))) %>% 
  mutate(across(c("leaf_number","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_pathogen * drought + (1|block))
summary(grandfir_pathogen_glm)
saveRDS(grandfir_pathogen_glm,"./Output/16S_Pathogen_Model_GrandFir.RDS")

snowbrush_pathogen_glm <- 
  pathogen_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_pathogen * drought + (1|block))
summary(snowbrush_pathogen_glm)
saveRDS(snowbrush_pathogen_glm,"./Output/16S_Pathogen_Model_Snowbrush.RDS")

# mutualist
grandfir_mutualist_glm <- 
  mutualist_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_mutualist * drought + (1|block))
summary(grandfir_mutualist_glm)
saveRDS(grandfir_mutualist_glm,"./Output/16S_mutualist_Model_GrandFir.RDS")

snowbrush_mutualist_glm <- 
  mutualist_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_mutualist * drought + (1|block))
summary(snowbrush_mutualist_glm)
saveRDS(snowbrush_mutualist_glm,"./Output/16S_mutualist_Model_Snowbrush.RDS")

full_guild_model_df <- 
  clean_model_df(snowbrush_pathogen_glm) %>% mutate(species="Ceanothus") %>% 
  full_join(clean_model_df(grandfir_pathogen_glm) %>% mutate(species="Abies")) %>% 
  full_join(clean_model_df(snowbrush_mutualist_glm) %>% mutate(species="Ceanothus")) %>% 
  full_join(clean_model_df(grandfir_mutualist_glm) %>% mutate(species="Abies")) %>% 
  dplyr::filter(effect=="fixed") %>% 
  select(-group)
full_guild_model_df


# PLOT ALL GUILD EFFECTS ####
full_guild_model_df <- full_guild_model_df %>% 
  mutate(PVAL = case_when(p.value == 0 ~ "P < 0.005",
                          TRUE ~ paste0("P = ",round(p.value,3) %>% as.character)))
gf <- guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(growth_response = scale(leaf_number))
sb <- guild_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(growth_response = scale(leaf_number))

full_guild_model_df %>% 
  dplyr::filter(p.value <= .05 & term != "(Intercept)")

gf_mutualist_plot <- 
gf %>% 
  mutate(Moisture = case_when(drought=="ND" ~ "High", drought == "D" ~ "Low")) %>% 
  ggplot(aes(x=proportion_mutualist,y=growth_response,color=Moisture)) +
  geom_point(alpha=.5,size=3) +
  geom_smooth(method='lm',se=FALSE,color='black') +
  annotate('text',
           x = .35,
           y = 3.5,
           label = "P < 0.005",
           fontface=2) +
  labs(x='',y="Plant growth response",title="Grand fir") +
  theme_bw() +
  theme(axis.title = element_text(face='bold',size=16),plot.title = element_text(face='bold',size=12,hjust=.5),
        legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12)) +
  scale_color_manual(values = drought_colors)
saveRDS(gf_mutualist_plot,"./Output/figs/16S_gf_mutualist_plot.RDS")

sb_mutualist_plot <- 
  sb %>% 
  mutate(Moisture = case_when(drought=="ND" ~ "High", drought == "D" ~ "Low")) %>% 
  ggplot(aes(x=proportion_mutualist,y=growth_response,color=Moisture)) +
  geom_point(alpha=.5,size=3) +
  geom_smooth(method='lm',se=FALSE) +
  # annotate('text',
  #          x = .6,
  #          y = 2,
  #          label = full_guild_model_df %>% 
  #            dplyr::filter(species == "Ceanothus" & term == "proportion_mutualist") %>%
  #            pluck("PVAL")) +
  labs(x="Proportion of mutualist ASVs",y="Plant growth response",title="Snowbrush") +
  theme_bw() +
  theme(axis.title = element_text(face='bold',size=16),plot.title = element_text(face='bold',size=12,hjust=.5))

gf_pathogen_plot <- 
  gf %>% 
  ggplot(aes(x=proportion_pathogen,y=growth_response)) +
  geom_point(alpha=.5,size=3) +
  geom_smooth(method='lm',se=FALSE) +
  # annotate('text',
  #          x = max(gf$proportion_pathogen) - .1,
  #          y = max(gf$growth_response),
  #          label = full_guild_model_df %>% 
  #            dplyr::filter(species == "Abies" & term == "proportion_pathogen") %>%
  #            pluck("PVAL")) +
  labs(x='',y="",title="Grand fir") +
  theme_bw() +
  theme(axis.title = element_text(face='bold',size=16),plot.title = element_text(face='bold',size=12,hjust=.5))

sb_pathogen_plot <- 
  sb %>% 
  mutate(Moisture = case_when(drought=="ND" ~ "High", drought == "D" ~ "Low")) %>% 
  ggplot(aes(x=proportion_pathogen,y=growth_response,color=Moisture)) +
  geom_point(alpha=.5,size=3) +
  geom_smooth(method='lm',se=FALSE) +
  annotate('text',
           x = .7,
           y = 2,
           label = "P = 0.001",
           fontface=2) +
  labs(color = "Drought") +
  labs(x="Proportion of pathogen ASVs",y="",title="Snowbrush") +
  theme_bw() +
  theme(axis.title = element_text(face='bold',size=16),
        legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12),plot.title = element_text(face='bold',size=12,hjust=.5)) +
  scale_color_manual(values = drought_colors)
saveRDS(sb_pathogen_plot,"./Output/figs/16S_sb_pathogen_plot.RDS")

pw <- 
(gf_mutualist_plot + gf_pathogen_plot) / (sb_mutualist_plot + sb_pathogen_plot) &
  ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
guild_plots <- 
wrap_elements(pw) +
  labs(tag = "Scaled plant growth measure\n") +
  theme(
    plot.tag = element_text(angle = 90,face='bold',size=16),
    plot.tag.position = "left"
  ) +
  plot_layout(guides = 'collect',axis_title="collect") +
  plot_annotation()
  
saveRDS(guild_plots,"./Output/figs/16S_Guild_Plots.RDS")

