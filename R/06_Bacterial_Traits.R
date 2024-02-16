# -----------------------------------------------------------------------------#
# Pulling and processinf bacterial traits from BacDive database 
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.42.0
#                     patchwork v 1.1.2
#                     BacDive v 0.8.0
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
#  install.packages("BacDive", repos="http://R-Forge.R-project.org")

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")
source("./R/palettes.R")

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
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_pathogen * drought * inoculum_site + (1|block))
summary(grandfir_pathogen_glm)
saveRDS(grandfir_pathogen_glm,"./Output/16S_Pathogen_Model_GrandFir.RDS")

snowbrush_pathogen_glm <- 
  pathogen_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_pathogen * drought * inoculum_site + (1|block))
summary(snowbrush_pathogen_glm)
saveRDS(snowbrush_pathogen_glm,"./Output/16S_Pathogen_Model_Snowbrush.RDS")


full_guild_model_df <- 
  clean_model_df(snowbrush_pathogen_glm) %>% mutate(species="Ceanothus") %>% 
  full_join(clean_model_df(grandfir_pathogen_glm) %>% mutate(species="Abies")) %>% 
  dplyr::filter(effect=="fixed") %>% 
  select(-group)

saveRDS(full_guild_model_df,"./Output/16S_Guild_Model_Output.RDS")


readRDS("./Output/16S_Guild_Model_Output.RDS")
readRDS("./Output/ITS_Guild_Model_Output.RDS")
