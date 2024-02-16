# -----------------------------------------------------------------------------#
# Pulling and processinf fungal traits from FungalTraits database 
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 2.0.0
#                     phyloseq v 1.46.0
#                     FUNGuildR v 0.2.0.9
#                     fungaltraits v 0.0.3
#                     broom v 1.0.5
#                     lmerTest v 3.1.3
#                     broom.mixed 0.2.9.4
# -----------------------------------------------------------------------------#

# SETUP ####

# load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(FUNGuildR); packageVersion("FUNGuildR")
library(fungaltraits); packageVersion("fungaltraits")
library(broom); packageVersion("broom")
library(lmerTest); packageVersion("lmerTest")
library(broom.mixed); packageVersion("broom.mixed")

# library(vegan); packageVersion("vegan")
# library(microbiome); packageVersion("microbiome")

source("./R/palettes.R")

drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]

# functions
clean_model_df <- function(x){
  broom.mixed::tidy(x) %>% 
    mutate(term=term %>% str_remove("indicator")) %>% 
    mutate(across(where(is.numeric),function(z){round(z,4)}))
}
'%ni%' <- Negate('%in%')

# load data
fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS")
fung@sam_data %>% row.names()

# BUILD DATABASES ####

## FunGuild ####

# make data frame of taxonomy, formatted for FUNGuildR package
tax_df <- 
  data.frame(
    Taxonomy = paste(
      tax_table(fung)[,1],
      tax_table(fung)[,2],
      tax_table(fung)[,3],
      tax_table(fung)[,4],
      tax_table(fung)[,5],
      tax_table(fung)[,6],
      tax_table(fung)[,7],
      sep=";"
    )
  )

# download current FunGuild database and save as RDS
guild_db <- FUNGuildR::get_funguild_db()

# save guild db as RDS
saveRDS(guild_db, "./Taxonomy/Funguild_Database.RDS")

# assign guild to fungal ASV taxonomy
guilds <- FUNGuildR::funguild_assign(tax_df)

# examine guild assignments
guilds %>% head

# add sequence column to guild assignments
guilds$ASV_seq <- phyloseq::taxa_names(fung)

# save assignments as RDS
saveRDS(guilds,"./Output/ITS_Fungal_Guild_Assignments.RDS")

# add guild assignments to "pseudo-tax_table"
# (this replaces and renames the "Kingdom" column)
fung@tax_table[,1] <- guilds$guild
attributes(fung@tax_table)$dimnames[[2]][1] <- "Guild"



## FungalTraits ####

# download traits metadata
traits_meta <- read_csv("https://github.com/traitecoevo/fungaltraits/releases/download/v0.0.3/funtothefun.csv")

# download FungalTraits database
traits_db <- fungaltraits::fungal_traits()
names(traits_db$species)
# match taxa at genus level
genera <- fung@tax_table[,6] %>% str_remove("^g__")
species <- fung@tax_table[,7] %>% str_remove("^s__")
fungal_traits <- 
  data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(traits_db,by=c("species","Genus"),multiple='all')

# need to condense/remove multiple matches
fungal_traits %>% 
  dplyr::filter(species != "NA_NA")

# remove traits not associated with biochem functional potential
traits_to_ignore <- c(
  "redChannel_mean","redChannel_sd","RNAHelicase_count","RNApolymerase_count","spore_length",
  "spore_size","spore_width","sporocarp_chitin","sporocarp_N","sporocarp_protein","sporocarp_resp",           
  "taxonomic_level_fg","tissue_c","tissue_cn","tissue_cp","tissue_n","tissue_np","tissue_p","total_genes",
  "trehalase_count","latitude","map","greenChannel_mean","greenChannel_sd","heatShockProtein_count",
  "extension_rate","fruiting_body_size","mat","longitude","melanin_content","melanin_count",
  "coldShockProtein_count","dsDNA","blueChannel_mean","blueChannel_sd","ifungorum_number",
  "sterol_type","studyName","substrate","trait_fg","trophic_mode_fg",'notes_fg',"source_funguild_fg",
  "growth_form_fg","guild_fg","higher_clade","culture_media","culture_notes","elevation","em_expl",
  "em_text","colour_mean","confidence_fg","ascoma_development","ascoma_type","ascus_dehiscence",
  "uuid","obj_id","speciesMatched"
)

# group by species; summarize to find mean values with na.omit=TRUE
summarized_traits <- 
  fungal_traits %>% 
  dplyr::select(-all_of(traits_to_ignore)) %>% 
  dplyr::group_by(species) %>% 
  summarize(across(where(is.numeric),function(x){mean(x,na.rm=TRUE)}))

names(summarized_traits)

# join traits with tax_table species 

traits <- 
data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(summarized_traits,by=c("species"))

# ANALYSES ####

## Traits and plant health ####

  # H1: Seedling performance will be enhanced when root microbiomes have 
      # a greater proportion of microbial mutualists (e.g., ectomycorrhizal fungi) 
      # relative to saprotrophs, endophytes, or pathogens

 
# make relative abundance version of phyloseq object
fung_ra <- transform_sample_counts(fung,function(x){x/sum(x)})


# identify "mutualist" taxa
unname(fung_ra@tax_table[,1]) %>% unique

# just using "mycorrhizal" as the keyword...
mutualist_guilds <- 
  grep("[M,m]ycorrhizal",(fung_ra@tax_table[,1]),value = TRUE) %>% 
  unique()

# identify "saprotrophs"
saprotroph_guilds <- 
  grep("[S,s]aprotroph",(fung_ra@tax_table[,1]),value = TRUE) %>% 
  grep(pattern="[M,m]ycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()

# identify "pathogens"
pathogen_guilds <- 
  grep("[P,p]athogen|[P,p]arasite",(fung_ra@tax_table[,1]),value = TRUE) %>% 
  grep(pattern="[M,m]ycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()


# subset taxa to only mutualists; get row sums; this will be proportion of mutualists
# in each sample

mutualist_proportions <- 
fung_ra %>% 
  subset_taxa(Guild %in% mutualist_guilds) %>% 
  sample_sums()

# build data frame for modeling
mutualism_df <- 
  microbiome::meta(fung_ra) %>% 
  mutate(proportion_mutualist = mutualist_proportions)

# subset taxa to only saprotrophs; get row sums; this will be proportion of mutualists
# in each sample

saprotroph_proportions <- 
  fung_ra %>% 
  subset_taxa(Guild %in% saprotroph_guilds) %>% 
  sample_sums()

# build data frame for modeling
saprotroph_df <- 
  microbiome::meta(fung_ra) %>% 
  mutate(proportion_saprotroph = saprotroph_proportions)

# subset taxa to only mutualists; get row sums; this will be proportion of mutualists
# in each sample

pathogen_proportions <- 
  fung_ra %>% 
  subset_taxa(Guild %in% pathogen_guilds) %>% 
  sample_sums()

# build data frame for modeling
pathogen_df <- 
  microbiome::meta(fung_ra) %>% 
  mutate(proportion_pathogen = pathogen_proportions)

# join together all 3 major guilds
guild_df <- 
mutualism_df %>% 
  full_join(saprotroph_df) %>% 
  full_join(pathogen_df)

# plot (just GrandFir)
(
grandfir_mutualist_plot <- 
guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_mutualist,y=value,color=drought)) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of mutualist fungi",y="Scaled/Centered Value",color="Drought") +
  scale_color_manual(values = pal.discrete[c(2,7)])
)
saveRDS(grandfir_mutualist_plot,"./Output/figs/ITS_Mutualist_Plot_grandfir.RDS")

grandfir_mutualist_plot2 <- 
  guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_mutualist,y=value,color=ordered(fire_freq,levels=c("0","1","3")))) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of mutualist fungi",y="Scaled/Centered Value",color="Fire frequency") +
  scale_color_manual(values = fire_colors)
saveRDS(grandfir_mutualist_plot2,"./Output/figs/ITS_Mutualist_Plot_grandfir_by_fire.RDS")

(
  snowbrush_mutualist_plot <- 
    guild_df %>% 
    dplyr::filter(species == "Snowbrush") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                   "leaf_length","height","shoot_dm","final_root_dm"),
                 names_to="indicator") %>% 
    mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
    ggplot(aes(x=proportion_mutualist,y=value,color=drought)) +
    geom_point(alpha=.5) +
    geom_smooth(method='lm',se=FALSE) +
    facet_wrap(~indicator,scales = 'free') +
    theme_minimal() +
    theme(strip.text = element_text(face="bold",size=12)) +
    labs(x="Proportion of mutualist fungi",y="Scaled/Centered Value",color="Drought") +
    scale_color_manual(values = pal.discrete[c(2,7)])
)
saveRDS(snowbrush_mutualist_plot,"./Output/figs/ITS_Mutualist_Plot_snowbrush.RDS")

snowbrush_mutualist_plot2 <- 
  guild_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
  ggplot(aes(x=proportion_mutualist,y=value,color=ordered(fire_freq,levels=c("0","1","3")))) +
  geom_point(alpha=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of mutualist fungi",y="Scaled/Centered Value",color="Fire frequency") +
  scale_color_manual(values = fire_colors)
saveRDS(snowbrush_mutualist_plot2,"./Output/figs/ITS_Mutualist_Plot_snowbrush_by_fire.RDS")


(
  grandfir_saprotroph_plot <- 
    guild_df %>% 
    dplyr::filter(species == "GrandFir") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                   "leaf_length","height","shoot_dm","final_root_dm"),
                 names_to="indicator") %>% 
    mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
    ggplot(aes(x=proportion_saprotroph,y=value,color=drought)) +
    geom_point(alpha=.5) +
    geom_smooth(method='lm',se=FALSE) +
    facet_wrap(~indicator,scales = 'free') +
    theme_minimal() +
    theme(strip.text = element_text(face="bold",size=12)) +
    labs(x="Proportion of saprotrophic fungi",y="Scaled/Centered Value",color="Drought") +
    scale_color_manual(values = pal.discrete[c(2,7)])
)
saveRDS(grandfir_saprotroph_plot,"./Output/figs/ITS_Saprotroph_Plot_grandfir.RDS")

(
  snowbrush_saprotroph_plot <- 
    guild_df %>% 
    dplyr::filter(species == "Snowbrush") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                   "leaf_length","height","shoot_dm","final_root_dm"),
                 names_to="indicator") %>% 
    mutate(indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence()) %>% 
    ggplot(aes(x=proportion_saprotroph,y=value,color=drought)) +
    geom_point(alpha=.5) +
    geom_smooth(method='lm',se=FALSE) +
    facet_wrap(~indicator,scales = 'free') +
    theme_minimal() +
    theme(strip.text = element_text(face="bold",size=12)) +
    labs(x="Proportion of saprotrophic fungi",y="Scaled/Centered Value",color="Drought") +
    scale_color_manual(values = pal.discrete[c(2,7)])
)
saveRDS(snowbrush_saprotroph_plot,"./Output/figs/ITS_Saprotroph_Plot_snowbrush.RDS")



(
  grandfir_pathogen_plot <- 
    guild_df %>% 
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
    labs(x="Proportion of pathogenic fungi",y="Scaled/Centered Value",color="Drought") +
    scale_color_manual(values = pal.discrete[c(2,7)])
)
saveRDS(grandfir_pathogen_plot,"./Output/figs/ITS_Pathogen_Plot_grandfir.RDS")

grandfir_pathogen_plot2 <- 
  guild_df %>% 
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
  labs(x="Proportion of pathogenic fungi",y="Scaled/Centered Value",color="Fire frequency") +
  scale_color_manual(values = fire_colors)
saveRDS(grandfir_pathogen_plot2,"./Output/figs/ITS_Pathogen_Plot_grandfir_by_fire.RDS")


(
  snowbrush_pathogen_plot <- 
    guild_df %>% 
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
    labs(x="Proportion of pathogenic fungi",y="Scaled/Centered Value",color="Drought") +
    scale_color_manual(values = pal.discrete[c(2,7)])
)
saveRDS(snowbrush_pathogen_plot,"./Output/figs/ITS_Pathogen_Plot_snowbrush.RDS")

snowbrush_pathogen_plot2 <- 
  guild_df %>% 
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
  labs(x="Proportion of pathogenic fungi",y="Scaled/Centered Value",color="Fire frequency") +
  scale_color_manual(values = fire_colors)

saveRDS(snowbrush_pathogen_plot2,"./Output/figs/ITS_Pathogen_Plot_snowbrush_by_fire.RDS")

### Modeling ####

# model: plant health ~ mutualism_% * drought * fire_freq + (1|block)
mutualism_glm_grandfir <- 
  guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
      formula=value ~ proportion_mutualist * drought * inoculum_site + (1|block))
summary(mutualism_glm_grandfir)
saveRDS(mutualism_glm_grandfir,"./Output/figs/ITS_Mutualist_Model_GrandFir.RDS")

mutualism_glm_snowbrush <- 
guild_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_mutualist * drought * inoculum_site + (1|block))
summary(mutualism_glm_snowbrush)
saveRDS(mutualism_glm_snowbrush,"./Output/ITS_Mutualist_Model_Snowbrush.RDS")

grandfir_pathogen_glm <- 
  guild_df %>% 
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
saveRDS(grandfir_pathogen_glm,"./Output/ITS_Pathogen_Model_GrandFir.RDS")

snowbrush_pathogen_glm <- 
  guild_df %>% 
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
saveRDS(snowbrush_pathogen_glm,"./Output/ITS_Pathogen_Model_Snowbrush.RDS")

# pull together all model outputs from guild predictions on plant traits
full_guild_model_df <- 
clean_model_df(snowbrush_pathogen_glm) %>% mutate(species="Ceanothus") %>% 
  full_join(clean_model_df(grandfir_pathogen_glm) %>% mutate(species="Abies")) %>% 
  full_join(
    clean_model_df(mutualism_glm_snowbrush) %>% mutate(species="Ceanothus") %>% 
      full_join(clean_model_df(mutualism_glm_grandfir) %>% mutate(species="Abies"))
  ) %>% 
  dplyr::filter(effect=="fixed") %>% 
  select(-group)

saveRDS(full_guild_model_df,"./Output/ITS_Guild_Model_Output.RDS")



# EXPORT PHYSEQ W GUILDS ####

fung@tax_table[,'Guild'] <- 
fung@tax_table %>% 
  as.data.frame() %>% 
  mutate(Guild = case_when(Guild %in% pathogen_guilds ~ "Pathogen",
                           Guild %in% mutualist_guilds ~ "Mutualist")) %>% 
  pluck("Guild")

saveRDS(fung,"./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")

## Functional diversity and plant health ###
# I don't trust this...
# H3: Functional diversity: Seedling performance will be enhanced in root microbiomes 
# comprised of taxa that represent high functional diversity
# 
# # calculate "functional diversity" into a single number for each taxon
# x <- 
# traits %>% 
#   dplyr::select(-c(Genus,species)) %>% 
#   mutate(across(everything(),abs)) # transform negative values
# 
# # convert NaN to 0 for diversity estimates
# for(i in names(x)){
#   pluck(x,i)[is.nan(pluck(x,i))] <- 0
# }
# # functional diversity, measured as number of annotated functional potentials
# # This is pretty sparse due to database lacking good info on lots of taxa
# functional_div1 <- 
#   x %>% 
#   as.matrix() %>% 
#   vegan::specnumber()
# 
# # Within each sample...
# # Each taxon relative abundance should be multiplied by that taxon's functional div #
# 
# scaled_func_div <- c()
# for(sampleid in sample_names(fung_ra)){
#   
#   single_sample <- # pull each sample, one at a time   
#     fung_ra %>% 
#     subset_samples(sample_names(fung_ra) == sampleid) 
#   single_sample <- # find present taxa
#     single_sample %>% 
#     subset_taxa(taxa_sums(single_sample) > 0)
#   present_taxa <- which(taxa_names(fung_ra) %in% taxa_names(single_sample)) # use their names
#   scaled_func_div[sampleid] <- sum(functional_div1[present_taxa] * taxa_sums(single_sample))
#   # ^ sumof: relative abundance * functional diversity score for that taxon
# }  
# 
# 
# # add to sample data frame
# functional_df <- 
#   microbiome::meta(fung_ra) %>% 
#   mutate(scaled_func_div=scaled_func_div)
# 
# # plot (just GrandFir)
# (
#   functional_plot <- 
#     functional_df %>% 
#     dplyr::filter(species == "GrandFir") %>% 
#     mutate(across(c("wilting_scale","bud_number","leaf_number",
#                     "leaf_length","height","shoot_dm","final_root_dm"),
#                   scale)) %>% # scale/center all indicators
#     pivot_longer(c("wilting_scale","bud_number","leaf_number",
#                    "leaf_length","height","shoot_dm","final_root_dm"),
#                  names_to="indicator") %>% 
#     ggplot(aes(x=scaled_func_div,y=value)) +
#     geom_point() +
#     geom_smooth(method='lm') +
#     facet_wrap(~indicator,scales = 'free') +
#     theme_minimal() +
#     theme(strip.text = element_text(face="bold",size=12)) +
#     labs(x="Scaled functional diversity",y="Scaled/Centered Value")
# )
# saveRDS(functional_plot,"./Output/ITS_Functional_Plot.RDS")
# 
# 
# ### Modeling ###
# # model: plant health ~ mutualism_% + block + (1|block)
# 
# functional_glm <- 
#   functional_df %>% 
#   dplyr::filter(species == "GrandFir") %>% 
#   mutate(across(c("wilting_scale","bud_number","leaf_number",
#                   "leaf_length","height","shoot_dm","final_root_dm"),
#                 scale)) %>% # scale/center all indicators
#   pivot_longer(c("wilting_scale","bud_number","leaf_number",
#                  "leaf_length","height","shoot_dm","final_root_dm"),
#                names_to="indicator") %>% 
#   glm(data=.,
#       formula=value ~ scaled_func_div + scaled_func_div:indicator)
# summary(functional_glm)
# saveRDS(functional_glm,"./Output/ITS_Functional_Model.RDS")
# 
