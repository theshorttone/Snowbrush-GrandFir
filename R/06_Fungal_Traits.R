# Fungal traits/guilds

# SETUP ####

# uncommon specialty package installation
           # devtools::install_github("brendanf/FUNGuildR")
           # devtools::install_github("ropenscilabs/datastorr")
           # devtools::install_github("traitecoevo/fungaltraits")

# load packages
library(tidyverse)
library(phyloseq)
library(FUNGuildR)
library(fungaltraits)
library(vegan)
library(microbiome)

# load data
fung <- readRDS("./Output/ITS_clean_phyloseq_object.RDS")
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
# guild_db <- readRDS("./Taxonomy/Funguild_Database.RDS")

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
mutualist_plot <- 
guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  ggplot(aes(x=proportion_mutualist,y=value)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_wrap(~indicator,scales = 'free') +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=12)) +
  labs(x="Proportion of mutualist fungi",y="Scaled/Centered Value")
)
saveRDS(mutualist_plot,"./Output/ITS_Mutualist_Plot.RDS")

(
  saprotroph_plot <- 
    guild_df %>% 
    dplyr::filter(species == "GrandFir") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                   "leaf_length","height","shoot_dm","final_root_dm"),
                 names_to="indicator") %>% 
    ggplot(aes(x=proportion_saprotroph,y=value)) +
    geom_point() +
    geom_smooth(method='lm') +
    facet_wrap(~indicator,scales = 'free') +
    theme_minimal() +
    theme(strip.text = element_text(face="bold",size=12)) +
    labs(x="Proportion of saprotrophic fungi",y="Scaled/Centered Value")
)
saveRDS(saprotroph_plot,"./Output/ITS_Saprotroph_Plot.RDS")

(
  pathogen_plot <- 
    guild_df %>% 
    dplyr::filter(species == "GrandFir") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                   "leaf_length","height","shoot_dm","final_root_dm"),
                 names_to="indicator") %>% 
    ggplot(aes(x=proportion_pathogen,y=value)) +
    geom_point() +
    geom_smooth(method='lm') +
    facet_wrap(~indicator,scales = 'free') +
    theme_minimal() +
    theme(strip.text = element_text(face="bold",size=12)) +
    labs(x="Proportion of pathogenic fungi",y="Scaled/Centered Value")
)
saveRDS(pathogen_plot,"./Output/ITS_Pathogen_Plot.RDS")

### Modeling ####
# model: plant health ~ mutualism_% + block + (1|block)

mutualism_glm <- 
guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  glm(data=.,
      formula=value ~ proportion_mutualist + proportion_mutualist:indicator)
summary(mutualism_glm)
saveRDS(mutualism_glm,"./Output/ITS_Mutualist_Model.RDS")

saprotroph_glm <- 
  guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  glm(data=.,
      formula=value ~ proportion_saprotroph + proportion_saprotroph:indicator)
summary(saprotroph_glm)
saveRDS(saprotroph_glm,"./Output/ITS_Saprotroph_Model.RDS")

pathogen_glm <- 
  guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  glm(data=.,
      formula=value ~ proportion_pathogen + proportion_pathogen:indicator)
summary(pathogen_glm)
saveRDS(pathogen_glm,"./Output/ITS_Pathogen_Model.RDS")

## Functional diversity and plant health ####

# H3: Functional diversity: Seedling performance will be enhanced in root microbiomes 
# comprised of taxa that represent high functional diversity

# calculate "functional diversity" into a single number for each taxon
x <- 
traits %>% 
  dplyr::select(-c(Genus,species)) %>% 
  mutate(across(everything(),abs)) # transform negative values

# convert NaN to 0 for diversity estimates
for(i in names(x)){
  pluck(x,i)[is.nan(pluck(x,i))] <- 0
}
# functional diversity, measured as number of annotated functional potentials
# This is pretty sparse due to database lacking good info on lots of taxa
functional_div1 <- 
  x %>% 
  as.matrix() %>% 
  vegan::specnumber()

# Within each sample...
# Each taxon relative abundance should be multiplied by that taxon's functional div #

scaled_func_div <- c()
for(sampleid in sample_names(fung_ra)){
  
  single_sample <- # pull each sample, one at a time   
    fung_ra %>% 
    subset_samples(sample_names(fung_ra) == sampleid) 
  single_sample <- # find present taxa
    single_sample %>% 
    subset_taxa(taxa_sums(single_sample) > 0)
  present_taxa <- which(taxa_names(fung_ra) %in% taxa_names(single_sample)) # use their names
  scaled_func_div[sampleid] <- sum(functional_div1[present_taxa] * taxa_sums(single_sample))
  # ^ sumof: relative abundance * functional diversity score for that taxon
}  


# add to sample data frame
functional_df <- 
  microbiome::meta(fung_ra) %>% 
  mutate(scaled_func_div=scaled_func_div)

# plot (just GrandFir)
(
  functional_plot <- 
    functional_df %>% 
    dplyr::filter(species == "GrandFir") %>% 
    mutate(across(c("wilting_scale","bud_number","leaf_number",
                    "leaf_length","height","shoot_dm","final_root_dm"),
                  scale)) %>% # scale/center all indicators
    pivot_longer(c("wilting_scale","bud_number","leaf_number",
                   "leaf_length","height","shoot_dm","final_root_dm"),
                 names_to="indicator") %>% 
    ggplot(aes(x=scaled_func_div,y=value)) +
    geom_point() +
    geom_smooth(method='lm') +
    facet_wrap(~indicator,scales = 'free') +
    theme_minimal() +
    theme(strip.text = element_text(face="bold",size=12)) +
    labs(x="Scaled functional diversity",y="Scaled/Centered Value")
)
saveRDS(functional_plot,"./Output/ITS_Functional_Plot.RDS")


### Modeling ####
# model: plant health ~ mutualism_% + block + (1|block)

functional_glm <- 
  functional_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  glm(data=.,
      formula=value ~ scaled_func_div + scaled_func_div:indicator)
summary(functional_glm)
saveRDS(functional_glm,"./Output/ITS_Functional_Model.RDS")

