library(tidyverse)
library(phyloseq)
library(zahntools)
library(vegan); packageVersion("vegan")
library(ggnewscale)

ra <- function(x){x/sum(x)}
# load data

# Data ####
inoc_16S <- readRDS("./Output/phyloseq_objects/16S_inoculum_samples_clean_phyloseq_object.RDS")
inoc_16S@sam_data$site <- inoc_16S@sam_data$other_frompreviouscolumn
inoc_16S_m <- inoc_16S %>% 
  merge_samples('site',fun = 'sum') 
inoc_16S_m <- inoc_16S_m %>% 
  subset_taxa(taxa_sums(inoc_16S_m) > 0) %>% 
  transform_sample_counts(ra)


inoc_ITS <- readRDS("./Output/phyloseq_objects/ITS_inoculum_samples_clean_phyloseq_object.RDS")
inoc_ITS@sam_data$site <- inoc_ITS@sam_data$other_frompreviouscolumn
inoc_ITS_m <- inoc_ITS %>% 
  merge_samples('site',fun = 'sum') 
inoc_ITS_m <- inoc_ITS_m %>% 
  subset_taxa(taxa_sums(inoc_ITS_m) > 0) %>% 
  transform_sample_counts(ra)

inoc_18S <- readRDS("./Output/phyloseq_objects/18S_inoculum_samples_clean_phyloseq_object.RDS")
inoc_18S@sam_data$site <- inoc_18S@sam_data$other_frompreviouscolumn
inoc_18S_m <- inoc_18S %>% 
  merge_samples('site',fun = 'sum') 
inoc_18S_m <- inoc_18S_m %>% 
  subset_taxa(taxa_sums(inoc_18S_m) > 0) %>% 
  transform_sample_counts(ra)

## FungalTraits ####

# download traits metadata
traits_meta <- read_csv("https://github.com/traitecoevo/fungaltraits/releases/download/v0.0.3/funtothefun.csv")

# download FungalTraits database
traits_db <- fungaltraits::fungal_traits()
head(traits_db$species)
# match taxa at species level
genera <- inoc_ITS@tax_table[,6] %>% str_remove("^g__")
species <- inoc_ITS@tax_table[,7] %>% str_remove("^s__")
fungal_traits <- 
  data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(traits_db,by=c("species","Genus"),multiple='all')%>% 
  dplyr::filter(species != "NA_NA")

# remove traits not associated with biochem functional potential
traits_to_ignore <- c(
  "redChannel_mean","redChannel_sd","RNAHelicase_count","RNApolymerase_count","spore_length",
  "spore_size","spore_width","sporocarp_chitin","sporocarp_N","sporocarp_protein","sporocarp_resp",           
  "taxonomic_level_fg","tissue_c","tissue_cn","tissue_cp","tissue_n","tissue_np","tissue_p","total_genes",
  "trehalase_count","latitude","map","greenChannel_mean","greenChannel_sd","heatShockProtein_count",
  "extension_rate","fruiting_body_size","mat","longitude","melanin_content","melanin_count",
  "coldShockProtein_count","dsDNA","blueChannel_mean","blueChannel_sd","ifungorum_number",
  "sterol_type","studyName","substrate","trait_fg","trophic_mode_fg",'notes_fg',"higher_clade","culture_media","culture_notes","elevation","em_expl",
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


# make relative abundance version of phyloseq object
fung_ra <- transform_sample_counts(inoc_ITS,function(x){x/sum(x)})


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

tax <- as.data.frame(tax_table(inoc_ITS))

tax$guild_fg         <- fungal_traits$guild_fg
tax$trophic_mode_fg  <- fungal_traits$trophic_mode_fg
tax$trait_fg         <- fungal_traits$trait_fg




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
