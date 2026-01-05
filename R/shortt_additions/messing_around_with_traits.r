
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(FUNGuildR); packageVersion("FUNGuildR")
library(fungaltraits); packageVersion("fungaltraits")
library(broom); packageVersion("broom")
library(lmerTest); packageVersion("lmerTest")
library(broom.mixed); packageVersion("broom.mixed")


source("./R/palettes.R")
source("./R/scale01.R")
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

# remove traits not associated with biochem functional potential
traits_to_ignore <- c(
  "redChannel_mean","redChannel_sd","RNAHelicase_count","RNApolymerase_count","spore_length",
  "spore_size","spore_width","sporocarp_chitin","sporocarp_N","sporocarp_protein","sporocarp_resp",           
  "tissue_c","tissue_cn","tissue_cp","tissue_n","tissue_np","tissue_p","total_genes",
  "trehalase_count","latitude","map","greenChannel_mean","greenChannel_sd","heatShockProtein_count",
  "extension_rate","fruiting_body_size","mat","longitude","melanin_content","melanin_count",
  "coldShockProtein_count","dsDNA","blueChannel_mean","blueChannel_sd","ifungorum_number",
  "sterol_type","studyName","higher_clade","elevation","em_expl","colour_mean","ascoma_development","ascoma_type","ascus_dehiscence",
  "uuid","obj_id", "Genus"
)

# functions ---------------------------------------------------------------

# load data
fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS")
fung@sam_data %>% row.names()

tax_df <- as.data.frame(tax_table(fung_ra))



traits_meta <- read_csv("https://github.com/traitecoevo/fungaltraits/releases/download/v0.0.3/funtothefun.csv")

# download FungalTraits database
traits_db <- fungaltraits::fungal_traits()
unique(traits_db$species)
# match taxa at genus level
genera <- fung@tax_table[,6] %>% str_remove("^g__")
species <- fung@tax_table[,7] %>% str_remove("^s__")

traits_db_norm <- traits_db %>%
  mutate(
    species = species %>%
      str_trim() %>%
      str_to_lower() %>%
      str_replace_all("\\s+", "_"),
    Genus = Genus %>%
      str_to_lower()
      )

#species level assignment db
fungal_traits_sp <- 
  data.frame(genus = genera) %>% 
  mutate(
    species = paste(genus, species, sep = "_") %>%
      str_to_lower(),
    genus = str_to_lower(genus)
    ) %>% 
    filter(species != "na_na") %>% 
  distinct(species, .keep_all = TRUE) %>% 
  left_join(traits_db_norm, by = "species", multiple = "all") %>%
  select(-all_of(traits_to_ignore))
         

#genus level assignment for those not matched at species 

fungal_traits_genus <- 
  data.frame(genus = genera) %>% 
  mutate(
    Genus_key = str_trim(genus) %>% str_to_lower()
  ) %>% 
  filter(!str_detect(Genus_key, "gen_incertae_sedis")) %>% 
  distinct(Genus_key, .keep_all = TRUE) %>% 
  left_join(traits_db_norm, by = c("Genus_key" = "Genus"), multiple = "all") %>% 
  select(-any_of(traits_to_ignore))

#function
collapse_chars <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(sort(x), collapse = "|")
}

traits_by_species <- fungal_traits_sp %>%
  group_by(species, genus) %>%
  summarise(
    n_trait_rows = n(),
    # numeric traits: mean
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    # character traits: list/union
    across(where(is.character), collapse_chars),
    .groups = "drop"
  )

traits_just_genus <- fungal_traits_genus %>%
  group_by(Genus_key) %>%
  summarise(
    n_trait_rows = n(),
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    across(where(is.character), collapse_chars),
    .groups = "drop"
  ) %>% filter(!is.na(speciesMatched))


traits_to_coalesce <- c(
  "guild_fg",
  "trophic_mode_fg",
  "culture_media",
  "notes_fg",
  "source_funguild_fg",
  "speciesMatched", 
  "confidence_fg"
)

cols <- traits_to_coalesce
cols_genus <- paste0(cols, ".genus")

traits_one <- traits_by_species %>%
  left_join(
    traits_just_genus %>% select(Genus_key, all_of(cols)),
    by = c("genus" = "Genus_key"),
    suffix = c("", ".genus")
  ) %>%
  mutate(
    guild_source = case_when(
      !is.na(guild_fg) ~ "species",
      is.na(guild_fg) & !is.na(guild_fg.genus) ~ "genus",
      TRUE ~ "none"
    )
  ) %>%
  mutate(
    across(all_of(cols), ~ na_if(.x, "")) # optional
  )

# do the coalesce using base R, then put back into the data frame
traits_one[cols] <- Map(
  coalesce,
  traits_one[cols],
  traits_one[cols_genus]
)

traits_one <- traits_one %>%
  select(-ends_with(".genus")) %>% 
  select(species, genus, all_of(traits_to_coalesce))






# join traits with tax_table species 

traits <- 
  data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(summarized_traits,by=c("species"))



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