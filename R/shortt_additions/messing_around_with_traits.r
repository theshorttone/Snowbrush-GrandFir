
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
    traits_just_genus %>% select(Genus_key, n_trait_rows, all_of(cols)),
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

traits_df <- traits_one %>%
  select(-ends_with(".genus")) %>% 
  select(species, genus, guild_source, n_trait_rows, all_of(traits_to_coalesce))
  

# adding fungal traits to taxa table --------------------------------------


# 1) Build a per-ASV lookup table from tax_table (rows = taxa_names)
tax_lookup <- tax_table(fung_ra) %>%
  as("matrix") %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("taxon_id") %>%
  mutate(
    genus_raw   = str_remove(Genus, "^g__") %>% str_trim(),
    species_raw = str_remove(Species, "^s__") %>% str_trim(),
    species_key = paste(genus_raw, coalesce(species_raw, "NA"), sep = "_") %>%
      str_to_lower()
  )

# 2) Join traits onto each ASV via species_key (rename to match your column name in traits_one)
#    This assumes traits_one has a column named "species" that is the species_key.
tax_with_traits <- tax_lookup %>%
  left_join(traits_df, by = c("species_key" = "species"))

# 3) Write into taxa_data (use taxon_id as rownames)
taxa_df <- tax_with_traits %>%
  select(taxon_id, everything(), -species_key) %>%  # keep species_key if you want
  tibble::column_to_rownames("taxon_id")

taxa_data(fung_ra) <- taxa_data(taxa_df)



####UNTESTED BELOW# join traits with tax_table species 

# phyloseq with traits ----------------------------------------------------

tt <- as(tax_table(fung_ra), "matrix") %>%
  as.data.frame(stringsAsFactors = FALSE)

taxon_id <- rownames(tt)

tax_key <- tibble(
  taxon_id = taxon_id,
  genus_raw   = str_remove(tt$Genus, "^g__") %>% str_trim(),
  species_raw = str_remove(tt$Species, "^s__") %>% str_trim()
) %>%
  mutate(
    species_clean = case_when(
      is.na(species_raw) | species_raw == "" ~ NA_character_,
      TRUE ~ species_raw
    ),
    species_key = paste(genus_raw, coalesce(species_clean, "na"), sep = "_") %>% str_to_lower()
  )

trait_cols <- c("guild_fg", "guild_source", "trophic_mode_fg",
                "notes_fg", "source_funguild_fg", "culture_media")

tax_traits <- tax_key %>%
  left_join(traits_df %>% select(species, any_of(trait_cols)),
            by = c("species_key" = "species"))

# 3) append traits into tax_table (as extra columns)
tt2 <- tt
for (nm in trait_cols) tt2[[nm]] <- tax_traits[[nm]]
rownames(tt2) <- tax_traits$taxon_id

fung_ra_traits <- fung_ra   # copy the phyloseq object

# ...build tt2 as you already did...

tax_table(fung_ra_traits) <- tax_table(as.matrix(tt2))


a <- as.data.frame(tax_table(fung_ra_traits))


# make relative abundance version of phyloseq object
fung_ra <- transform_sample_counts(fung,function(x){x/sum(x)})

unique(as.character(tax_table(fung_ra_traits)[, "guild_fg"]))

# just using "mycorrhizal" as the keyword...
mutualist_guilds <- 
  grep("Ectomycorrhizal",(fung_ra_traits@tax_table[,8]),value = TRUE) %>% 
  unique()

# identify "saprotrophs"
saprotroph_guilds <- 
  grep("[S,s]aprotroph",(fung_ra_traits@tax_table[,8]),value = TRUE) %>% 
  grep(pattern="Ectomycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()

# identify "pathogens"
pathogen_guilds <- 
  grep("Plant [P,p]athogen|Fungal Parasite",(fung_ra@tax_table[,1]),value = TRUE) %>% 
  grep(pattern="Ectomycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()

  

mutualist_proportions <- 
  fung_ra %>% 
  subset_taxa(Guild %in% mutualist_guilds) %>% 
  sample_sums()


