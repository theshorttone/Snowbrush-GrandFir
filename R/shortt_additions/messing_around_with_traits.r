
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

fung_ra <- transform_sample_counts(fung,function(x){x/sum(x)})
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
  

# adding fungal traits to tax table --------------------------------------


#Build a per-ASV lookup table from tax_table (rows = taxa_names)
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

#Join traits onto each ASV via species_key 
tax_with_traits <- tax_lookup %>%
  left_join(traits_df, by = c("species_key" = "species"))

# taxa_df <- tax_with_traits %>%
#   select(taxon_id, everything(), -species_key) %>%  
#   tibble::column_to_rownames("taxon_id")


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
fung_traits_ra <- transform_sample_counts(fung_ra_traits,function(x){x/sum(x)})

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
  grep("Plant [P,p]athogen|Fungal Parasite",(fung_ra_traits@tax_table[,8]),value = TRUE) %>% 
  grep(pattern="Ectomycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()

  

mutualist_proportions <- 
  fung_ra %>% 
  subset_taxa(fung_ra_traits@tax_table[,8] %in% mutualist_guilds) %>% 
  sample_sums()

mutualism_df <- 
  microbiome::meta(fung_ra_traits) %>% 
  mutate(proportion_mutualist = mutualist_proportions)

saprotroph_proportions <- 
  fung_ra %>% 
  subset_taxa(fung_ra_traits@tax_table[,8] %in% saprotroph_guilds) %>% 
  sample_sums()

saprotroph_df <- 
  microbiome::meta(fung_ra_traits) %>% 
  mutate(proportion_saprotroph = saprotroph_proportions)

pathogen_proportions <- 
  fung_ra %>% 
  subset_taxa(fung_ra_traits@tax_table[,8] %in% pathogen_guilds) %>% 
  sample_sums()

pathogen_df <- 
  microbiome::meta(fung_ra_traits) %>% 
  mutate(proportion_pathogen = pathogen_proportions)

# Geoff code below unedited so far ----------------------------------------



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
saveRDS(grandfir_mutualist_plot,"R/shortt_additions/figures/test/ITSgrandfir_mutualist_plot.RDS")

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
saveRDS(grandfir_mutualist_plot2,"R/shortt_additions/figures/test/ITS_Mutualist_Plot_grandfir_by_fire.RDS")

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
saveRDS(snowbrush_mutualist_plot,"R/shortt_additions/figures/test/ITS_Mutualist_Plot_snowbrush.RDS")

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
snowbrush_mutualist_plot2
saveRDS(snowbrush_mutualist_plot2,"R/shortt_additions/figures/test/ITS_Mutualist_Plot_snowbrush_by_fire.RDS")


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
saveRDS(grandfir_saprotroph_plot,"R/shortt_additions/figures/test/ITS_Saprotroph_Plot_grandfir.RDS")

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
grandfir_pathogen_plot2
saveRDS(grandfir_pathogen_plot2,"R/shortt_additions/figures/ITS_Pathogen_Plot_grandfir_by_fire.RDS")


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
snowbrush_pathogen_plot2
saveRDS(snowbrush_pathogen_plot2,"R/shortt_additions/figures/ITS_Pathogen_Plot_snowbrush_by_fire.RDS")

### Modeling ####

# model: plant health ~ mutualism_% * drought + (1|block)
mutualism_glm_grandfir <- 
  guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  dplyr::select(all_of(c("leaf_number","shoot_dm","final_root_dm","proportion_mutualist","drought","block"))) %>% 
  mutate(across(c("leaf_number","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_mutualist * drought + (1|block))
summary(mutualism_glm_grandfir)
saveRDS(mutualism_glm_grandfir,"./Output/figs/ITS_Mutualist_Model_GrandFir.RDS")

mutualism_glm_snowbrush <- 
  guild_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  dplyr::select(all_of(c("leaf_number","shoot_dm","final_root_dm","proportion_mutualist","drought","block"))) %>% 
  mutate(across(c("leaf_number","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_mutualist * drought + (1|block))
summary(mutualism_glm_snowbrush)
saveRDS(mutualism_glm_snowbrush,"./Output/ITS_Mutualist_Model_Snowbrush.RDS")

grandfir_pathogen_glm <- 
  guild_df %>% 
  dplyr::filter(species == "GrandFir") %>% 
  dplyr::select(all_of(c("leaf_number","shoot_dm","final_root_dm","proportion_pathogen","drought","block"))) %>% 
  mutate(across(c("leaf_number","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_pathogen * drought + (1|block))
summary(grandfir_pathogen_glm)
saveRDS(grandfir_pathogen_glm,"./Output/ITS_Pathogen_Model_GrandFir.RDS")

snowbrush_pathogen_glm <- 
  guild_df %>% 
  dplyr::filter(species == "Snowbrush") %>% 
  dplyr::select(all_of(c("leaf_number","shoot_dm","final_root_dm","proportion_pathogen","drought","block"))) %>% 
  mutate(across(c("leaf_number","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("leaf_number","shoot_dm","final_root_dm"),
               names_to="indicator") %>% 
  lmer(data=.,
       formula=value ~ proportion_pathogen * drought + (1|block))
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
full_guild_model_df
saveRDS(full_guild_model_df,"R/shortt_additions/figures/ITS_Guild_Model_Output.RDS")


# EXPORT PHYSEQ W GUILDS ####

tt <- as.data.frame(tax_table(fung_ra_traits))

# 2. Recode Guild column
tt$guild_fg<- case_when(
  tt$guild_fg %in% pathogen_guilds  ~ "Pathogen",
  tt$guild_fg %in% mutualist_guilds ~ "Mutualist",
  TRUE                          ~ tt$guild_fg   # keep everything else unchanged
)

# 3. Write it back into the phyloseq object
tax_table(fung_ra_traits) <- tax_table(as.matrix(tt))

# 4. Save the modified object (THIS is the one you just edited)
saveRDS(
  fung_ra_traits,
  "./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS"
)

saveRDS(fung,"./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")

table(tax_table(fung_ra_traits)[,"guild_fg"])




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
  mutate(Moisture = case_when(drought == "D" ~ "Low",
                              drought == "ND" ~ "High")) %>% 
  ggplot(aes(x=proportion_mutualist,y=growth_response,color=Moisture)) +
  geom_point(alpha=.5,size=3) +
  geom_smooth(method='lm',se=FALSE,color='black') +
  annotate('text',
           x = .85,
           y = 3.5,
           label = "P = 0.044",
           fontface=2) +
  labs(x='Proportion of putative fungal mutualists',y="Plant growth response",title="Grand fir",color="Moisture") +
  theme_bw() +
  theme(axis.title = element_text(face='bold',size=16),plot.title = element_text(face='bold',size=12,hjust=.5)) +
  scale_color_manual(values = drought_colors)
saveRDS(gf_mutualist_plot,"./Output/figs/ITS_gf_mutualist_plot.RDS")

sb_mutualist_plot <- 
  sb %>% 
  ggplot(aes(x=proportion_mutualist,y=growth_response)) +
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
  dplyr::mutate(
    proportion_pathogen = as.numeric(proportion_pathogen),
    prop_pathogen_log10 = log10(proportion_pathogen + 1e-6)
  ) %>% 
  ggplot(aes(x = prop_pathogen_log10, y = growth_response)) +
  geom_point(alpha = .5, size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_continuous(
    breaks = log10(c(1e-6, 1e-4, 1e-2, 1e-1, 1)),
    labels = c("0", "0.0001", "0.01", "0.1", "1")
  ) +
  labs(x = "Proportion pathogen (log scale)", y = "", title = "Grand fir") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 12, hjust = .5)
  )
gf_pathogen_plot
saveRDS(gf_pathogen_plot,"R/shortt_additions/figures/test/gf_pathogen_logscale.RDS")


sb_pathogen_plot <- 
  sb %>% 
  dplyr::mutate(
    proportion_pathogen = as.numeric(proportion_pathogen),
    prop_pathogen_log10 = log10(proportion_pathogen + 1e-6)
  ) %>% 
  ggplot(aes(x = prop_pathogen_log10, y = growth_response)) +
  geom_point(alpha = .5, size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_continuous(
    breaks = log10(c(1e-6, 1e-4, 1e-2, 1e-1, 1)),
    labels = c("0", "0.0001", "0.01", "0.1", "1")
  ) +
  labs(x = "Proportion of pathogen ASVs (log scale)", y = "", title = "Snowbrush") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 14),
    legend.text  = element_text(face = "bold", size = 12),
    plot.title   = element_text(face = "bold", size = 12, hjust = .5)
  )
sb_pathogen_plot
saveRDS(sb_pathogen_plot,"R/shortt_additions/figures/test/sb_pathogen_logscale.RDS")



# prolly delete stuff below  ----------------------------------------------


# build 3-panel plot for manuscript
gf_fungal_mutualist_plot <- readRDS("./Output/figs/ITS_gf_mutualist_plot.RDS") + labs(x="Proportion of putative\nbacterial mutualists") + theme(legend.title = element_text(face='bold',size=14),
                                                                                                                                                legend.text = element_text(face='bold',size=12),
                                                                                                                                                legend.position = 'none', plot.title = element_text(face='bold',size=22), axis.title = element_text(face='bold',size=20)) + geom_point(size=6) + geom_smooth(method = 'lm',color='black')
gf_bacterial_mutualist_plot <- readRDS("./Output/figs/16S_gf_mutualist_plot.RDS") + labs(x='Proportion of putative\nbacterial mutualists',color="Moisture") + theme(legend.position = 'none', plot.title = element_text(face='bold',size=22), axis.title = element_text(face='bold',size=20)) + geom_point(size=6) + geom_smooth(method='lm',color='black')
sb_bacterial_pathogen_plot <- readRDS("./Output/figs/16S_sb_pathogen_plot.RDS") + labs(x='Proportion of putative\nbacterial pathogens',color="Moisture") + theme(plot.title = element_text(face='bold',size=22), axis.title = element_text(face='bold',size=20)) + geom_point(size=6) + geom_smooth(method='lm')

pw <- 
  (gf_fungal_mutualist_plot + gf_bacterial_mutualist_plot + sb_bacterial_pathogen_plot) &
  ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 5.5, 0))
guild_plots <- 
  wrap_elements(pw) +
  labs(tag = "Scaled plant growth measure\n") +
  theme(
    plot.tag = element_text(angle = 90,face='bold',size=16),
    plot.tag.position = "left"
  ) +
  plot_layout(guides = 'collect',axis_title="collect") +
  plot_annotation(title = "Fungal ASVs")
guild_plots
saveRDS(guild_plots,"./Output/figs/ALL_guild_plots.RDS")
ggsave("./Output/figs/ALL_guild_plots.tiff",device = 'tiff',width = 24,height = 8,dpi=400)
ggsave("./Output/figs/ALL_guild_plots.png",device = 'png',width = 24,height = 8,dpi=400)

# saveRDS(guild_plots,"./Output/figs/ITS_Guild_Plots.RDS")
bact_guild_plots <- readRDS("./Output/figs/16S_Guild_Plots.RDS") + plot_annotation(title = "Bacterial ASVs")
saveRDS(bact_guild_plots,"./Output/figs/16S_Guild_Plots.RDS")

guild_plots
ggsave("./Output/figs/ITS_Guild_Plots.tiff",device = 'tiff',height = 8,width = 8,dpi=350)

bact_guild_plots
ggsave("./Output/figs/16S_Guild_Plots.tiff",device = 'tiff',height = 8,width = 8,dpi=350)
