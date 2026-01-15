
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

# calculate proportions of major guilds ------------------------------------
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


# join together all 3 major guilds
guild_df <- 
  mutualism_df %>% 
  full_join(saprotroph_df) %>% 
  full_join(pathogen_df)

indicators <- c("wilting_scale","bud_number","leaf_number",
                "leaf_length","height","shoot_dm","final_root_dm")

guild_plot_df <- guild_df %>% 
  mutate(across(all_of(indicators), scale)) %>% 
  pivot_longer(all_of(indicators),
               names_to = "indicator",
               values_to = "value") %>% 
  mutate(
    indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence(),
    
    proportion_mutualist = as.numeric(proportion_mutualist),
    proportion_pathogen  = as.numeric(proportion_pathogen),
    proportion_saprotroph = as.numeric(proportion_saprotroph),
    
    prop_mut_log10  = log10(proportion_mutualist + 1e-6),
    prop_path_log10 = log10(proportion_pathogen  + 1e-6),
    prop_sap_log10  = log10(proportion_saprotroph + 1e-6),
    
    # nicer drought label if you want it
    Moisture = case_when(drought == "D"  ~ "Low",
                         drought == "ND" ~ "High",
                         TRUE ~ as.character(drought))
  )

#plotting function:
plot_guild_effects <- function(df,
                               species_sel = c("GrandFir","Snowbrush"),
                               guild = c("mutualist","pathogen","saprotroph"),
                               color_by = c("drought","fire"),
                               log_x = TRUE) {
  
  species_sel <- match.arg(species_sel)
  guild    <- match.arg(guild)
  color_by <- match.arg(color_by)
  
  # choose x column
  x_col <- dplyr::case_when(
    guild == "mutualist"  & log_x ~ "prop_mut_log10",
    guild == "mutualist"  & !log_x ~ "proportion_mutualist",
    guild == "pathogen"   & log_x ~ "prop_path_log10",
    guild == "pathogen"   & !log_x ~ "proportion_pathogen",
    guild == "saprotroph" & log_x ~ "prop_sap_log10",
    guild == "saprotroph" & !log_x ~ "proportion_saprotroph"
  )
  
  df2 <- df %>% dplyr::filter(.data$species == species_sel)
  
  if (color_by == "drought") {
    col_aes <- rlang::sym("Moisture")   # or "drought" for D/ND
    color_scale <- scale_color_manual(values = drought_colors)
    col_title <- "Moisture"
  } else {
    df2 <- df2 %>% dplyr::mutate(fire_freq = ordered(fire_freq, levels = c("0","1","3")))
    col_aes <- rlang::sym("fire_freq")
    color_scale <- scale_color_manual(values = fire_colors)
    col_title <- "Fire frequency"
  }
  
  p <- ggplot(df2, aes(x = .data[[x_col]], y = value, color = !!col_aes)) +
    geom_point(alpha = .5) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~indicator, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 12)) +
    labs(
      title = species_sel,
      x = paste0("Proportion of ", guild, " fungi", ifelse(log_x, " (log scale)", "")),
      y = "Scaled/Centered Value",
      color = col_title
    ) +
    color_scale
  
  # add log tick labels only when using prelogged x
  if (log_x) {
    p <- p + scale_x_continuous(
      breaks = log10(c(1e-6, 1e-4, 1e-2, 1e-1, 1)),
      labels = c("0", "0.0001", "0.01", "0.1", "1")
    )
  }
  
  p
}

p1 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="mutualist", color_by="drought", log_x=TRUE)
p2 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="pathogen", color_by="drought", log_x=TRUE)
p3 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="mutualist", color_by="drought", log_x=TRUE)
p4 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="pathogen", color_by="drought", log_x=TRUE)
p5 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="mutualist", color_by="fire", log_x=TRUE)
p6 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="pathogen", color_by="fire", log_x=TRUE)
p7 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="mutualist", color_by="fire", log_x=TRUE)
p8 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="pathogen", color_by="fire", log_x=TRUE)
p9 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="saprotroph", color_by="fire", log_x=TRUE)
p10 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="saprotroph", color_by="fire", log_x=TRUE)


# p1;ggsave("ITS_guild_all_growth/log_scale_GF_mutualist_drought", p1)
# p2;ggsave("ITS_guild_all_growth/log_scale_GF_pathogen_drought", p2)
# p3;ggsave("ITS_guild_all_growth/log_scale_SB_mutualist_drought", p3)
# p4;ggsave("ITS_guild_all_growth/log_scale_SB_pathogen_drought", p4)
# p5;ggsave("ITS_guild_all_growth/log_scale_GF_mutualist_fire", p5)
# p6;ggsave( "ITS_guild_all_growth/log_scale_GF_pathogen_fire", p6)
# p7;ggsave("ITS_guild_all_growth/log_scale_SB_mutualist_fire", p7)
# p8;ggsave("ITS_guild_all_growth/log_scale_SB_pathogen_fire", p8)
# p9;ggsave("ITS_guild_all_growth/log_scale_SB_saprotroph_fire", p9)
# p10;ggsave("ITS_guild_all_growth/log_scale_GF_saprotroph_fire", p10)



# Modeling ####
clean_model_df <- function(model) {
  broom.mixed::tidy(model, effects = "fixed")  # with lmerTest loaded, includes df + p.value
}

fit_one_lmer <- function(df, species_sel, indicator_sel,
                         predictor, interact = c("drought","fire_freq"),
                         log_x = FALSE, eps = 1e-6) {
  
  interact <- match.arg(interact)
  
  # choose x column name
  x_col <- dplyr::case_when(
    predictor == "mutualist"  & log_x ~ "prop_mut_log10",
    predictor == "mutualist"  & !log_x ~ "proportion_mutualist",
    predictor == "pathogen"   & log_x ~ "prop_path_log10",
    predictor == "pathogen"   & !log_x ~ "proportion_pathogen",
    predictor == "saprotroph" & log_x ~ "prop_sap_log10",
    predictor == "saprotroph" & !log_x ~ "proportion_saprotroph"
  )
  
  dat <- df %>%
    dplyr::filter(.data$species == species_sel) %>%
    dplyr::select(all_of(c(indicator_sel, x_col, interact, "block"))) %>%
    dplyr::mutate(
      value = scale01(.data[[indicator_sel]]),
      x = as.numeric(.data[[x_col]])
    )
  
  if (interact == "fire_freq") {
    dat <- dat %>% dplyr::mutate(fire_freq = factor(fire_freq, levels = c("0","1","3")))
    mod <- lmerTest::lmer(value ~ x * fire_freq + (1|block), data = dat)
  } else {
    dat <- dat %>% dplyr::mutate(drought = factor(drought))
    mod <- lmerTest::lmer(value ~ x * drought + (1|block), data = dat)
  }
  
  mod
}


prep_guild_df <- function(guild_df, eps = 1e-6) {
  
}

eps <-  1e-6
guild_df2 <-guild_df %>%
  dplyr::mutate(
    proportion_mutualist  = as.numeric(proportion_mutualist),
    proportion_pathogen   = as.numeric(proportion_pathogen),
    proportion_saprotroph = as.numeric(proportion_saprotroph),
    prop_mut_log10  = log10(proportion_mutualist  + eps),
    prop_path_log10 = log10(proportion_pathogen   + eps),
    prop_sap_log10  = log10(proportion_saprotroph + eps)
  )

# combined traits  --------------------------------------------------------



indicators <- c("leaf_number","shoot_dm","final_root_dm", "height")  # or your full set

make_results_df <- function(df, interact = c("drought","fire_freq")) {
  interact <- match.arg(interact)
  
  grid <- tidyr::expand_grid(
    species_sel = c("GrandFir","Snowbrush"),
    indicator   = indicators,
    predictor   = c("mutualist","pathogen","saprotroph"),
    log_x       = c(FALSE, TRUE)
  )
  
  res <- grid %>%
    dplyr::mutate(
      model = pmap(
        list(species_sel, indicator, predictor, log_x),
        ~ fit_one_lmer(df,
                       species_sel = ..1,
                       indicator_sel = ..2,
                       predictor = ..3,
                       interact = interact,
                       log_x = ..4)
      ),
      tidy = map(model, clean_model_df)
    ) %>%
    tidyr::unnest(tidy) %>%
    dplyr::mutate(
      interact = interact,
      scale = ifelse(log_x, "log10(x+1e-6)", "linear")
    ) %>%
    dplyr::filter(effect == "fixed")
  res
}

results_drought <- make_results_df(guild_df2, interact = "drought")
results_fire    <- make_results_df(guild_df2, interact = "fire_freq")

# combine if you want
results_all <- dplyr::bind_rows(results_drought, results_fire) %>% 
  filter(term != "(Intercept)") %>% 
  filter(p.value <= .05)

saveRDS(results_all, "R/shortt_additions/stats/p_hacked_ITS_results_all_individual_growth(including log).RDS")



# scaled fig --------------------------------------------------------------
growth_traits <- c("leaf_number","shoot_dm","final_root_dm", "height")

guild_df_comp <- guild_df %>%
  dplyr::mutate(
    proportion_mutualist  = as.numeric(proportion_mutualist),
    proportion_pathogen   = as.numeric(proportion_pathogen),
    proportion_saprotroph = as.numeric(proportion_saprotroph),
    
    prop_mut_log10  = log10(proportion_mutualist  + eps),
    prop_path_log10 = log10(proportion_pathogen   + eps),
    prop_sap_log10  = log10(proportion_saprotroph + eps)
  ) %>%
  dplyr::mutate(across(all_of(growth_traits), scale01)) %>%
  dplyr::mutate(
    growth_index = rowMeans(dplyr::across(all_of(growth_traits)), na.rm = TRUE)
  )



plot_guild_growth_index_log <- function(df,
                                        species_sel = c("GrandFir","Snowbrush"),
                                        guild = c("mutualist","pathogen","saprotroph"),
                                        color_by = c("none","drought","fire_freq")) {
  
  species_sel <- match.arg(species_sel)
  guild       <- match.arg(guild)
  color_by    <- match.arg(color_by)
  
  df2 <- df %>%
    dplyr::filter(.data$species == species_sel) %>%
    dplyr::mutate(
      x_log = dplyr::case_when(
        guild == "mutualist"  ~ .data$prop_mut_log10,
        guild == "pathogen"   ~ .data$prop_path_log10,
        guild == "saprotroph" ~ .data$prop_sap_log10
      ),
      Moisture = dplyr::case_when(
        drought == "D"  ~ "Low",
        drought == "ND" ~ "High",
        TRUE ~ as.character(drought)
      ),
      fire_freq = factor(fire_freq, levels = c("0","1","3"))
    )
  
  base <- ggplot2::ggplot(df2, ggplot2::aes(x = x_log, y = growth_index)) +
    ggplot2::geom_point(alpha = .5, size = 3) +
    ggplot2::geom_smooth(method = "lm", se = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = log10(c(1e-6, 1e-4, 1e-2, 1e-1, 1)),
      labels = c("0", "0.0001", "0.01", "0.1", "1")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 16),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = .5),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text  = ggplot2::element_text(face = "bold", size = 12)
    ) +
    ggplot2::labs(
      x = paste0("Proportion of ", guild, " fungi (log scale)"),
      y = "Combined scaled growth index",
      title = ifelse(species_sel == "GrandFir", "Grand fir", "Snowbrush")
    )
  
  if (color_by == "none") {
    return(base)
  }
  
  if (color_by == "drought") {
    return(
      base +
        ggplot2::aes(color = Moisture) +
        ggplot2::scale_color_manual(values = drought_colors) +
        ggplot2::labs(color = "Moisture")
    )
  }
  
  # color_by == "fire_freq"
  base +
    ggplot2::aes(color = fire_freq) +
    ggplot2::scale_color_manual(values = fire_colors) +
    ggplot2::labs(color = "Fire frequency")
}





# Grand fir, pathogens, color by drought
p1 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "pathogen", "fire_freq")

p1# Snowbrush, mutualists, color by fire
p2 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "mutualist", "fire_freq")

# No color
p3 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "saprotroph", "fire_freq")

p4 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "pathogen", "fire_freq")

p5 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "mutualist", "fire_freq")

p6 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "saprotroph", "fire_freq")

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave(combined_plot,
       filename = "R/shortt_additions/figures/ITS_guild_all_growth/log_scale_combined_growth_index_fire.pdf",
       width = 12,
       height = 8,
       units = "in",
       dpi = 300)

p1 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "pathogen", "drought")

p1# Snowbrush, mutualists, color by fire
p2 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "mutualist", "drought")

# No color
p3 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "saprotroph", "drought")

p4 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "pathogen", "drought")

p5 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "mutualist", "drought")

p6 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "saprotroph", "drought")

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6)
combined_plot
ggsave(combined_plot,
       filename = "R/shortt_additions/figures/ITS_guild_all_growth/log_scale_combined_growth_index_drought.pdf",
       width = 12,
       height = 8,
       units = "in",
       dpi = 300)


# Scaled all growth model ------------------------------------------------------------

#
fit_one_lmer_composite <- function(df, species_sel,
                                   predictor = c("mutualist","pathogen","saprotroph"),
                                   interact = c("drought","fire_freq"),
                                   log_x = FALSE) {
  
  predictor <- match.arg(predictor)
  interact  <- match.arg(interact)
  
  # choose the predictor column
  x_col <- dplyr::case_when(
    predictor == "mutualist"  & log_x ~ "prop_mut_log10",
    predictor == "mutualist"  & !log_x ~ "proportion_mutualist",
    predictor == "pathogen"   & log_x ~ "prop_path_log10",
    predictor == "pathogen"   & !log_x ~ "proportion_pathogen",
    predictor == "saprotroph" & log_x ~ "prop_sap_log10",
    predictor == "saprotroph" & !log_x ~ "proportion_saprotroph"
  )
  
  dat <- df %>%
    dplyr::filter(.data$species == species_sel) %>%
    dplyr::mutate(
      x = .data[[x_col]]
    )
  
  if (interact == "fire_freq") {
    dat <- dat %>% dplyr::mutate(fire_freq = factor(fire_freq, levels = c("0","1","3")))
    mod <- lmerTest::lmer(growth_index ~ x * fire_freq + (1|block), data = dat)
  } else {
    dat <- dat %>% dplyr::mutate(drought = factor(drought))
    mod <- lmerTest::lmer(growth_index ~ x * drought + (1|block), data = dat)
  }
  
  mod
}


make_results_df_composite <- function(df, interact = c("drought","fire_freq")) {
  interact <- match.arg(interact)
  
  grid <- tidyr::expand_grid(
    species_sel = c("GrandFir","Snowbrush"),
    predictor   = c("mutualist","pathogen","saprotroph"),
    log_x       = c(FALSE, TRUE)
  )
  
  res <- grid %>%
    mutate(
      model = purrr::pmap(
        list(species_sel, predictor, log_x),
        ~ fit_one_lmer_composite(df,
                                 species_sel = ..1,
                                 predictor   = ..2,
                                 interact    = interact,
                                 log_x       = ..3)
      ),
      tidy = purrr::map(model, clean_model_df)
    ) %>%
    tidyr::unnest(tidy) %>%
    mutate(
      interact = interact,
      scale = ifelse(log_x, "log10(x+1e-6)", "linear")
    )
  
  res
}

results_drought_comp <- make_results_df_composite(guild_df_comp, interact = "drought")
results_fire_comp    <- make_results_df_composite(guild_df_comp, interact = "fire_freq")
results_all_comp     <- bind_rows(results_drought_comp, results_fire_comp) %>% 
  filter(term != "(Intercept)") %>%
  filter(p.value <= .1)
saveRDS(results_all_comp, "R/shortt_additions/stats/p_hacked_ITS_results_combined_growth(including log).RDS")






# proportion mutualist bar graph ------------------------------------------


