library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(fungaltraits); packageVersion("fungaltraits")
library(broom); packageVersion("broom")
library(lmerTest); packageVersion("lmerTest")
library(broom.mixed); packageVersion("broom.mixed")
library(ggthemes)
library(patchwork)




source("./R/shortt_additions/funguild.R")
# source("./R/palettes.R")
# source("./R/scale01.R")
# drought_colors <- pal.discrete[c(2,5)]
# host_colors <- pal.discrete[c(7,10)] 

# EcM ---------------------------------------------------------------------

fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS")

site1.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site1.inoc.full.RDS")
site2.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site2.inoc.full.RDS")
site3.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site3.inoc.full.RDS")
site4.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site4.inoc.full.RDS")
site5.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site5.inoc.full.RDS")
site6.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site6.inoc.full.RDS")

inoc.full<- merge_phyloseq(site1.inoc.full,site2.inoc.full,site3.inoc.full,
                           site4.inoc.full,site5.inoc.full,site6.inoc.full)

#tax_df <- as.data.frame(tax_table(fung)) %>% rownames_to_column("taxon_id")
#tax_df <- as.data.frame(tax_table(fung_glom)) %>% rownames_to_column("taxon_id")



# functions ---------------------------------------------------------------

clean_func <- function (df){
  df_samples <- df %>%
    mutate(
      plot_nmds_color = case_when(
        !is.na(.data$inoculum) & .data$inoculum != "" ~ str_replace(.data$inoculum, "^.*_", ""),
        TRUE ~ recode(
          .data$other_frompreviouscolumn,
          "Site1" = "Site 1",
          "Site2" = "Site 2",
          "Site3" = "Site 3",
          "Site4" = "Site 4",
          "Site5" = "Site 5",
          "Site6" = "Site 6",
          .default = NA_character_
        )),
      plot_nmds_color = recode(plot_nmds_color,
                               "W1" = "Site 1",
                               "W2" = "Site 2",
                               "W3" = "Site 3",
                               "W4" = "Site 4",
                               "W5" = "Site 5",
                               "W6" = "Site 6",
                               .default = .data$plot_nmds_color
      ),
      fire_freq = if_else(
        .data$other_frompreviouscolumn %in% paste0("Site", 1:6),
        case_when(
          .data$other_frompreviouscolumn %in% c("Site1", "Site2") ~ 0,
          .data$other_frompreviouscolumn %in% c("Site3", "Site4") ~ 1,
          .data$other_frompreviouscolumn %in% c("Site5", "Site6") ~ 3
        ),
        as.numeric(.data$fire_freq)
      )
    )
}

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

# fungi -------------------------------------------------------------------

fungal_traits <- assign_funguild_to_phyloseq(fung)

final_ecm_ps <- fungal_traits[[1]]

## EcM taxa table ----------------------------------------------------------
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange", "#DA5724", "#508578", "#CD9BCD",
  "#8a592f", "#673770", "#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861",
  
  # Added 5 new distinct colors
  "#1FA187",  # teal
  "#F6E8C3",  # light sand
  "#B8E186",  # soft green (not overlapping your darker greens)
  "#E66101",  # warm orange-brown, distinct from pure orange
  "#4D4D4D"   # neutral charcoal for grounding
)
full_fung <- inoc.full

fung_traits <- assign_funguild_to_phyloseq(full_fung)

ecm_ps_all <- fung_traits[[1]]


#get only EcM
mutualist_guilds <- 
  grep("Ectomycorrhizal",(ecm_ps_all@tax_table[,8]),value = TRUE) %>% 
  unique()

#ignore snowbrush
ecm_ps <-subset_samples(ecm_ps_all, is.na(species) | species == "GrandFir" ) 


ecm_ps <- ecm_ps%>%subset_taxa(guild_fg %in% mutualist_guilds)

ecm_taxa <- as.data.frame(tax_table(ecm_ps)) %>% unique()

ecm_taxa$Genus %>% unique()

ps_colors <- ecm_ps %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.001) %>%                        
  mutate(Family = ifelse(Abundance < 0.01, "Other", Genus))

## more color stuff ---------------------------------------------------


# global genus levels (all figs should use the same order and colors for consistency)
genus_levels <- ps_colors %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Genus)

pal16 <- unname(phylum_colors)

genus_colors <- setNames(
  rep(pal16, length.out = length(genus_levels)),
  genus_levels
)


missing <- setdiff(genus_levels, names(genus_colors))
if (length(missing) > 0) {
  genus_colors[missing] <- "grey90"  # fallback palette
}

## 7 plots -----------------------------------------------------------------
mutualist_guilds <- 
  grep("Ectomycorrhizal",(final_ecm_ps@tax_table[,8]),value = TRUE) %>% 
  unique()

final_ecm_ps <- final_ecm_ps%>%subset_taxa(guild_fg %in% mutualist_guilds)

s_data <- data.frame(sample_data(final_ecm_ps))

tot_reads <- sample_sums(final_ecm_ps)

s_data <- s_data %>% 
  mutate(inoculum_site = replace_na(inoculum_site, "control"),
         total_reads = as.factor(tot_reads[rownames(.)]))
  
sample_data(final_ecm_ps) <- sample_data(s_data)

final_ecm_ps <- final_ecm_ps %>% 
  subset_samples(species == "GrandFir")

a <- data.frame(sample_data(final_ecm_ps))


# plot for all fung -------------------------------------------------------
psa <- fung %>% subset_samples(species == "GrandFir")


s_data <- data.frame(sample_data(psa))

s_data <- s_data %>% 
  mutate(total_reads = as.factor(tot_reads[rownames(.)]))

sample_data(psa) <- sample_data(s_data)

###Remove after test
#final_ecm_ps <- psa
###

ps_genus_20 <- final_ecm_ps %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.02, "Other", Genus))

# make sure Genus has fixed global levels (important for consistent legend/colors)
ps_genus_20 <- ps_genus_20 %>%
  filter(!is.na(Genus), Genus != "") %>%
  mutate(Genus = factor(Genus, levels = genus_levels))

  
#ordered names for ggplot by site 
sites <- sort(unique(ps_genus_20$inoculum_site))

# combine 7 plots
plots_by_site <- setNames(lapply(sites, function(s) {
  df_s <- ps_genus_20 %>% filter(inoculum_site == s)
  
  ggplot(df_s, aes(x = total_reads, y = Abundance, fill = Genus)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = genus_colors, drop = FALSE, na.translate = FALSE) +
    labs(
      x = "Sample",
      y = "Proportion of Genera",
      fill = "Genus",
      title = s
    ) +
    theme_few() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 6),
      legend.key.size = unit(0.5, "lines"),
      legend.spacing = unit(0.5, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1)
)
}), sites)


# --- 1) Reference plot for legend only
p_ref <- ggplot(ps_genus_20, aes(x = sample_name, y = Abundance, fill = Genus)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = genus_colors, drop = FALSE, na.translate = FALSE) +
  labs(fill = "Genus") +
  theme_few() +
  theme(legend.position = "right")

# extract legend from reference
leg <- cowplot::get_legend(p_ref)
legend_panel <- wrap_elements(full = leg)



## fungal inoc data input -------------------------------------------------------------

# Pull sample data
df_samples <- sample_data(inoc.full) %>%
  data.frame()

# Clean 
df_samples <- clean_func(df_samples)

# Keep only inoculum
df_samples <- df_samples %>%
  dplyr::filter(tolower(trimws(community)) == "inoculum")

# Push back into phyloseq
sample_data(inoc.full) <- sample_data(df_samples)

# Prune the phyloseq object to match
inoc_fung_ps <- prune_samples(rownames(df_samples), inoc.full)

fung_inoc_traits <- assign_funguild_to_phyloseq(inoc_fung_ps)

inoc_ecm_ps <- fung_inoc_traits[[1]]
inoc_ecm_guild_df <- fung_inoc_traits[[2]]


## fungal inoc fig --------------------------------------------------------

#get only EcM
mutualist_guilds_inoc <- 
  grep("Ectomycorrhizal",(inoc_ecm_ps@tax_table[,8]),value = TRUE) %>% 
  unique()

a <- data.frame(sample_data(inoc_ecm_ps), check.names = FALSE)
a <- a[order(a$plot_nmds_color), , drop = FALSE]
a$plot_nmds_color_order <- seq_len(nrow(a))

sample_data(inoc_ecm_ps) <- sample_data(a)



inoc_ecm_ps <- inoc_ecm_ps%>%subset_taxa(guild_fg %in% mutualist_guilds_inoc)

ecm_taxa_inoc <- as.data.frame(tax_table(inoc_ecm_ps)) %>% unique()


ecm_taxa$Genus %>% unique()
inoc_genus_20 <- inoc_ecm_ps %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.03, "Other", Genus))



panel_8 <- ggplot(inoc_genus_20)+
  aes(x = plot_nmds_color, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
  geom_col(position = "fill") +  
  scale_fill_manual(values = genus_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  guides(alpha = "none")+
  labs(
    x = "Inoculum site",
    y = "Proportion of Genera",
    fill = "Genus",
    title = "EcM commuinity from Inoculum"
  )+
  theme_few()
panel_8


ggplot(inoc_genus_20)+
  aes(x = fire_freq, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
  #facet_wrap(~plot_nmds_color)+
  geom_col(position = "fill") +  
  scale_fill_manual(values = genus_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  guides(alpha = "none")+
  labs(
    x = "",
    y = "Proportion of Genera",
    fill = "Genus",
    title = "EcM general in final community"
  )+
  theme_few()


## combine for final fig ---------------------------------------------------

plots_by_site[[8]] <- panel_8

# --- 2) Turn off legends in the 7 site plots ---
plots_noleg <- lapply(plots_by_site, \(p) p + theme(legend.position = "none"))

# --- 3) 3x3 grid (7 + 2 blanks) + legend once ---
grid_3x3 <- wrap_plots(
  c(plots_noleg, list(plot_spacer(), plot_spacer())),
  ncol = 3
)

final_plot <- grid_3x3 | legend_panel +
  plot_layout(widths = c(1, 0.30))  # adjust legend column width

final_plot

class(plots_by_site)

# 18S ---------------------------------------------------------------------                     site4.inoc.am,site5.inoc.am,site6.inoc.am)

am_full <- readRDS("./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")

#am_ps <- readRDS("./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")

am_full <- am_full %>%subset_taxa(Phylum %in% c("Glomeromycota", "p__Glomeromycota"))
b<- data.frame(tax_table(am_full))


tot_reads <- sample_sums(am_full)

# add to sample_data
sam <- data.frame(sample_data(am_full)) %>%
  mutate(total_reads = as.factor(tot_reads[rownames(.)]))

sample_data(am_full) <- sample_data(sam)

a1 <- data.frame(sample_data(am_full))


#am_full <-subset_samples(am_full, is.na(species) | species == "Snowbrush" ) 


am_taxa <- as.data.frame(tax_table(am_full)) %>% unique()
am_taxa$Genus %>% unique()

am_subset <- am_full %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.001) %>%                        
  mutate(Family = ifelse(Abundance < 0.01, "Other", Genus))

## 18S colors ---------------------------------------------------

# global genus levels (all figs should use the same order and colors for consistency)
am_genus_levels <- am_subset %>%
  group_by(Genus) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Genus)

am_genus_colors <- setNames(
  rep(pal16, length.out = length(am_genus_levels)),
  am_genus_levels
)


missing <- setdiff(am_genus_levels, names(am_genus_colors))
if (length(missing) > 0) {
  am_genus_colors[missing] <- "grey90"  # fallback palette
}

## 7 plots -----------------------------------------------------------------
am_sam_data <- data.frame(sample_data(am_full))

am_sam_data <- am_sam_data %>% 
  mutate(inoculum_site = replace_na(inoculum_site, "control"))

sample_data(am_full) <- sample_data(am_sam_data)



am_final_subset <- am_full %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.02, "Other", Genus))

# make sure Genus has fixed global levels (important for consistent legend/colors)
am_final_subset <- am_final_subset %>%
  filter(!is.na(Genus), Genus != "") %>%
  mutate(Genus = factor(Genus, levels = am_genus_levels))


#ordered names for ggplot by site 
sites <- sort(unique(am_final_subset$inoculum_site))

# combine 7 plots
plots_by_site <- setNames(lapply(sites, function(s) {
  df_s <- am_final_subset %>% filter(inoculum_site == s)
  
  ggplot(df_s, aes(x = total_reads, y = Abundance, fill = Genus)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = am_genus_colors, drop = FALSE, na.translate = FALSE) +
    labs(
      x = "# reads",
      y = "Proportion of Genera",
      fill = "Genus",
      title = s
    ) +
    theme_few() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 6),
      legend.key.size = unit(0.5, "lines"),
      legend.spacing = unit(0.5, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1))
}), sites)
plots_by_site

# --- 1) Reference plot for legend only
p_ref <- ggplot(am_final_subset, aes(x = sample_name, y = Abundance, fill = Genus)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = am_genus_colors, drop = FALSE, na.translate = FALSE) +
  labs(fill = "Genus") +
  theme_few() +
  theme(legend.position = "right")

# extract legend from reference
leg <- cowplot::get_legend(p_ref)
legend_panel <- wrap_elements(full = leg)



# AM  inoc data input -------------------------------------------------------------
site1.inoc.am <- readRDS("Output/phyloseq_objects/18S_site1.inoc.full.RDS")
site2.inoc.am <- readRDS("Output/phyloseq_objects/18S_site2.inoc.full.RDS")
site3.inoc.am <- readRDS("Output/phyloseq_objects/18S_site3.inoc.full.RDS")
site4.inoc.am <- readRDS("Output/phyloseq_objects/18S_site4.inoc.full.RDS")
site5.inoc.am <- readRDS("Output/phyloseq_objects/18S_site5.inoc.full.RDS")
site6.inoc.am<- readRDS("Output/phyloseq_objects/18S_site6.inoc.full.RDS")

am_inoc <- merge_phyloseq(site1.inoc.am,site2.inoc.am,site3.inoc.am,
                           site4.inoc.am,site5.inoc.am,site6.inoc.am)

am_inoc <- am_inoc %>%subset_taxa(Phylum %in% c("Glomeromycota", "p__Glomeromycota")) %>% 
  subset_samples(community == "inoculum")

tot_reads <- sample_sums(am_inoc)

# add to sample_data
sam <- data.frame(sample_data(am_inoc)) %>%
  mutate(total_reads = as.factor(tot_reads[rownames(.)]))

sample_data(am_inoc) <- sample_data(sam)

a1 <- data.frame(sample_data(am_inoc))
# #fungal inoc fig --------------------------------------------------------


inoc_df <- am_inoc %>%
  tax_glom(taxrank = "Genus") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.03, "Other", Genus))



panel_8 <- ggplot(inoc_df)+
  aes(x = other_frompreviouscolumn, y = Abundance, fill = fct_rev(fct_infreq(Genus))) +
  geom_col(position = "fill") +  
  scale_fill_manual(values = am_genus_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  guides(alpha = "none")+
  labs(
    x = "Inoculum site",
    y = "Proportion of Genera",
    fill = "Genus",
    title = "Commuinity from Inoculum"
  )+
  theme_few()
panel_8



# combine for final fig ---------------------------------------------------

plots_by_site[[8]] <- panel_8

# --- 2) Turn off legends in the 7 site plots ---
plots_noleg <- lapply(plots_by_site, \(p) p + theme(legend.position = "none"))

# --- 3) 3x3 grid (7 + 2 blanks) + legend once ---
grid_3x3 <- wrap_plots(
  c(plots_noleg, list(plot_spacer(), plot_spacer())),
  ncol = 3
)

final_plot <- grid_3x3 | legend_panel +
  plot_layout(widths = c(1, 0.30))  # adjust legend column width

final_plot+
  labs(
    title = "AM community composition by inoculum site"
  )

class(plots_by_site)



# useless ITS AM analysis below, keep for now I guess  ---------------------------------------------------------
ITS_am <- fung %>% subset_samples(species == "Snowbrush")

just_am <- ITS_am %>%
  subset_taxa(Phylum %in% c("Glomeromycota", "p__Glomeromycota"))

a <- data.frame(tax_table(just_am))


am_colors <- just_am %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.001) %>%                        
  mutate(Family = ifelse(Abundance < 0.01, "Other", Family))

## am colors ---------------------------------------------------


# global genus levels (all figs should use the same order and colors for consistency)
am_fam_levels <- am_colors %>%
  group_by(Family) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Family)

am_colors <- setNames(
  rep(pal16, length.out = length(am_fam_levels)),
  am_fam_levels
)

am_toplot <- ITS_am %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x / sum(x)}) %>% 
  psmelt() %>%                                         
  filter(Abundance > 0.01) %>%                        
  mutate(Family = ifelse(Abundance < 0.02, "Other", Family))


ggplot(am_toplot, aes(x = sample_name, y = Abundance, fill = Family))+
  geom_col(position = "fill") +  
  facet_wrap(~inoculum_site)+
  scale_fill_manual(values = am_colors) +
  theme(legend.position = "right",                 
        legend.title = element_text(size = 8),      
        legend.text = element_text(size = 6),      
        legend.key.size = unit(0.5, "lines"),       
        legend.spacing = unit(0.5, "lines")) +  
  guides(alpha = "none")+
  theme_few()
  
  
  
  
  
  
  
