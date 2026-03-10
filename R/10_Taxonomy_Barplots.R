# Setup ####
library(tidyverse)
library(phyloseq)
library(zahntools)
library(vegan); packageVersion("vegan")
library(ggnewscale)


fungal_traits <- readRDS("./Output/ITS_Fungal_Guild_Assignments.RDS")

bact_traits <- readRDS("./Output/16S_Bacterial_Trait_Database.RDS")

# 6 site colors (your extremes)
site_colors <- c(
  "Site1" = "#4cbfe6", "Site2" = "#2443f0",
  "Site3" = "#f0c424", "Site4" = "#db7e04",
  "Site5" = "#cc6866", "Site6" = "#9e0402"
)

# 3 mid-blend colors (averages) for ellipses
burn_colors <- c(
  "unburned" = "#3881F3",
  "1 burn"   = "#E59D14",
  "3 burn"   = "#A73634"
)
# functions ####
'%ni' <- Negate('%in%')

ra <- function(x){x/sum(x)}

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




# Prep for barplots ####
melt_16S <- inoc_16S %>% 
  tax_glom("Phylum") %>% 
  merge_samples('site', fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  psmelt()

phyorder_16S <- 
melt_16S %>% 
  group_by(Phylum) %>% 
  summarize(mean_ra = mean(Abundance)) %>% 
  arrange(desc(mean_ra)) %>% 
  pluck("Phylum")

melt_16S <- 
  melt_16S %>% 
  mutate(Phylum = factor(Phylum,levels=phyorder_16S))

melt_ITS <- inoc_ITS %>% 
  tax_glom("Phylum") %>% 
  merge_samples('site',fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  psmelt()

melt_ITS <- 
melt_ITS %>% 
  mutate(Phylum = Phylum %>% str_remove("p__"))
phyorder_ITS <- 
  melt_ITS %>% 
  group_by(Phylum) %>% 
  summarize(mean_ra = mean(Abundance)) %>% 
  arrange(desc(mean_ra)) %>% 
  pluck("Phylum")
melt_ITS <- 
  melt_ITS %>% 
  mutate(Phylum = factor(Phylum,levels=phyorder_ITS))

melt_18S <- inoc_18S %>% 
  subset_taxa(Phylum == "p__Glomeromycota") %>% 
  tax_glom("Family") %>% 
  merge_samples('site',fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% 
  mutate(Family = Family %>% str_remove("f__"))

phyorder_18S <- 
  melt_18S %>% 
  group_by(Family) %>% 
  summarize(mean_ra = mean(Abundance)) %>% 
  arrange(desc(mean_ra)) %>% 
  pluck("Family")
melt_18S <- 
  melt_18S %>% 
  mutate(Family = factor(Family,levels=phyorder_18S))


melt_am <- ps_am %>% 
  transform_sample_counts(ra) %>%     # convert to relative abundance
  psmelt()                      # long tidy format

# barplot of each individual sample
ggplot(melt_am, aes(x = sample_name, y = Abundance, fill = Genus)) +
  geom_col() +
  facet_wrap(~ sample_name, scales = "free_x") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Genus")




# Zac prep for Barplots ------------------------------------------------------------
melt_16S_fam <- inoc_16S %>% 
  tax_glom("Family") %>% 
  merge_samples('site', fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  psmelt()

keep_fams_16S <- melt_16S_fam %>% 
  group_by(Family) %>% 
  summarize(mean_ra = mean(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  filter(mean_ra > 0.03) %>% 
  arrange(desc(mean_ra)) %>% 
  pull(Family)

plot_df_16S <- melt_16S_fam %>% 
  filter(Family %in% keep_fams_16S) %>% 
  mutate(
    Family = factor(Family, levels = keep_fams_16S),
    site   = factor(site, levels = unique(site))
  )

melt_ITS <- inoc_ITS %>% 
  tax_glom("Family") %>% 
  merge_samples('site', fun = 'sum') %>% 
  transform_sample_counts(ra) %>% 
  psmelt()

keep_fams_ITS <- melt_ITS %>% 
  group_by(Family) %>% 
  summarize(mean_ra = mean(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  filter(mean_ra > 0.03) %>% 
  arrange(desc(mean_ra)) %>% 
  pull(Family)
       
plot_df_ITS <- melt_ITS %>% 
  filter(Family %in% keep_fams_ITS) %>% 
  mutate(
    Family = factor(Family, levels = keep_fams_ITS),
    site   = factor(site, levels = unique(site))
  )



# family_bargraphs --------------------------------------------------------

#16S
plot_df_16S %>% 
  ggplot(aes(x=Sample,y=Abundance,fill = Family)) +
  geom_col() +
  theme_bw() +
  labs(x="Inoculum",y="Relative abundance") +
  scale_fill_viridis_d(end=.9,option = 'turbo') +
  theme(legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face='bold',size=12))
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_family_16S.png",dpi=300,height = 6,width = 8)
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_family_16S.tiff",dpi=600,height = 6,width = 8)

#ITS
plot_df_ITS %>% 
  ggplot(aes(x=Sample,y=Abundance,fill = Family)) +
  geom_col() +
  theme_bw() +
  labs(x="Inoculum",y="Relative abundance") +
  scale_fill_viridis_d(end=.9,option = 'turbo') +
  theme(legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face='bold',size=12))
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_family_ITS.png",dpi=300,height = 6,width = 8)
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_family_ITS.tiff",dpi=600,height = 6,width = 8)



# Build barplots ####
melt_16S %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Phylum)) +
  geom_col() +
  theme_bw() +
  labs(x="Inoculum",y="Relative abundance") +
  scale_fill_viridis_d(end=.9,option = 'turbo') +
  theme(legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face='bold',size=12))
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_16S.png",dpi=300,height = 6,width = 8)
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_16S.tiff",dpi=600,height = 6,width = 8)


melt_ITS %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Phylum)) +
  geom_col() +
  theme_bw() +
  labs(x="Inoculum",y="Relative abundance") +
  scale_fill_viridis_d(end=.9,option = 'turbo') +
  theme(legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face='bold',size=12))
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_ITS.png",dpi=300,height = 6,width = 8)
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_ITS.tiff",dpi=600,height = 6,width = 8)


melt_18S %>% 
  ggplot(aes(x=,y=Abundance,fill=Family)) +
  geom_col() +
  theme_bw() +
  labs(x="Inoculum",y="Relative abundance") +
  scale_fill_viridis_d(end=.9,option = 'turbo') +
  theme(legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text.x = element_text(face='bold',size=12,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face='bold',size=12))
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_18S.png",dpi=500,height = 6,width = 8)
ggsave("./Output/figs/manuscript_versions/Inoculum_Barplot_18S.tiff",dpi=600,height = 6,width = 8)

# Quick diversity stats ####

data.frame(bact_shannon = inoc_16S %>% 
             estimate_richness() %>% 
             pluck('Shannon'),
           fung_shannon = inoc_ITS %>% 
             estimate_richness() %>% 
             pluck('Shannon')
           ) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(x=name,y=value)) +
  geom_boxplot()

bact_shannon = inoc_16S %>% 
             estimate_richness() %>% 
             pluck('Shannon')
fung_shannon = inoc_ITS %>% 
             estimate_richness() %>% 
             pluck('Shannon')

t.test(bact_shannon,fung_shannon)

# Hub taxa ####
hubs <- read_csv("./Output/Full_HubTaxa_List.csv") %>% 
  dplyr::filter(Amplicon != "AMF")
fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")
bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS")

# subset full physeqs to hub taxa only
fung_hubs <- 
fung %>% 
  subset_taxa(taxa_names(fung) %in% hubs$ASV)

bact_hubs <- 
  bact %>% 
  subset_taxa(taxa_names(bact) %in% hubs$ASV)

# reduce those to only samples with those taxa
fung_hubs <- 
  fung_hubs %>% 
  subset_samples(sample_sums(fung_hubs) > 0)
bact_hubs <- 
  bact_hubs %>% 
  subset_samples(sample_sums(bact_hubs) > 0)
taxa_sums(bact_hubs)


plot_bar2(bact_hubs,fill = "Genus")
sample_sums(bact)



# NMDS figs ---------------------------------------------------------------

## Load Data ####

ps_bact <- readRDS("Output/phyloseq_objects/16S_inoculum_samples_clean_phyloseq_object.RDS")
ps_fung <- readRDS("Output/phyloseq_objects/ITS_inoculum_samples_clean_phyloseq_object.RDS")
ps_am <- readRDS("Output/phyloseq_objects/18S_inoculum_samples_clean_phyloseq_object.RDS")



# seed
set.seed(666)

# SETUP ####
plot_inoc_nmds <- function(ps, plot_title = "") {
  
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  
  ord <- ordinate(ps_rel, method = "NMDS", distance = "bray")
  
  
  
  ord_df <- cbind(
    as.data.frame(ord$points),
    as.data.frame(sample_data(ps_rel))
  )
  
  # rename first two columns (the ordination axes)
  names(ord_df)[1:2] <- c("NMDS1", "NMDS2")
  
  
  # Burn grouping
  ord_df <- ord_df %>%
    mutate(burn_severity = case_when(
      other_frompreviouscolumn %in% c("Site1", "Site2") ~ "unburned",
      other_frompreviouscolumn %in% c("Site3", "Site4") ~ "1 burn",
      other_frompreviouscolumn %in% c("Site5", "Site6") ~ "3 burn",
      TRUE ~ NA_character_
    )) %>%
    mutate(burn_severity = factor(burn_severity, 
                                  levels = c("unburned","1 burn","3 burn")))
  
  p <- ggplot(ord_df, aes(NMDS1, NMDS2)) +
    geom_point(aes(color = other_frompreviouscolumn), size = 3, alpha = 0.8) +
    scale_color_manual(values = site_colors, name = "Site") +
    ggnewscale::new_scale_color() +
    stat_ellipse(aes(color = burn_severity), level = 0.95, linewidth = 1) +
    scale_color_manual(values = burn_colors, name = "Burn severity") +
    theme_classic(base_size = 14) +
    labs(title = plot_title, x = "NMDS1", y = "NMDS2")
  
  
  # Add stress (if available)
  if (!is.null(ord$stress)) {
    p <- p + annotate(
      "text", x = Inf, y = -Inf,
      label = paste0("Stress = ", round(ord$stress, 3)),
      hjust = 1.1, vjust = -0.8, size = 4
    )
    
  }
  
  dist_mat <- phyloseq::distance(ps_rel, method = "bray")
  
  adon_burn <- vegan::adonis2(
    dist_mat ~ burn_severity,
    data = ord_df,
    permutations = 999
  )
  
  disp_burn <- vegan::betadisper(dist_mat, ord_df$burn_severity)
  disp_burn_anova <- anova(disp_burn)
  
  # return everything
  list(
    plot = p,
    adonis = adon_burn,
    betadisper = disp_burn,
    betadisper_anova = disp_burn_anova
  )
  
}
# plotting ----------------------------------------------------------------
bact <- plot_inoc_nmds(ps_bact, "Bacteria (16S)"); bact
fung <- plot_inoc_nmds(ps_fung, "Fungal (ITS)"); fung
#am <- plot_inoc_nmds(ps_am, "AM (18S)"); am


ps_am_genus <- tax_glom(ps_am, "Genus", NArm = FALSE)

ord_am_genus <- ps_am_genus %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  ordinate(method = "NMDS", distance = "bray", trymax = 200)

am <- plot_inoc_nmds(ps_am_genus, plot_title = "AM (18S) NMDS (Genus level)"); am







