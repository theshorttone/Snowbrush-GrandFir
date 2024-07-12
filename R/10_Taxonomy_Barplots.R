# Setup ####
library(tidyverse)
library(phyloseq)
library(zahntools)

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
  merge_samples('site',fun = 'sum') %>% 
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
  ggplot(aes(x=Sample,y=Abundance,fill=Family)) +
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
