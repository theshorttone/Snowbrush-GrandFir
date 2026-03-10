# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(janitor); packageVersion("janitor")
library(microbiome); packageVersion("microbiome")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(ecodist); packageVersion("ecodist")
library(corncob); packageVersion("corncob")
library(patchwork); packageVersion("patchwork")
library(zahntools); packageVersion('zahntools') 

# seed
set.seed(666)

## Load Data ####

ps_bact <- readRDS("Output/phyloseq_objects/16S_inoculum_samples_clean_phyloseq_object.RDS")
ps_fung <- readRDS("Output/phyloseq_objects/ITS_inoculum_samples_clean_phyloseq_object.RDS")
ps_am <- readRDS("Output/phyloseq_objects/18S_inoculum_samples_clean_phyloseq_object.RDS")

library(phyloseq)
library(ggplot2)

### 16S (Bacteria)
ord_bact <- ps_bact %>% 
  transform_sample_counts(function(x) x/sum(x)) %>% 
  ordinate(method = "NMDS",distance = "bray")

bact_plot_inoc <- plot_ordination(ps_bact, ord_bact, color = "other_frompreviouscolumn") +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, lwd = 1) +
  theme_classic(base_size = 14)

bact_plot_inoc
saveRDS(bact_plot_inoc, "./Output/figs/16S_inoc_nmds_1.RDS")


### ITS (Fungi)
ord_fung <- ps_fung %>% 
  transform_sample_counts(function(x) x/sum(x)) %>% 
  ordinate(method = "NMDS",distance = "bray")

fung_plot_inoc <- plot_ordination(ps_fung, ord_fung, color = "other_frompreviouscolumn") +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, lwd = 1) +
  theme_classic(base_size = 14)

fung_plot_inoc
saveRDS(fung_plot_inoc, "./Output/figs/ITS_nmds_all.RDS")

ord_am <- 
  ps_am %>% 
  ordinate(method = "NMDS",distance = "bray")
am_plot_inoc <- plot_ordination(ps_am, ord_am, color = "other_frompreviouscolumn") +
  geom_point(size = 3, alpha = 0.8) +        # bigger dots
  stat_ellipse(level = 0.95, lwd = 1) +      # elliptical group circles
  theme_classic(base_size = 14)
am_plot_inoc; saveRDS(bact_plot_inoc,"./Output/figs/18S_inoc_nmds_all.RDS")
  
  
