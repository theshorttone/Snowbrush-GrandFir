# Explore ITS2 data

# this had pretty weak 5' ends and had to be cut off at 150nt

# -----------------------------------------------------------------------------#
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
# -----------------------------------------------------------------------------#

# SETUP ####

# packages 
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(corncob); packageVersion("corncob")
library(patchwork)
library(igraph)
# Options
options(scipen=999)

# Data
ps <- readRDS("./Output/ITS_clean_phyloseq_object.RDS") %>% 
  subset_samples(species == "GrandFir") 

# ALPHA-DIV ESTIMATES ####
alpha <- estimate_richness(ps) %>% 
  select(Observed,Shannon, Simpson)

variables <- names(sample_data(ps))

plot_richness(ps,
              measures = c("Shannon"),
              color="inoculum_site",sortby = "Shannon") +
  theme(legend.title = element_text(hjust=.5),
        legend.position = "bottom") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
ggsave("./Output/ITS_richness_by_inoculum.png",height = 4,width = 4)

# add alpha diversity estimates to metadata for modeling
alpha_df <- microbiome::meta(ps) %>%
  bind_cols(alpha)

# plot effect of alpha div on leaf number
alpha_df %>% 
  ggplot(aes(x=Observed,y=as.numeric(leaf_number))) +
  geom_point()

# custom plot
p <- 
  alpha_df %>% 
  pivot_longer(c(Observed,Shannon,Simpson),names_to = "Measure") %>% 
  ggplot(aes(x=inoculum_burn_freq,y=value)) +
  geom_boxplot() +
  facet_wrap(~Measure,scales = 'free') +
  theme_minimal()
p
saveRDS(p,"./Output/figs/16S_burn_frequency_and_alpha-div.RDS")
ggsave("./Output/ITS_burn_frequency_and_alpha-div.png",height = 4,width = 6,dpi=300)

p <- 
  alpha_df %>% 
  pivot_longer(c(Observed,Shannon,Simpson),names_to = "Measure") %>% 
  ggplot(aes(x=inoculum_site,y=value)) +
  geom_boxplot() +
  facet_wrap(~Measure,scales = 'free') +
  theme_minimal()
p
saveRDS(p,"./Output/figs/16S_inoc_source_and_alpha-div.RDS")
ggsave("./Output/ITS_inoc_source_and_alpha-div.png",height = 4,width = 6,dpi=300)

# regression on alpha diversity with burn freq predictor
tuk <- alpha_df %>% 
  aov(data=.,formula= Simpson ~ inoculum_site) %>% 
  TukeyHSD()
tuk$inoculum_site %>%
  as.data.frame() %>% 
  dplyr::filter(`p adj` < 0.05) %>% 
  ggplot(aes(x=row.names(.),y=diff)) +
  geom_hline(yintercept = 0, linetype=2, color='red') +
  geom_point() +
  geom_errorbar(aes(ymax=upr,ymin=lwr)) +
  labs(x="Inoculum Pair Comparison",y="Alpha Diversity Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(face='bold',size=12))
ggsave("Output/figs/ITS_Simpson_Div_Difference_Between_Inoc_Source.png",height = 4,width = 4)

# Beta-diversity ####
ord <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS")

(
  ord1 <- 
    plot_ordination(ps,ord,color="inoculum_site") +
    stat_ellipse() +
    theme_minimal() +
    labs(color="Inoculum source")
)
ggsave("./Output/ITS_Ordination_by_inoculum_source.png",height = 4,width = 4)
saveRDS(ord1,"./Output/ITS_Ordination_by_inoculum_source.RDS")

(
  ord2 <- 
    plot_ordination(ps,ord,color="inoculum_burn_freq") +
    stat_ellipse() +
    theme_minimal() +
    labs(color="Inoculum burn\nfrequency")
)
ggsave("./Output/ITS_Ordination_by_inoculum_burnfreq.png",height = 4,width = 4)
saveRDS(ord2,"./Output/ITS_Ordination_by_inoculum_burnfreq.RDS")

# Overlay environmental soil variables with ordinations ####

# get dataframe of just soil variables
soilvars <- microbiome::meta(ps) %>% 
  dplyr::select(starts_with("mean_"))

# fit environmental vars to NMDS
env <- envfit(ord,soilvars,na.rm=TRUE)
plot(env)
plot(ord)
plot(env)

# extract info from ordination and envfit
ord_scores <- scores(ord)$sites %>% as.data.frame
ord_df <- 
  ord_scores %>% 
  bind_cols(soilvars)
ord_df$inoculum_burn_freq <- ps@sam_data$inoculum_burn_freq

env_scores <- 
  as.data.frame(scores(env,"vectors")) * ordiArrowMul((env))

envfit_plot_ITS <- 
  ggplot(data = ord_df, aes(x=NMDS1,y=NMDS2,color=inoculum_burn_freq)) +
  geom_point(alpha=.7) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = env_scores, linewidth =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = env_scores, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(env_scores), colour = "navy", size=2) +
  coord_cartesian(c(-2.5,3.5)) +
  theme_minimal()
saveRDS(envfit_plot_ITS,"./Output/figs/envfit_plot_ITS.RDS")


# covariation between inoculum site and burn freq?
table(alpha_df$inoculum_burn_freq, alpha_df$inoculum_site) %>% 
  saveRDS("./Output/ITS_inoculum_nestedness.RDS")

# definitely yes!

# Check with permanova
# relabund table
ra_table <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% 
  as.data.frame()
adonis2(data=alpha_df,
        formula = ra_table ~ ps@sam_data$inoculum_burn_freq * ps@sam_data$inoculum_site) %>% 
  broom::tidy() %>% 
  saveRDS("./Output/ITS_PermANOVA_inoculum_source.RDS")
# inoculum site explains about 27% of community variation!

# These are all the same for bacteria and fungi ....
# removing burn freq as a variable...
# inoc_leaf_mod
# glm(data = alpha_df,
#     formula = leaf_number ~ inoculum_site) %>% summary 
#   saveRDS("./Output/ITS_inoc_leafnumber_mod.RDS")
# 
# # inoc_rootdm_mod <- 
# glm(data = alpha_df,
#     formula = final_root_dm ~ inoculum_site) %>% 
#   saveRDS("./Output/ITS_inoc_rootdm_mod.RDS")
# 
# glm(data = alpha_df,
#     formula = shoot_dm ~ inoculum_site) %>% 
#   saveRDS("./Output/ITS_inoc_shootdm_mod.RDS")
# 
# glm(data = alpha_df,
#     formula = height ~ inoculum_site) %>% 
#   saveRDS("./Output/ITS_inoc_height_mod.RDS")

saveRDS(alpha_df,"./Output/ITS_alpha_df.RDS")

# direct effect of community on leaf number? (backwards)
adonis2(data=alpha_df,
        formula = ra_table ~ ps@sam_data$leaf_number, strata = ps@sam_data$block) %>% 
  broom::tidy() %>% 
  saveRDS("./Output/ITS_PermANOVA_leaf_number.RDS")
# community explains about 3% of leaf number vairance within blocks

# make taxa names...taxonomy
taxa_print <- corncob::otu_to_taxonomy(taxa_names(ps),ps)
colnames(ra_table) <- janitor::make_clean_names(taxa_print)
# subset relabund table to present taxa
ra_table <- ra_table[,colSums(ra_table)>0]
df <- ra_table
df$leaf_number <- ps@sam_data$leaf_number

colnames(df)

rf_mod <- ranger::ranger(formula = leaf_number ~ .,
                         data=df,importance = 'permutation',num.trees = 999)

summary(rf_mod)
its_vip_plot <- vip::vip(rf_mod) + 
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(title="Top ten most important taxa detected by Random-Forest model",
       subtitle = "Explains 7.6% of leaf number variation.",
       caption = "k_fungi most likely refers to an undescribed species")
ggsave("./Output/ITS_RF_mod_VIP_leaf_count.png",width = 18,height = 5)
saveRDS(its_vip_plot,"./Output/ITS_RF_mod_VIP_leaf_count.RDS")

# do same as above, but with taxa grouped at genus/species level
ra_table2 <- ps %>% 
  tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% 
  as.data.frame()

taxa_print <- corncob::otu_to_taxonomy(taxa_names(ps %>% tax_glom("Genus")),ps %>% tax_glom("Genus"))
colnames(ra_table2) <- janitor::make_clean_names(taxa_print)


# subset relabund table to present taxa
ra_table2 <- ra_table2[,colSums(ra_table2)>0]
df <- ra_table2
df$leaf_number <- ps@sam_data$leaf_number

colnames(df)

rf_mod <- ranger::ranger(formula = leaf_number ~ .,
                         data=df,importance = 'permutation',num.trees = 999)

(rf_mod)
its_vip_plot2 <- 
vip::vip(rf_mod) + 
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(title="Top ten most important genera detected by Random-Forest model",
       subtitle = "Explains 18.7% of leaf number variation.")
ggsave("./Output/ITS_RF_mod_VIP_genus-leaf_count.png",width = 18,height = 5)
saveRDS(its_vip_plot2,"./Output/ITS_RF_mod_VIP_genus-leaf_count.RDS")

# what about the good-looking inoculum #4?

# build new rf table with logical outcome for site 4
df <- ra_table
df$site4 <- ps@sam_data$inoculum_site %in% c("4")

# run rf model
rf_mod <- ranger::ranger(formula = site4 ~ .,
                         data=df,importance = 'permutation',num.trees = 999,classification = TRUE)

(rf_mod)
vip::vip(rf_mod) + 
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(title="Top ten most important genera detected by Random-Forest model")
ggsave("./Output/ITS_RF_mod_VIP_site4.png",width = 18,height = 5)


# Heatmap ####
# full relabund transformation
ps_ra_genus <- 
  ps %>% 
  tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)})

(
  p3 <- 
    ps_ra_genus %>% 
    subset_taxa(taxa_sums(ps_ra_genus) > 0.1) %>% 
    plot_heatmap(sample.order = "inoculum_site",sample.label = "inoculum_site") +
    scale_fill_viridis_c(option = "magma") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),axis.text.x = element_blank(),
          strip.text = element_text(face='bold',size=14),
          axis.title.x = element_blank()) +
    labs(fill="Relative\nabundance",y="Most abundant ASVs",title="Comparison of most abundant bacterial ASVs\nin samples from each inoculum source") +
    facet_wrap(~inoculum_site,scales = 'free_x')
  
)
saveRDS(p3,"./Output/ITS_heatmap_mostabund.RDS")


# Corncob analysis ####

# merge taxa at genus level
ps_genus <- ps %>% 
  tax_glom(taxrank = "Genus")

# Clean up ASV names to show taxonomy
ASV_names <- otu_table(ps) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Phylum","Class","Order","Family","Genus"))

genus_names <- otu_table(ps_genus) %>% colnames()
genus_taxa <- otu_to_taxonomy(genus_names,ps_genus,level = c("Phylum","Class","Order","Family","Genus"))

# CORNCOB DIFFABUND ####
ps_genus@sam_data$inoculum_site %>% unique

# set levels of inoculum site so "sterile" is intercept
ps_genus@sam_data$inoculum_site <- 
  ps_genus@sam_data$inoculum_site %>% factor(levels = c("Sterile","1","2","3","4","5","6"))

# use raw count data for corncob
da_analysis_inocsite <- differentialTest(formula = ~ inoculum_site, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_genus,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
p4 <- 
  plot(da_analysis_inocsite) +
  theme(legend.position = 'none') 
saveRDS(p4,"./Output/ITS_diffabund_inocsite.RDS")

# quick look at Rhizobium and Phenylobacterium and Neorhizobium
# 
# rhizobium_asv <- 
#   ps_genus@tax_table[which(ps_genus@tax_table[,6] == "Rhizobium"),] %>% row.names()
# phenylobacterium_asv <- 
#   ps_genus@tax_table[which(ps_genus@tax_table[,6] == "Phenylobacterium"),] %>% row.names()
# neorhizobium_asv <- 
#   ps_genus@tax_table[which(ps_genus@tax_table[,6] == "Neorhizobium"),] %>% row.names()
# 
sig_taxa <- da_analysis_inocsite$significant_taxa %>% otu_to_taxonomy(data=ps_genus)
bbdml_obj <- multi_bbdml(da_analysis_inocsite,
                         ps_object = ps_genus,
                         mu_predictor = "inoculum_site",
                         phi_predictor = "inoculum_site",
                         taxlevels = 4:6)

# names(bbdml_obj)

# plot 2 significant taxa, for now (rhizobium and phenylobacterium)
plot_multi_bbdml(bbdml_obj,color = "inoculum_site",pointsize = 3)  
library(patchwork)
(bbdml_plot_1  + theme(axis.title.y = element_blank(),legend.position = 'none')) / 
  (bbdml_plot_2  + theme(axis.title.y = element_blank(),legend.position = 'none')) / 
  (bbdml_plot_3  + theme(axis.title.y = element_blank(),legend.position = 'none')) / 
  (bbdml_plot_4  + theme(axis.title.y = element_blank(),legend.position = 'right')) /
  (bbdml_plot_5  + theme(axis.title.y = element_blank(),legend.position = 'none')) / 
  (bbdml_plot_6  + theme(axis.title.y = element_blank(),legend.position = 'none')) / 
  (bbdml_plot_7  + theme(axis.title.y = element_blank(),legend.position = 'none')) / 
  (bbdml_plot_8  + theme(axis.title.y = element_blank(),legend.position = 'none'))
ggsave("./Output/ITS_bbdml_plots.png",height = 16,width = 8)


# NETWORK ANALYSIS ####
set.seed(666)
fung_network <- ps %>% 
  make_network(keep.isolates = TRUE)

# quick plot
plot_net(ps,distance = 'jaccard',color = 'inoculum_site',rescale = TRUE,maxdist = .7,point_label = "inoculum_burn_freq")

fung_network

igraph::alpha_centrality()
igraph::assortativity_degree(fung_network,directed = FALSE)
deg <- igraph::degree(fung_network)
tmax <- igraph::centr_degree_tmax((fung_network),loops = FALSE)
igraph::centralize(deg, tmax)
igraph::cliques(fung_network)
igraph::assortativity(fung_network,
              types1 = ps@sam_data$inoculum_site %>% factor %>% as.numeric)
igraph::largest_cliques(fung_network)
