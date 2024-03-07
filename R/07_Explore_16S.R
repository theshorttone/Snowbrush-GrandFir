# -----------------------------------------------------------------------------#
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
#                     zahntools v 0.1.0 (github: gzahn/zahntools)
#                     compositions v 2.0.6
#                     broom.mixed v 0.2.9.4
#                     broom v 1.0.5
#                     ranger v 0.15.1
#                     vip v 0.4.1
#                     janitor v 2.2.0
# -----------------------------------------------------------------------------#

# SETUP ####

## packages ####
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(lmerTest); packageVersion("lmerTest")
library(corncob); packageVersion("corncob")
library(zahntools); packageVersion("zahntools")
library(compositions); packageVersion("compositions")
library(broom.mixed); packageVersion("broom.mixed")
library(broom); packageVersion("broom")
library(ranger); packageVersion("ranger")
library(vip); packageVersion("vip")
library(janitor); packageVersion("janitor")

## options and functions ####
options(scipen=999)
'%ni%' <- Negate('%in%')
theme_set(theme_minimal())

source("./R/palettes.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]

set.seed(666)

plot_topten_relabund <- function(x){
  dat <- 
    x %>% 
    psmelt() %>% 
    mutate(full_taxa_name= paste(Phylum,Class,Order,Family,Genus)) 
  
  # set factor levels for plotting
  dat$full_taxa_name <- factor(dat$full_taxa_name,
                               levels = dat %>% 
                                 group_by(full_taxa_name) %>% 
                                 summarize(A = sum(Abundance)) %>% 
                                 arrange(A) %>% 
                                 pluck("full_taxa_name")
  )
    
  dat %>% 
    ggplot(aes(x=full_taxa_name,y=Abundance,fill=Sample)) +
    geom_col(position='dodge')+
    coord_flip() +
    labs(x="Taxon",y="Relative abundance") +
    theme(axis.text.y = element_text(face='bold.italic'))
}



multi_bbdml_st <- function(da_analysis, ps_object, mu_predictor, phi_predictor, 
    seed = 123, taxlevels = 1:7,sig_taxa=NA) 
{
  if (!class(da_analysis) == "differentialTest") {
    stop("da_analysis must be of 'differentialTest' class")
  }
  else if (!class(ps_object) == "phyloseq") {
    stop("ps_object must be of 'phyloseq' class")
  }
  else if (!class(mu_predictor) == "character" & length(mu_predictor) == 
           1) {
    stop("mu_predictor must be 'character' vector of length 1")
  }
  else if (!class(phi_predictor) == "character" & length(phi_predictor) == 
           1) {
    stop("phi_predictor must be 'character' vector of length 1")
  }
  else if (!class(seed) %in% c("numeric", "integer") & length(seed) == 
           1) {
    stop("random seed must be 'numeric' or 'integer' vector of length 1")
  }
  else if (!length(taxlevels) > 0 & !class(taxlevels) %in% 
           c("integer", "numeric")) {
    stop("'taxlevels' must be a numeric vector of length > 0 that indicated which taxonomic levels to include as a name for significant sequences. Default is 1:7")
  }
  else {
    
    if(any(is.na(sig_taxa))){
      sig_taxa <- unlist(da_analysis["significant_taxa"])
    } 
    forms <- paste0(sig_taxa, " ~ ", mu_predictor)
    forms <- lapply(forms, as.formula)
    phi.form <- as.formula(paste0(" ~ ", phi_predictor))
    bbdml_list <- list()
    for (i in 1:length(forms)) {
      set.seed(seed)
      bb <- bbdml(formula = forms[[i]], phi.formula = phi.form, 
                  data = ps_object)
      bbdml_list[[i]] <- bb
      tax_name <- paste(tax_table(ps_object)[da_analysis[["significant_taxa"]], 
                                             taxlevels][i], sep = "_")
      taxon_name <- paste(tax_name[1], tax_name[2], tax_name[3], 
                          tax_name[4], tax_name[5], tax_name[6], tax_name[7], 
                          sep = "_")
      taxon_name <- str_remove_all(taxon_name, "_NA")
      names(bbdml_list)[[i]] <- taxon_name
    }
    return(bbdml_list)
  }
}

plot_multi_bbdml <- function (bbdml_list, color = "none", obj_basename = "bbdml_plot_", 
                              pointsize = 1, whichtaxa = 1:length(bbdml_list)) 
{
  pal.discrete <- c("#c1593c", "#688e52", "#643d91", "#894e7d", 
                    "#477887", "#12aa91", "#705f36", "#8997b2", "#c4a113", 
                    "#753c2b", "#3c3e44", "#b3bf2d", "#82b2a4", "#820616", 
                    "#a17fc1", "#262a8e", "#abb5b5", "#000000", "#493829", 
                    "#816C5B", "#A9A18C", "#613318", "#855723", "#B99C6B", 
                    "#8F3B1B", "#D57500", "#DBCA69", "#404F24", "#668D3C", 
                    "#BDD09F", "#4E6172", "#83929F", "#A3ADB8")
  if (!class(color) %in% c("character", "NULL") & length(color) == 
      1) {
    stop("variable to color by must be character vector of length 1")
  }
  else if (!class(bbdml_list[[1]]) == "bbdml") {
    stop("bbdml_list must be a list of bbdml objects; typically the output from multi_bbdml()")
  }
  else if (!color %in% c(as.character(bbdml_list[[1]]$formula), 
                         "none")) {
    stop("varible to color by is not found in bbdml model outputs")
  }
  else if (!class(whichtaxa) %in% c("numeric", "integer") & 
           length(whichtaxa) > 0) {
    stop("'whichtaxa' should be a numeric vector with length > 0, indicating which significant taxa to plot. Default is all of them.")
  }
  else if (color == "none") {
    plot_list <- list()
    for (i in whichtaxa) {
      p <- plot(bbdml_list[[i]], size = pointsize) + ggtitle(names(bbdml_list)[i])
      plot_list[[i]] <- p
    }
  }
  else if (length(unique(color)) > length(pal.discrete)) {
    stop("You probably have too many discrete levels of your grouping variable to color by. \nThis isn't a hard no, but if you want to make it work, you'll have to do it manually.")
  }
  else {
    plot_list <- list()
    for (i in 1:length(bbdml_list)) {
      p <- plot(bbdml_list[[i]], color = color, size = pointsize) + 
        ggtitle(names(bbdml_list)[i]) + scale_color_manual(values = pal.discrete) + 
        theme_bw() + theme(plot.title = element_text(face = "italic"), 
                           axis.title.y = element_text(face = "bold"), axis.text.x = element_blank())
      plot_list[[i]] <- p
    }
  }
  return(plot_list)
}



extract_bbdml_model <- function(da_analysis,sig.taxa=NA){
  names(da_analysis[["significant_models"]]) <- da_analysis[["significant_taxa"]]
  da_models <- da_analysis[["significant_models"]]
  
  # if provided a list of significant taxa, subset to just those present in list
  # requires a named list where names are the asv sequence
  if(!any(is.na(sig.taxa))){
    da_models <- da_models[names(da_models) %in% names(sig.taxa)]
  }
  
  da_coefs <- da_models %>% map(coef)
  
  x <- purrr::reduce(da_coefs,rbind)
  # remove intercept and mu rows
  keeper_rows <- row.names(x) %>% grep(pattern="phi.|Intercept", value=TRUE,invert = TRUE)
  x <- x[row.names(x) %in% keeper_rows,] %>% as.data.frame()
  x$asv <- names(da_coefs)
  x$taxonomy <- corncob::otu_to_taxonomy(x$asv,ps,level = c('Kingdom','Phylum',"Class",'Order', 'Family',"Genus"))
  
  
  x <- 
  x %>% 
    janitor::clean_names() %>% 
    mutate(est_min = estimate - std_error,
           est_max = estimate + std_error,
           model_factor = keeper_rows)
  return(x)
}

## data ####
ps <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS") 
inoc <- readRDS("./Output/phyloseq_objects/16S_inoculum_samples_clean_phyloseq_object.RDS")

sample_names(inoc)

variables <- names(sample_data(ps))
plant_measures <- c("bud_number","leaf_number","leaf_length","height")
soil_variables <- grep("mean_",variables,value=TRUE)  
predictors <- c("drought","inoculum_site","fire_freq","host")


# PLANT MEASURES ####

## plot ####
plant_measures_plot <- microbiome::meta(ps) %>% 
  mutate(across(all_of(plant_measures),compositions::scale)) %>%
  pivot_longer(all_of(plant_measures),names_to = "plant_measure") %>% 
  dplyr::filter(plant_measure != "bud_number") %>% 
  ggplot(aes(x=inoculum_site,y=value,fill=drought)) +
  geom_boxplot() +
  facet_wrap(~host*plant_measure,scales = 'free',nrow = 2) +
  scale_fill_manual(values=pal.discrete[c(2,5)]) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold')) +
  labs(fill="Drought")
saveRDS(plant_measures_plot,"./Output/figs/Plant_Health_Measures_Plot.RDS")

## model ####
microbiome::meta(ps) %>% 
  mutate(across(all_of(plant_measures),compositions::scale)) %>%
  pivot_longer(all_of(plant_measures),names_to = "plant_measure") %>% 
  dplyr::filter(plant_measure != "bud_number") %>% 
  dplyr::filter(plant_measure != "wilting_scale") %>% 
  # mutate(inoculum_site = inoculum_site %>% factor(levels=c("Sterile",as.character(1:6)))) %>%
  lmer(data = .,
       formula = value ~ drought * host * inoculum_site + (1|block)) %>% 
  broom.mixed::tidy() %>% 
  saveRDS("./Output/Plant_Health_Measures_Overview_Model_table.RDS")


# ALPHA-DIV ESTIMATES ####
alpha <- estimate_richness(ps) %>% 
  select(Observed,Shannon, Simpson)

alpha_df <- microbiome::meta(ps) %>%
  bind_cols(alpha)




# build long-format df with CLR transformed plant measure values
alpha_long_transformed <- 
  alpha_df %>% 
  mutate(across(all_of(plant_measures),function(x){x %>% clr %>% as.numeric})) %>% 
  pivot_longer(all_of(plant_measures),names_to = "plant_measure",values_to = "transformed_value")


## plots ####

# Barplots
# make samples ordered by site
x <- data.frame(site=ps@sam_data$inoculum_site,
                sample=ps@sam_data$sample_name) %>% 
  arrange(site)

sample_order <- 
  x %>% 
  pluck("sample")

break_points <- 
  x %>% group_by(site) %>% 
  summarize(break_point = tail(sample,1)) %>% 
  pluck("break_point")

# melt to data frame for plotting
rank_names(ps)
melted <- 
  ps %>% 
  # subset_taxa(Kingdom == "Bacteria") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt()

# plot bar chart
p <- 
  melted %>% 
  mutate(Sample = factor(Sample,levels=sample_order),
         Genus = ifelse(is.na(Genus),"Undetermined",Genus),
         inoculum_site = factor(inoculum_site, levels= c("Sterile",as.character(1:6)))) %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Phylum)) +
  geom_col() +
  facet_wrap(~inoculum_site,nrow = 1,scales = 'free_x') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5,size=6),
        strip.background = element_blank(),
        strip.text = element_text(face='bold',size=12)) +
  scale_fill_viridis_d() 
p
saveRDS(p,"./Output/figs/16S_Barplot_Phylum.RDS")

### drought ####
p <- 
alpha_long_transformed %>%
  ggplot(aes(x=Shannon,y=transformed_value,color=drought)) +
  geom_point(size=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~host*plant_measure,scales = 'free',nrow = 2) +
  scale_color_manual(values=drought_colors) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold')) +
  labs(x="Shannon diversity",y="CLR-transformed plant measure values",color="Drought")
saveRDS(p,"./Output/figs/16S_Shannon_vs_plant_drought.RDS")
p
ggsave("./Output/figs/16S_Shannon_vs_plant_drought.png",height = 4,width = 10)

### burn freq ####
p <- 
  alpha_long_transformed %>% 
  ggplot(aes(x=Shannon,y=transformed_value,color=ordered(fire_freq,levels=c("0","1","3",'NA')))) +
  geom_point(size=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~host*plant_measure,scales = 'free',nrow = 2) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold')) +
  labs(x="Shannon diversity",y="CLR-transformed plant measure values",color="Fire\nfrequency")
saveRDS(p,"./Output/figs/16S_Shannon_vs_plant_burnfreq.RDS")
p
ggsave("./Output/figs/16S_Shannon_vs_plant_burnfreq.png",height = 4,width = 10)

### inoc site ####
p <- 
alpha_long_transformed %>% 
  ggplot(aes(x=Shannon,y=transformed_value,color=inoculum_site)) +
  geom_point(size=.5) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(~host*plant_measure,scales = 'free',nrow = 2) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold')) +
  labs(x="Shannon diversity",y="CLR-transformed plant measure values",color="Inoculum source") +
  scale_color_viridis_d()
saveRDS(p,"./Output/figs/16S_Shannon_vs_plant_site.RDS")
p
ggsave("./Output/figs/16S_Shannon_vs_plant_site.png",height = 4,width = 10)


## models ####
alpha_long_transformed$overall_transformed_plant_value <- compositions::scale(alpha_long_transformed$transformed_value)
alpha_long_transformed %>% 
  ggplot(aes(x=Shannon,y=overall_transformed_plant_value)) + geom_point()

alpha_mod <- alpha_long_transformed %>% 
  lmer(formula = overall_transformed_plant_value ~ Shannon * host * drought * inoculum_site + (1|block),
       data = .)
summary(alpha_mod)
saveRDS(alpha_mod,"./Output/16S_alpha_diversity_lmermod_model.RDS")
saveRDS(broom.mixed::tidy(alpha_mod),"./Output/16S_alpha_diversity_lmermod_table.RDS")

alpha_df %>% 
  mutate(height = scale(height)) %>% 
  lmer(formula = height ~ Shannon * host * drought * inoculum_site + (1|block), data = .) %>% 
  summary()


# custom plot
# p <- 
#   alpha_df %>% 
#   pivot_longer(c(Observed,Shannon,Simpson),names_to = "Measure") %>% 
#   ggplot(aes(x=inoculum_burn_freq,y=value)) +
#   geom_boxplot() +
#   facet_wrap(~Measure,scales = 'free') +
#   theme_minimal()
# p
# saveRDS(p,"./Output/figs/16S_burn_frequency_and_alpha-div.RDS")
# ggsave("./Output/16S_burn_frequency_and_alpha-div.png",height = 4,width = 6,dpi=300)
# 
# p <- 
# alpha_df %>% 
#   pivot_longer(c(Observed,Shannon,Simpson),names_to = "Measure") %>% 
#   ggplot(aes(x=inoculum_site,y=value)) +
#   geom_boxplot() +
#   facet_wrap(~Measure,scales = 'free') +
#   theme_minimal()
# p
# saveRDS(p,"./Output/figs/16S_inoc_source_and_alpha-div.RDS")
# ggsave("./Output/16S_inoc_source_and_alpha-div.png",height = 4,width = 6,dpi=300)


# BETA-DIV ESTIMATES ####

## ordinations ####

### unifrac nmds ####
ord_nmds <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",distance = "unifrac")

### unifrac rda ####
ord_rda <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "RDA",distance = "unifrac")

### unifrac pcoa ####
ord_pca <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "PCoA",distance = "unifrac")

# build data frame
ord_values <- 
data.frame(
  nmds_1=ord_nmds$points[,1],
  nmds_2=ord_nmds$points[,1],
  rda_1=ord_rda$CA$u[,"PC1"],
  rda_2=ord_rda$CA$u[,"PC2"],
  pcoa_1=ord_pca$vectors[,"Axis.1"],
  pcoa_2=ord_pca$vectors[,"Axis.2"]
)

ord_df <- 
alpha_df %>% 
  select(sample_name,block,drought,fire_freq,host,inoculum_site,Shannon) %>% 
  bind_cols(ord_values) %>% 
  mutate(case_when(inoculum_site == "Sterile" ~ "Sterile"))

ord_df_long <- 
  ord_df %>% 
  pivot_longer(contains("_1"),names_to = "method_x",values_to = "X") %>% 
  pivot_longer(contains("_2"),names_to = "method_y",values_to = "Y") %>% 
  arrange(method_x,method_y,sample_name) %>% 
  mutate(method = method_x %>% str_remove("_1") %>% str_to_upper) %>% 
  select(-starts_with("method_"))


## plots ####
p1 <- 
ord_df_long %>% 
  ggplot(aes(X,Y,color=drought)) +
  geom_point() +
  stat_ellipse() +
  facet_wrap(~host*method,scales = 'free',nrow=2) +
  scale_color_manual(values=drought_colors) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold')) +
  labs(color="Drought",title="Bacterial ordinations")
p1; saveRDS(p1,"./Output/figs/16S_Ordination_Plots_Drought.RDS")
    
p2 <- 
  ord_df_long %>% 
  ggplot(aes(X,Y,color=inoculum_site)) +
  geom_point() +
  stat_ellipse() +
  facet_wrap(~host*method,scales = 'free',nrow=2) +
  scale_color_viridis_d() +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold')) +
  labs(color="Inoculum source",title="Bacterial ordinations")
p2; saveRDS(p2,"./Output/figs/16S_Ordination_Plots_Site.RDS")

p3 <- 
ord_df_long %>% 
  ggplot(aes(X,Y,color=host)) +
  geom_point() +
  stat_ellipse() +
  facet_wrap(~method,scales = 'free',nrow=1) +
  scale_color_manual(values=host_colors) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold'),
        legend.text = element_text(face='italic')) +
  labs(color="Inoculum source",title="Bacterial ordinations") 
p3; saveRDS(p3,"./Output/figs/16S_Ordination_Plots_Host.RDS")

p4 <- 
  ord_df_long %>% 
  ggplot(aes(X,Y,color=ordered(fire_freq,levels=c("0","1","3")))) +
  geom_point() +
  stat_ellipse() +
  facet_wrap(~method,scales = 'free',nrow=1) +
  scale_color_manual(values=fire_colors) +
  theme(strip.text.x = element_text(face='bold.italic'),
        axis.title = element_text(face='bold',size=12),
        legend.title = element_text(face='bold'),
        legend.text = element_text(face='italic')) +
  labs(color="Fire frequency",title="Bacterial ordinations") 
p4; saveRDS(p4,"./Output/figs/16S_Ordination_Plots_Fire.RDS")


## models ####

# make tidy asv table

# have to remove any rows with NA in predictor variables
# sterile samples have NA for fire_freq... Should that be 0 or really NA? (probably NA)
complete.rows <- ord_df %>% dplyr::filter(inoculum_site != "Sterile") %>% row.names()
adonis_df <- ord_df[complete.rows,]

ra_table <- ps %>% 
  subset_samples(sample_names(ps) %in% complete.rows) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% 
  as.data.frame()
  
permanova_results <- 
adonis2(data = adonis_df,
        formula = ra_table ~ adonis_df$host * adonis_df$inoculum_site * adonis_df$drought, strata = adonis_df$block) %>% 
  broom::tidy() %>% 
  mutate(term = term %>% str_remove_all("adonis_df\\$"))
saveRDS(permanova_results,"./Output/16S_Permanova_Table.RDS")


# FIND IMPORTANT TAXA ####

### Do this at the genus level to account for multiple ASVs per species

# make backup of full taxonomy ps
ps_full_taxonomy <- ps
ps_genus <- ps %>% tax_glom("Genus")

# make taxa names...taxonomy
ra_table <- ps_genus %>% 
  subset_samples(sample_names(ps_genus) %in% complete.rows) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% 
  as.data.frame()

taxa_print <- corncob::otu_to_taxonomy(taxa_names(ps_genus),ps_genus)
colnames(ra_table) <- janitor::make_clean_names(taxa_print)
taxa_dictionary <- data.frame(assignment=colnames(ra_table),asv=taxa_names(ps_genus))

# subset relabund table to present taxa
df <- ra_table[,colSums(ra_table)>0]

# make data frame of relative abundances with variables of interest
df <- 
  adonis_df %>% 
  dplyr::select(all_of(predictors)) %>% 
  bind_cols(df)
# convert to logicals for random forest models
df <- 
df %>% 
  mutate(drought=case_when(drought == "D" ~ TRUE,TRUE ~ FALSE),
         grandfir=case_when(host == "Abies grandis" ~ TRUE, TRUE ~ FALSE),
         burned=case_when(fire_freq > 0 ~ TRUE, TRUE ~ FALSE))

## random forest models ####

### drought ####
rf_drought_mod <- 
  df %>% dplyr::select(drought,starts_with("bacteria")) %>%
  ranger::ranger(formula = drought ~ .,
               data=.,importance = 'permutation',num.trees = 999)
saveRDS(rf_drought_mod,"./Output/16S_RF_model_drought.RDS")

### fire ####
rf_fire_mod <- 
  df %>% dplyr::select(burned,starts_with("bacteria")) %>%
  ranger::ranger(formula = burned ~ .,
                 data=.,importance = 'permutation',num.trees = 999)
saveRDS(rf_fire_mod,"./Output/16S_RF_model_fire.RDS")

### host ####
rf_host_mod <- 
  df %>% dplyr::select(grandfir,starts_with("bacteria")) %>%
  ranger::ranger(formula = grandfir ~ .,
                 data=.,importance = 'permutation',num.trees = 999)
saveRDS(rf_host_mod,"./Output/16S_RF_model_host.RDS")

## VIP plots ####
vip::vip(rf_drought_mod) + 
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(title="Top ten most important taxa detected by Random-Forest model",
       subtitle = "Taxa indicative of drought.")
ggsave("./Output/figs/16S_VIP_Plot_drought.png",width = 8,height = 4)

drought_topten <- 
  vip::vi_model(rf_drought_mod) %>% 
  arrange(desc(Importance)) %>% 
  head(10) %>% 
  pluck("Variable")

vip::vip(rf_fire_mod) + 
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(title="Top ten most important taxa detected by Random-Forest model",
       subtitle = "Taxa indicative of burned inoculum.")
ggsave("./Output/figs/16S_VIP_Plot_fire.png",width = 8,height = 4)

fire_topten <- 
  vip::vi_model(rf_fire_mod) %>% 
  arrange(desc(Importance)) %>% 
  head(10) %>% 
  pluck("Variable")

vip::vip(rf_host_mod) + 
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(title="Top ten most important taxa detected by Random-Forest model",
       subtitle = "Taxa indicative of Grand Fir.")
ggsave("./Output/figs/16S_VIP_Plot_host.png",width = 8,height = 4)

# pull together top ten taxa for each model
topten_taxa_rf <- 
vip::vi_model(rf_host_mod) %>% 
  arrange(desc(Importance)) %>% 
  head(10) %>% mutate(model="host") %>% 
  full_join(
    vip::vi_model(rf_fire_mod) %>% 
      arrange(desc(Importance)) %>% 
      head(10) %>% mutate(model="fire")
  ) %>% 
  full_join(
    vip::vi_model(rf_drought_mod) %>% 
      arrange(desc(Importance)) %>% 
      head(10) %>% mutate(model="drought")
  )
    
  
  
# now, find out there relative abundance in each of the model classes (burned vs unburned, etc)

# subset physeq to top ten taxa from each RF model
topten_dictionary <- 
  taxa_dictionary %>% 
  dplyr::filter(assignment %in% unique(topten_taxa_rf$Variable))

ps_genus_topten <- 
ps_genus %>% 
  subset_taxa(taxa_names(ps_genus) %in% topten_dictionary$asv)
# ps_genus_topten %>% tax_table() %>% View
host_topten <- 
topten_dictionary$assignment %in%
  (topten_taxa_rf %>% 
  dplyr::filter(model=="host") %>% 
  pluck("Variable"))

host_topten_asvs <- topten_dictionary[host_topten,"asv"]

fire_topten <- 
  topten_dictionary$assignment %in%
  (topten_taxa_rf %>% 
     dplyr::filter(model=="fire") %>% 
     pluck("Variable"))

fire_topten_asvs <- topten_dictionary[fire_topten,"asv"]

drought_topten <- 
  topten_dictionary$assignment %in%
  (topten_taxa_rf %>% 
     dplyr::filter(model=="drought") %>% 
     pluck("Variable"))

drought_topten_asvs <- topten_dictionary[drought_topten,"asv"]


ps_genus_host_topten <- 
  ps_genus_topten %>% 
  subset_taxa(taxa_names(ps_genus_topten) %in% host_topten_asvs)

ps_genus_fire_topten <- 
  ps_genus_topten %>% 
  subset_taxa(taxa_names(ps_genus_topten) %in% fire_topten_asvs)

ps_genus_drought_topten <- 
  ps_genus_topten %>% 
  subset_taxa(taxa_names(ps_genus_topten) %in% drought_topten_asvs)

## comparison plots ####

# maybe these should be boxplots?
p <- ps_genus_host_topten %>% 
  merge_samples("host") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_topten_relabund() +
  scale_fill_manual(values = host_colors) +
  labs(fill="Host")
saveRDS(p,"./Output/figs/16S_Important_Taxa_Abundances_host.RDS")
# ggsave("./Output/figs/16S_Important_Taxa_Abundances_host.png",width = 14,height = 6)
  
p <- 
ps_genus_fire_topten %>% 
  merge_samples("fire_freq") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_topten_relabund() +
  scale_fill_viridis_d() +
  labs(fill="Fire frequency")
saveRDS(p,"./Output/figs/16S_Important_Taxa_Abundances_fire.RDS")
# ggsave("./Output/figs/16S_Important_Taxa_Abundances_fire.png",width = 14,height = 6)

p <- 
ps_genus_drought_topten %>% 
  merge_samples("drought") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_topten_relabund() +
  scale_fill_manual(values = drought_colors) +
  labs(fill="Drought")
saveRDS(p,"./Output/figs/16S_Important_Taxa_Abundances_drought.RDS")
# ggsave("./Output/figs/16S_Important_Taxa_Abundances_drought.png",width = 8,height = 5)


# HEATMAPS ####
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
  scale_fill_viridis_c(option = "mako") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),axis.text.x = element_blank(),
        strip.text = element_text(face='bold',size=14),
        axis.title.x = element_blank()) +
  labs(fill="Relative\nabundance",y="Most abundant ASVs",title="Comparison of most abundant bacterial ASVs\nin root samples from each inoculum source") +
  facet_wrap(~inoculum_site,scales = 'free_x')

)
saveRDS(p3,"./Output/figs/16S_heatmap_mostabund.RDS")


# Corncob analysis ####

# Clean up ASV names to show taxonomy
ASV_names <- otu_table(ps) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Phylum","Class","Order","Family","Genus"))

genus_names <- otu_table(ps %>% tax_glom("Genus")) %>% colnames()
genus_taxa <- otu_to_taxonomy(genus_names,ps_genus,level = c("Phylum","Class","Order","Family","Genus"))


# set levels of inoculum site so "sterile" is intercept
ps_genus@sam_data$inoculum_site <- 
  ps_genus@sam_data$inoculum_site %>% factor(levels = c("Sterile","1","2","3","4","5","6"))
ps_genus@sam_data$inoculum_burn_freq_uo <- ps_genus@sam_data$inoculum_burn_freq %>% factor(ordered=FALSE,levels = 
                                                                                             c("Sterile","0","1","3"))
ps_genus@sam_data$inoculum_burn_freq_uo[is.na(ps_genus@sam_data$inoculum_burn_freq_uo)] <- "Sterile"
# use raw count data for corncob
da_analysis_host <- differentialTest(formula = ~ host, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_genus,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
da_analysis_inocburnfreq <- differentialTest(formula = ~ inoculum_burn_freq_uo, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_genus,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
da_analysis_drought <- differentialTest(formula = ~ drought, #abundance
                                             phi.formula = ~ 1, #dispersion
                                             formula_null = ~ 1, #mean
                                             phi.formula_null = ~ 1,
                                             test = "Wald", boot = FALSE,
                                             data = ps_genus,
                                             fdr_cutoff = 0.05,
                                             full_output = TRUE)

# TAXA FOUND BY RF AND CORNCOB ####

# find taxa that were identified by corncob and random forest
host_both_significant <- da_analysis_host$significant_taxa[da_analysis_host$significant_taxa %in% host_topten_asvs]
host_both_significant_taxa <- corncob::otu_to_taxonomy(host_both_significant,ps,level = c('Kingdom','Phylum',"Class",'Order', 'Family',"Genus"))

fire_both_significant <- da_analysis_inocburnfreq$significant_taxa[da_analysis_inocburnfreq$significant_taxa %in% fire_topten_asvs]
fire_both_significant_taxa <- corncob::otu_to_taxonomy(fire_both_significant,ps,level = c('Kingdom','Phylum',"Class",'Order', 'Family',"Genus"))

drought_both_significant <- da_analysis_drought$significant_taxa[da_analysis_drought$significant_taxa %in% drought_topten_asvs]
drought_both_significant_taxa <- corncob::otu_to_taxonomy(drought_both_significant,ps,level = c('Kingdom','Phylum',"Class",'Order', 'Family',"Genus"))


# pull out for custom plotting
host_model_df <- extract_bbdml_model(da_analysis = da_analysis_host,sig.taxa = host_both_significant_taxa) 
host_model_df <- host_model_df %>% 
  mutate(model_factor = model_factor %>% str_remove("mu.host"),
         guild = ps@tax_table[host_model_df$asv,"Guild"] %>% as.character())
fire_model_df <- extract_bbdml_model(da_analysis = da_analysis_inocburnfreq,sig.taxa = fire_both_significant_taxa) 
fire_model_df <- fire_model_df %>% 
  mutate(model_factor = model_factor %>% str_replace("mu.inoculum_burn_freq_uo","Burn frequency: "),
         guild = ps@tax_table[fire_model_df$asv,"Guild"] %>% as.character())
drought_model_df <- extract_bbdml_model(da_analysis = da_analysis_drought,sig.taxa = drought_both_significant_taxa) 
drought_model_df <- drought_model_df %>% 
  mutate(model_factor = model_factor %>% str_replace("mu.droughtND","No drought"),
         guild = ps@tax_table[drought_model_df$asv,"Guild"] %>% as.character())

## Overview plots ####
p <- host_model_df %>%
  ggplot(aes(y=taxonomy,x=estimate)) +
  geom_vline(xintercept=0,linetype=2,color='gray') +
  geom_errorbarh(aes(xmin=est_min,xmax=est_max,height=.3)) +
  labs(subtitle = "Intercept: A. grandis",color="Guild") +
  facet_wrap(~model_factor) +
  theme(strip.text = element_text(face='bold.italic',size=12),
        axis.text.y = element_text(face="bold.italic")) +
  scale_color_manual(values = "brown4")

saveRDS(p,"./Output/figs/16S_DiffAbund_Overview_host.RDS")
# ggsave("./Output/figs/16S_DiffAbund_Overview_host.png",dpi=300,height = 4,width = 12)


p <- fire_model_df %>%
  ggplot(aes(y=taxonomy,x=estimate)) +
  geom_vline(xintercept=0,linetype=2,color='gray') +
  geom_errorbarh(aes(xmin=est_min,xmax=est_max,height=.3)) +
  labs(subtitle = "Intercept: A. grandis",color="Guild") +
  facet_wrap(~model_factor) +
  theme(strip.text = element_text(face='bold.italic',size=12),
        axis.text.y = element_text(face="bold.italic")) +
  scale_color_manual(values = "brown4")
saveRDS(p,"./Output/figs/16S_DiffAbund_Overview_fire.RDS")
# ggsave("./Output/figs/16S_DiffAbund_Overview_fire.png",dpi=300,height = 4,width = 14)


p <- drought_model_df %>% 
  ggplot(aes(y=taxonomy,x=estimate)) +
  geom_vline(xintercept=0,linetype=2,color='gray') +
  geom_errorbarh(aes(xmin=est_min,xmax=est_max,height=.3)) +
  labs(subtitle = "Intercept: Drought") +
  facet_wrap(~model_factor) +
  theme(strip.text = element_text(face='bold.italic',size=12),
        axis.text.y = element_text(face="bold.italic")) +
  scale_color_manual(values = "brown4")

saveRDS(p,"./Output/figs/16S_DiffAbund_Overview_drought.RDS")
# ggsave("./Output/figs/16S_DiffAbund_Overview_drought.png",dpi=300,height = 4,width = 12)


########## Stopped here Jan 24th, 2024 ...

# Might want to get rid of perfectly discriminant taxa (found only in one condition)  ???
# (pull top 20 from RF, match those found with corncob, remove those that are discriminant, then use that list) ???

## bbdml plots ####
host_multi_bbdml <- multi_bbdml(da_analysis_host,ps,"host",'host')
host_bbdml_plots <- plot_multi_bbdml(host_multi_bbdml,color='host')
saveRDS(host_bbdml_plots,"./Output/figs/host_bbdml_plots.RDS")

fire_multi_bbdml <- multi_bbdml(da_analysis_fire,ps,"inoculum_burn_freq",'inoculum_burn_freq')
fire_bbdml_plots <- plot_multi_bbdml(fire_multi_bbdml,color='inoculum_burn_freq')
saveRDS(fire_bbdml_plots,"./Output/figs/fire_bbdml_plots.RDS")

drought_multi_bbdml <- multi_bbdml(da_analysis_drought,ps,"drought",'drought')
drought_bbdml_plots <- plot_multi_bbdml(drought_multi_bbdml,color='drought')
saveRDS(drought_bbdml_plots,"./Output/figs/drought_bbdml_plots.RDS")


#







p4 <- 
plot(da_analysis_host) +
  theme(legend.position = 'none') 
saveRDS(p4,"./Output/16S_diffabund_inocsite.RDS")

# quick look at Rhizobium and Phenylobacterium and Neorhizobium

rhizobium_asv <- 
  ps_genus@tax_table[which(ps_genus@tax_table[,6] == "Rhizobium"),] %>% row.names()
phenylobacterium_asv <- 
  ps_genus@tax_table[which(ps_genus@tax_table[,6] == "Phenylobacterium"),] %>% row.names()
neorhizobium_asv <- 
  ps_genus@tax_table[which(ps_genus@tax_table[,6] == "Neorhizobium"),] %>% row.names()

sig_taxa <- da_analysis_inocsite$significant_taxa %>% otu_to_taxonomy(data=ps_genus)
bbdml_obj <- multi_bbdml(da_analysis_inocsite,
                         ps_object = ps_genus,
                         mu_predictor = "inoculum_site",
                         phi_predictor = "inoculum_site",
                         taxlevels = 4:6)

names(bbdml_obj)

# plot 2 significant taxa, for now (rhizobium and phenylobacterium)
plot_multi_bbdml(bbdml_obj[c(2,9)],color = "inoculum_site",pointsize = 3)  

bbdml_plot_1
bbdml_plot_2
saveRDS(bbdml_plot_1,"./Output/16S_bbdml_plot_1.RDS")
saveRDS(bbdml_plot_2,"./Output/16S_bbdml_plot_2.RDS")
