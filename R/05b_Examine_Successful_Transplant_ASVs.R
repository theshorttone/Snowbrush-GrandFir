# -----------------------------------------------------------------------------#
# Cleaning up the phyloseq object (fungal)
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.46.0
#                     corncob v 0.4.1
# -----------------------------------------------------------------------------#

# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(corncob); packageVersion("corncob")
library(ranger); packageVersion("ranger")
library(vip); packageVersion("vip")
ra <- function(x){x/sum(x)}
'%ni%' <- Negate('%in%')

# DATA ####

# Load physeq objects 
bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS")
fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")
amf <- readRDS("./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")
bact_full <- bact
fung_full <- fung
amf_full <- amf

# load successful asv lists
asvs_bact <- readRDS("./Output/Successful_ASVs_16S.RDS")
asvs_fung <- readRDS("./Output/Successful_ASVs_ITS.RDS")
asvs_amf <- readRDS("./Output/Successful_ASVs_18S.RDS")


# subset physeqs to those taxa only
bact <- bact %>% 
  subset_taxa(taxa_names(bact) %in% asvs_bact)
fung <- fung %>% 
  subset_taxa(taxa_names(fung) %in% asvs_fung)
amf <- amf %>% 
  subset_taxa(taxa_names(amf) %in% asvs_amf)
bact@sam_data %>% names

# diffabund test
bact_da <- corncob::differentialTest(formula = ~ species + drought,
                          phi.formula = ~ 1, #dispersion
                          formula_null = ~ 1, #mean
                          phi.formula_null = ~ 1,
                          test = "Wald", boot = FALSE,
                          data = bact,
                          fdr_cutoff = 0.05,
                          full_output = TRUE)

fung_da <- corncob::differentialTest(formula = ~ species + drought,
                                     phi.formula = ~ 1, #dispersion
                                     formula_null = ~ 1, #mean
                                     phi.formula_null = ~ 1,
                                     test = "Wald", boot = FALSE,
                                     data = fung,
                                     fdr_cutoff = 0.05,
                                     full_output = TRUE)

amf_da <- corncob::differentialTest(formula = ~ drought,
                                     phi.formula = ~ 1, #dispersion
                                     formula_null = ~ 1, #mean
                                     phi.formula_null = ~ 1,
                                     test = "Wald", boot = FALSE,
                                     data = amf,
                                     fdr_cutoff = 0.05,
                                     full_output = TRUE)
if(length(bact_da$significant_taxa) > 0){
  plot(bact_da)
}

if(length(fung_da$significant_taxa) > 0){
  plot(fung_da)
}

if(length(amf_da$significant_taxa) > 0){
  plot(amf_da)
}

# nothing in exp. design was predictive of successful taxa :(


# Inoc properties predictive of success?

# reload
# Load physeq objects 
# bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS")
# fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")
# amf <- readRDS("./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")


# find properties of inoculum communities
# find 'rate of success' for each final sample
# model whether those inoc properties predict rate of success


# rate of success
bact_mat <- bact@otu_table %>% as("matrix")
bact_mat[bact_mat > 0] <- 1
bact_success_rate <- rowSums(bact_mat) / ntaxa(readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS"))
bact@sam_data$success_rate <- bact_success_rate

fung_mat <- fung@otu_table %>% as("matrix")
fung_mat[fung_mat > 0] <- 1
fung_success_rate <- rowSums(fung_mat) / ntaxa(readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS"))
fung@sam_data$success_rate <- fung_success_rate

amf_mat <- amf@otu_table %>% as("matrix")
amf_mat[amf_mat > 0] <- 1
amf_success_rate <- rowSums(amf_mat) / ntaxa(readRDS("./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS"))
amf@sam_data$success_rate <- amf_success_rate


# bacteria
# pull success rates for samples from each inoculum
bact_inoc_success <- list()
for(i in bact@sam_data$inoculum_site %>% unique %>% sort){
  bact_inoc_success[[i]] <- 
  bact %>% 
    subset_samples(inoculum_site == "1") %>% 
    microbiome::meta() %>% 
    pluck("success_rate")
}


# import inoculum data
site1.inoc.full <- readRDS("Output/phyloseq_objects/16S_site1.inoc.full.RDS")
site1.inoc.full <- subset_samples(site1.inoc.full,sample_names(site1.inoc.full) %ni% grep("^[A,B,C]",sample_names(site1.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site1.inoc.full <- site1.inoc.full %>% subset_taxa(taxa_sums(site1.inoc.full) > 0)

site2.inoc.full <- readRDS("Output/phyloseq_objects/16S_site2.inoc.full.RDS")
site2.inoc.full <- subset_samples(site2.inoc.full,sample_names(site2.inoc.full) %ni% grep("^[A,B,C]",sample_names(site2.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site2.inoc.full <- site2.inoc.full %>% subset_taxa(taxa_sums(site2.inoc.full) > 0)

site3.inoc.full <- readRDS("Output/phyloseq_objects/16S_site3.inoc.full.RDS")
site3.inoc.full <- subset_samples(site3.inoc.full,sample_names(site3.inoc.full) %ni% grep("^[A,B,C]",sample_names(site3.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site3.inoc.full <- site3.inoc.full %>% subset_taxa(taxa_sums(site3.inoc.full) > 0)

site4.inoc.full <- readRDS("Output/phyloseq_objects/16S_site4.inoc.full.RDS")
site4.inoc.full <- subset_samples(site4.inoc.full,sample_names(site4.inoc.full) %ni% grep("^[A,B,C]",sample_names(site4.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site4.inoc.full <- site4.inoc.full %>% subset_taxa(taxa_sums(site4.inoc.full) > 0)

site5.inoc.full <- readRDS("Output/phyloseq_objects/16S_site5.inoc.full.RDS")
site5.inoc.full <- subset_samples(site5.inoc.full,sample_names(site5.inoc.full) %ni% grep("^[A,B,C]",sample_names(site5.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site5.inoc.full <- site5.inoc.full %>% subset_taxa(taxa_sums(site5.inoc.full) > 0)

site6.inoc.full <- readRDS("Output/phyloseq_objects/16S_site6.inoc.full.RDS")
site6.inoc.full <- subset_samples(site6.inoc.full,sample_names(site6.inoc.full) %ni% grep("^[A,B,C]",sample_names(site6.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site6.inoc.full <- site6.inoc.full %>% subset_taxa(taxa_sums(site6.inoc.full) > 0)
sample_names(site1.inoc.full)



# randomforest models
bact1_rf <- data.frame(rate=bact_inoc_success[["1"]]) %>% 
  bind_cols(site1.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

bact2_rf <- data.frame(rate=bact_inoc_success[["2"]]) %>% 
  bind_cols(site2.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

bact3_rf <- data.frame(rate=bact_inoc_success[["3"]]) %>% 
  bind_cols(site3.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

bact4_rf <- data.frame(rate=bact_inoc_success[["4"]]) %>% 
  bind_cols(site4.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

bact5_rf <- data.frame(rate=bact_inoc_success[["5"]]) %>% 
  bind_cols(site5.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

bact6_rf <- data.frame(rate=bact_inoc_success[["6"]]) %>% 
  bind_cols(site6.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

bact1_rf$r.squared
bact2_rf$r.squared
bact3_rf$r.squared
bact4_rf$r.squared
bact5_rf$r.squared
bact6_rf$r.squared

# fungi
# pull success rates for samples from each inoculum
fung_inoc_success <- list()
for(i in fung@sam_data$inoculum_site %>% unique %>% sort){
  fung_inoc_success[[i]] <- 
    fung %>% 
    subset_samples(inoculum_site == "1") %>% 
    microbiome::meta() %>% 
    pluck("success_rate")
}


# import inoculum data
site1.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site1.inoc.full.RDS")
site1.inoc.full <- subset_samples(site1.inoc.full,sample_names(site1.inoc.full) %ni% grep("^F-[A,B,C]",sample_names(site1.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site1.inoc.full <- site1.inoc.full %>% subset_taxa(taxa_sums(site1.inoc.full) > 0)

site2.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site2.inoc.full.RDS")
site2.inoc.full <- subset_samples(site2.inoc.full,sample_names(site2.inoc.full) %ni% grep("^F-[A,B,C]",sample_names(site2.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site2.inoc.full <- site2.inoc.full %>% subset_taxa(taxa_sums(site2.inoc.full) > 0)

site3.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site3.inoc.full.RDS")
site3.inoc.full <- subset_samples(site3.inoc.full,sample_names(site3.inoc.full) %ni% grep("^F-[A,B,C]",sample_names(site3.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site3.inoc.full <- site3.inoc.full %>% subset_taxa(taxa_sums(site3.inoc.full) > 0)

site4.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site4.inoc.full.RDS")
site4.inoc.full <- subset_samples(site4.inoc.full,sample_names(site4.inoc.full) %ni% grep("^F-[A,B,C]",sample_names(site4.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site4.inoc.full <- site4.inoc.full %>% subset_taxa(taxa_sums(site4.inoc.full) > 0)

site5.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site5.inoc.full.RDS")
site5.inoc.full <- subset_samples(site5.inoc.full,sample_names(site5.inoc.full) %ni% grep("^F-[A,B,C]",sample_names(site5.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site5.inoc.full <- site5.inoc.full %>% subset_taxa(taxa_sums(site5.inoc.full) > 0)

site6.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site6.inoc.full.RDS")
site6.inoc.full <- subset_samples(site6.inoc.full,sample_names(site6.inoc.full) %ni% grep("^F-[A,B,C]",sample_names(site6.inoc.full),value=TRUE)) %>% transform_sample_counts(ra)
site6.inoc.full <- site6.inoc.full %>% subset_taxa(taxa_sums(site6.inoc.full) > 0)
sample_names(site1.inoc.full)



# randomforest models
fung1_rf <- data.frame(rate=fung_inoc_success[["1"]]) %>% 
  bind_cols(site1.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

fung2_rf <- data.frame(rate=fung_inoc_success[["2"]]) %>% 
  bind_cols(site2.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

fung3_rf <- data.frame(rate=fung_inoc_success[["3"]]) %>% 
  bind_cols(site3.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

fung4_rf <- data.frame(rate=fung_inoc_success[["4"]]) %>% 
  bind_cols(site4.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

fung5_rf <- data.frame(rate=fung_inoc_success[["5"]]) %>% 
  bind_cols(site5.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

fung6_rf <- data.frame(rate=fung_inoc_success[["6"]]) %>% 
  bind_cols(site6.inoc.full@otu_table %>% as('matrix') %>% as.data.frame) %>% 
  ranger(formula = rate ~ .,importance = 'permutation')

fung1_rf$r.squared
fung2_rf$r.squared
fung3_rf$r.squared
fung4_rf$r.squared
fung5_rf$r.squared
fung6_rf$r.squared


fung1_important_taxa <- corncob::otu_to_taxonomy(OTU = vi(fung1_rf)[1:10,"Variable"] %>% unlist,data = site1.inoc.full)
site1.inoc.full %>% transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(site1.inoc.full) %in% names(fung1_important_taxa)) %>% 
  sample_sums()


# frankia vs plant growth
frankia <- bact_full %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(Genus == "Frankia") %>% 
  subset_samples(species == "Snowbrush")

rhizobium <- bact_full %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(Genus == "Rhizobium") %>% 
  subset_samples(species == "Snowbrush") 
  
bact_mutualists <- 
data.frame(leaf_number=frankia@sam_data$leaf_number,
           frankia_ra = sample_sums(frankia),
           rhizobium_ra = sample_sums(rhizobium),
           block = frankia@sam_data$block,
           nodules = frankia@sam_data$nodule_num,
           drought = frankia@sam_data$drought)

bact_mutualists %>% 
  pivot_longer(ends_with("ra")) %>% 
  ggplot(aes(x=value,y=leaf_number,color=name)) +
  geom_point() +
  geom_smooth(method = 'lm')

bact_mutualists %>% 
  # pivot_longer(ends_with("ra")) %>% 
  ggplot(aes(x=nodules,y=frankia_ra)) +
  geom_point()

bact_mutualists %>% 
  lmerTest::lmer(data=.,
                 formula= leaf_number ~ (frankia_ra + rhizobium_ra) * drought + (1|block)) %>% 
  summary
  

frankia <- bact_full %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(Genus == "Frankia") %>% 
  subset_samples(species == "GrandFir")

rhizobium <- bact_full %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(Genus == "Rhizobium") %>% 
  subset_samples(species == "GrandFir") 

bact_mutualists <- 
  data.frame(leaf_number=frankia@sam_data$leaf_number,
             frankia_ra = sample_sums(frankia),
             rhizobium_ra = sample_sums(rhizobium),
             block = frankia@sam_data$block,
             nodules = frankia@sam_data$nodule_num,
             drought = frankia@sam_data$drought)

bact_mutualists %>% 
  pivot_longer(ends_with("ra")) %>% 
  ggplot(aes(x=value,y=leaf_number,color=name)) +
  geom_point() +
  geom_smooth(method = 'lm')

bact_mutualists %>% 
  # pivot_longer(ends_with("ra")) %>% 
  ggplot(aes(x=nodules,y=frankia_ra)) +
  geom_point()

bact_mutualists %>% 
  lmerTest::lmer(data=.,
                 formula= leaf_number ~ (frankia_ra + rhizobium_ra) * drought + (1|block)) %>% 
  summary
