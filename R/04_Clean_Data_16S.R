# -----------------------------------------------------------------------------#
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.46.0
#                     ShortRead v 1.60.0
#                     Biostrings v 2.70.1
#                     adegenet v 2.1.10
#                     readxl v 1.4.3
#                     janitor v 2.2.0
#                     microbiome v 1.24.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#        Remove non-bacteria, chloroplast and mitochondrial sequences           #
#                                                                               #
#################################################################################

# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(adegenet); packageVersion("adegenet")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")
library(microbiome); packageVersion("microbiome")

# seed
set.seed(666)

# DATA ####
ps <- readRDS("./Output/16S_ps_not-cleaned_w_tree.RDS") # change to non-phylogeny stuff
sample_names(ps)

# POSITIVE CONTROL ####
# check positive control efficacy...

pos <- ps %>% 
  subset_samples(sample_names(ps) == "P1")
y <- tax_table(pos)[which(taxa_sums(pos) > 0),5:7] %>% as.data.frame()

x <- as.data.frame(otu_table(pos))[,which(taxa_sums(pos) > 0)]
y$abund <- as.numeric(x[1,])

y %>% group_by(Genus) %>% summarize(total=sum(abund)) %>% arrange(desc(total)) %>% mutate(relabund=total/sum(total)) %>% 
  mutate(Genus=factor(Genus,levels=pluck(.,"Genus"))) %>% 
  ggplot(aes(x=Genus,y=relabund)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,hjust=1)) +
  labs(y="Relative abundance",title="Bacterial positive control sample")
ggsave("./Output/figs/16S_Bacterial_positive_control_taxonomy.png")

# Remove postive control
ps <- ps %>% subset_samples(sample_names(ps) != "P1")

# SEPARATE INOCULUM SAMPLES ####
# back-up full ps object
ps_full <- ps
# pull out inoculum reads
inoc <- subset_samples(ps, !is.na(other_frompreviouscolumn))
# remove non-bacteria
inoc <- subset_taxa(inoc, Kingdom == "Bacteria")
inoc <- subset_taxa(inoc,Class != "Chloroplast")
inoc <- subset_taxa(inoc, taxa_sums(inoc) > 0)
inoc <- subset_samples(inoc, sample_sums(inoc) > 0)
sample_names(inoc)

# REMOVE NON-BACTERIA FROM FULL PS ####
# Remove inoculum samples first
ps <- subset_samples(ps, is.na(other_frompreviouscolumn))
ps <- subset_taxa(ps, Kingdom == "Bacteria")
ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)
sample_names(ps)

# add variable for initial/final
ps_full@sam_data$community <- 
factor(ifelse(is.na(ps_full@sam_data$inoculum_site),
       "inoculum","final"),levels=c("inoculum","final"))


# Pull out full ps objects for each inoculum
site1.inoc.full <- ps_full %>% subset_samples(other_frompreviouscolumn == "Site1" | inoculum_site == "1")
site2.inoc.full <- ps_full %>% subset_samples(other_frompreviouscolumn == "Site2" | inoculum_site == "2")
site3.inoc.full <- ps_full %>% subset_samples(other_frompreviouscolumn == "Site3" | inoculum_site == "3")
site4.inoc.full <- ps_full %>% subset_samples(other_frompreviouscolumn == "Site4" | inoculum_site == "4")
site5.inoc.full <- ps_full %>% subset_samples(other_frompreviouscolumn == "Site5" | inoculum_site == "5")
site6.inoc.full <- ps_full %>% subset_samples(other_frompreviouscolumn == "Site6" | inoculum_site == "6")

# remove empty taxa
site1.inoc <- 
  site1.inoc.full %>% 
  subset_taxa(taxa_sums(site1.inoc.full) > 0)
site2.inoc <- 
  site2.inoc.full %>% 
  subset_taxa(taxa_sums(site2.inoc.full) > 0)
site3.inoc <- 
  site3.inoc.full %>% 
  subset_taxa(taxa_sums(site3.inoc.full) > 0)
site4.inoc <- 
  site4.inoc.full %>% 
  subset_taxa(taxa_sums(site4.inoc.full) > 0)
site5.inoc <- 
  site5.inoc.full %>% 
  subset_taxa(taxa_sums(site5.inoc.full) > 0)
site6.inoc <- 
  site6.inoc.full %>% 
  subset_taxa(taxa_sums(site6.inoc.full) > 0)

# SAVE CLEAN PHYLOSEQ OBJECTS ####
saveRDS(site1.inoc.full,"Output/phyloseq_objects/16S_site1.inoc.full.RDS")
saveRDS(site2.inoc.full,"Output/phyloseq_objects/16S_site2.inoc.full.RDS")
saveRDS(site3.inoc.full,"Output/phyloseq_objects/16S_site3.inoc.full.RDS")
saveRDS(site4.inoc.full,"Output/phyloseq_objects/16S_site4.inoc.full.RDS")
saveRDS(site5.inoc.full,"Output/phyloseq_objects/16S_site5.inoc.full.RDS")
saveRDS(site6.inoc.full,"Output/phyloseq_objects/16S_site6.inoc.full.RDS")
saveRDS(site1.inoc,"Output/phyloseq_objects/16S_site1.inoc.RDS")
saveRDS(site2.inoc,"Output/phyloseq_objects/16S_site2.inoc.RDS")
saveRDS(site3.inoc,"Output/phyloseq_objects/16S_site3.inoc.RDS")
saveRDS(site4.inoc,"Output/phyloseq_objects/16S_site4.inoc.RDS")
saveRDS(site5.inoc,"Output/phyloseq_objects/16S_site5.inoc.RDS")
saveRDS(site6.inoc,"Output/phyloseq_objects/16S_site6.inoc.RDS")
# full objects
saveRDS(ps, file = "./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")
saveRDS(inoc, file = "./Output/phyloseq_objects/16S_inoculum_samples_clean_phyloseq_object.RDS")





# 
# # clean up newly missing taxa ... after performing MRM using ALL taxa
# ps <- subset_taxa(ps,taxa_sums(ps) > 0)
# # inoc %>% microbiome::meta() %>% View
# 
# # Clean metadata ####
# microbiome::meta(ps) %>% apply(.,2,class)
# 
# # clean variable classes
# ps@sam_data$fire_freq <- as.numeric(ps@sam_data$fire_freq)
# ps@sam_data$wilting_scale <- as.numeric(ps@sam_data$wilting_scale)
# ps@sam_data$bud_number <- as.numeric(ps@sam_data$bud_number)
# ps@sam_data$leaf_number <- as.numeric(ps@sam_data$leaf_number)
# ps@sam_data$leaf_length <- as.numeric(ps@sam_data$leaf_length)
# ps@sam_data$height <- as.numeric(ps@sam_data$height)
# ps@sam_data$nodule_num <- as.numeric(ps@sam_data$nodule_num)
# ps@sam_data$shoot_dm <- as.numeric(ps@sam_data$shoot_dm)
# ps@sam_data$final_root_dm <- as.numeric(ps@sam_data$final_root_dm)
# 
# # add inoculum "site" and "burn_freq"
# ps@sam_data$inoculum_site <- ps@sam_data$inoculum %>% str_remove(".Burn_W")
# ps@sam_data$inoculum_burn_freq <- ps@sam_data$inoculum %>% str_remove("Burn_.*")
# ps_full@sam_data$inoculum_site <- ps_full@sam_data$inoculum %>% str_remove(".Burn_W")
# ps_full@sam_data$inoculum_burn_freq <- ps_full@sam_data$inoculum %>% str_remove("Burn_.*")
# 
# 
# # Add Soil Chemistry ####
# chem <- read_csv("./Data/SoilChem.csv")
# 
# # make id column match ps
# chem$inoculum <- 
# chem$Field %>% 
#   str_remove(pattern = "\\.[1]") %>% 
#   str_remove(pattern = "\\.[2]") %>% 
#   str_remove(pattern = "\\.[3]") %>% 
#   str_remove(pattern = "Unmixed") %>% 
#   str_remove(pattern = "Mixed") %>% 
#   str_replace(pattern = "-",replacement = "_") %>% 
#   str_remove(pattern = "_$")
# # find mean values of all variables for each inoc site
# chem <- 
# chem %>% 
#   group_by(inoculum) %>% 
#   summarize(mean_NO3 = mean(NO3,na.rm=TRUE),
#             mean_NH4 = mean(NH4,na.rm=TRUE),
#             mean_P = mean(P,na.rm=TRUE),
#             mean_K = mean(K,na.rm=TRUE),
#             mean_SO4_S = mean(`SO4-S`,na.rm=TRUE),
#             mean_B = mean(B,na.rm=TRUE),
#             mean_OM = mean(OM,na.rm=TRUE),
#             mean_pH = mean(pH,na.rm=TRUE),
#             mean_EC = mean(EC,na.rm=TRUE),
#             mean_Zn = mean(Zn,na.rm=TRUE),
#             mean_Mn = mean(Mn,na.rm=TRUE),
#             mean_Cu = mean(Cu,na.rm=TRUE),
#             mean_Fe = mean(Fe,na.rm=TRUE),
#             mean_Ca = mean(Ca,na.rm=TRUE),
#             mean_Mg = mean(Mg,na.rm=TRUE),
#             mean_Na = mean(Na,na.rm=TRUE),
#             mean_TKN = mean(TKN,na.rm=TRUE),
#             mean_C_N_ratio = mean(`C:N`,na.rm=TRUE),
#             mean_Sand = mean(Sand,na.rm=TRUE),
#             mean_Silt = mean(Silt,na.rm=TRUE),
#             mean_Clay = mean(Clay,na.rm=TRUE))
# 
# 
# # join soil chemistry data with sample metadata
# meta <- microbiome::meta(ps) %>% full_join(chem,by="inoculum") %>% sample_data()
# sample_names(meta) <- meta$sample_name
# ps@sam_data <- sample_data(meta)
# 
# meta <- microbiome::meta(ps_full) %>% full_join(chem,by="inoculum") %>% sample_data()
# sample_names(meta) <- meta$sample_name
# ps_full@sam_data <- sample_data(meta)  
# 
# # fix "burn frequency" variable to ordered factor
# ps@sam_data$inoculum_burn_freq <- 
#   ps@sam_data$inoculum_burn_freq %>% as.numeric() %>% ordered()
# ps_full@sam_data$inoculum_burn_freq <- 
#   ps_full@sam_data$inoculum_burn_freq %>% as.numeric() %>% ordered()
# 
# # pull out full ps for each inoculum site (recipients)
# for(i in as.character(1:6)){
#   x <- ps_full %>% subset_samples(inoculum_site == i)
#   saveRDS(x,file.path(paste0("./Output/inoculum_ps_from_site_",i,".RDS")))
#   assign(paste0("inoculum_ps_site_",i),x,envir = .GlobalEnv)
# }
# # site1 <- ps_full %>% subset_samples(inoculum_site == "1")
# # site2 <- ps_full %>% subset_samples(inoculum_site == "2")
# # site3 <- ps_full %>% subset_samples(inoculum_site == "3")
# # site4 <- ps_full %>% subset_samples(inoculum_site == "4")
# # site5 <- ps_full %>% subset_samples(inoculum_site == "5")
# # site6 <- ps_full %>% subset_samples(inoculum_site == "6")
# 
# # subset to just taxa found in either source or donor
# # test1 <- 
# #   site1 %>% subset_taxa(
# #   taxa_names(ps_full) %in% taxa_names(ps_full)[which(taxa_sums(site1) > 0 | taxa_sums(site1.inoc) > 0)])
# # test1.inoc <- 
# #   site1.inoc %>% subset_taxa(
# #   taxa_names(ps_full) %in% taxa_names(ps_full)[which(taxa_sums(site1) > 0 | taxa_sums(site1.inoc) > 0)])
# # 
# # # compare site1.inoc community with site2 recipient community
# # test2 <- 
# #   site1 %>% subset_taxa(
# #     taxa_names(ps_full) %in% taxa_names(ps_full)[which(taxa_sums(site2) > 0 | taxa_sums(site1.inoc) > 0)])
# # test2.inoc <- 
# #   site1.inoc %>% subset_taxa(
# #     taxa_names(ps_full) %in% taxa_names(ps_full)[which(taxa_sums(site2) > 0 | taxa_sums(site1.inoc) > 0)])
# # 
# # dist.site1.inoc <- dist(t(otu_table(test1.inoc))) 
# # dist.site1 <- dist(t(otu_table(test2)))
# # dist.site1.inoc <- dist(t(otu_table(test2.inoc))) 
# # 
# # # MRM
# # mrm.test <- ecodist::MRM(dist.site1 ~ dist.site1.inoc)
# # as.data.frame(mrm.test)
# 
# 
# 
# Save DNA sequences apart from rownames (from subsetted ps object)
# seqs <- taxa_names(ps)
# seqs <- DNAStringSet(seqs)
# saveRDS(seqs,"./Output/16S_ASV_reference_sequences.RDS")
# 
# seqs <- taxa_names(inoc)
# seqs <- DNAStringSet(seqs)
# saveRDS(seqs,"./Output/16S_Inoculum_ASV_reference_sequences.RDS")
# 
# saveRDS(ps, file = "./Output/16S_clean_phyloseq_object.RDS")
# saveRDS(inoc, file = "./Output/16S_inoculum_samples_clean_phyloseq_object.RDS")
# 
