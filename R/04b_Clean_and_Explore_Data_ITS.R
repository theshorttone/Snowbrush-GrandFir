# -----------------------------------------------------------------------------#
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.40.0
#                     ShortRead v 1.54.0
#                     Biostrings v 2.64.0
#                     adegenet v 2 .1 10
#                     readxl v 1.4.1
#                     janitor v 2.1.0
#                     microbiome v 1.20.0
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
set.seed(789)

# DATA ####
ps <- readRDS("./Output/ITS_ps_not-cleaned.RDS") # change to non-phylogeny stuff

# Remove Non-Fungi
ps <- ps %>% subset_taxa(Kingdom == "k__Fungi")
ps@tax_table[,1] %>% unique

# check positive control efficacy...
pos <- ps %>% 
  subset_samples(control == "Positive")
y <- tax_table(pos)[which(taxa_sums(pos) > 0),5:7] %>% as.data.frame()


x <- as.data.frame(otu_table(pos))[,which(taxa_sums(pos) > 0)]
y$abund <- as.numeric(x[1,])

y %>% group_by(Genus) %>% summarize(total=sum(abund)) %>% arrange(desc(total)) %>% mutate(relabund=total/sum(total)) %>% 
  mutate(Genus=factor(Genus,levels=pluck(.,"Genus"))) %>% 
  ggplot(aes(x=Genus,y=relabund)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,hjust=1)) +
  labs(y="Relative abundance",title="Fungal positive control sample")
ggsave("./Output/ITS_Fungal_positive_control_taxonomy.png",height = 4,width = 6)


# SEPARATE INOCULUM SAMPLES ####

inoc <- ps %>% 
  subset_samples(other_frompreviouscolumn %in% c("Site1","Site2","Site3","Site4","Site5","Site6"))
inoc <- subset_taxa(inoc,taxa_sums(inoc) > 0)

ps <- ps %>% 
  subset_samples(block != "NA")
ps <- subset_taxa(ps,taxa_sums(ps) > 0)

# Clean metadata ####
microbiome::meta(ps) %>% apply(.,2,class)

# clean variable classes
ps@sam_data$fire_freq <- as.numeric(ps@sam_data$fire_freq)
ps@sam_data$wilting_scale <- as.numeric(ps@sam_data$wilting_scale)
ps@sam_data$bud_number <- as.numeric(ps@sam_data$bud_number)
ps@sam_data$leaf_number <- as.numeric(ps@sam_data$leaf_number)
ps@sam_data$leaf_length <- as.numeric(ps@sam_data$leaf_length)
ps@sam_data$height <- as.numeric(ps@sam_data$height)
ps@sam_data$nodule_num <- as.numeric(ps@sam_data$nodule_num)
ps@sam_data$shoot_dm <- as.numeric(ps@sam_data$shoot_dm)
ps@sam_data$final_root_dm <- as.numeric(ps@sam_data$final_root_dm)

# add inoculum "site" and "burn_freq"
ps@sam_data$inoculum_site <- ps@sam_data$inoculum %>% str_remove(".Burn_W")
ps@sam_data$inoculum_burn_freq <- ps@sam_data$inoculum %>% str_remove("Burn_.*")

# Add Soil Chemistry ####
chem <- read_csv("./Data/SoilChem.csv")

# make id column match ps
chem$inoculum <- 
  chem$Field %>% 
  str_remove(pattern = "\\.[1]") %>% 
  str_remove(pattern = "\\.[2]") %>% 
  str_remove(pattern = "\\.[3]") %>% 
  str_remove(pattern = "Unmixed") %>% 
  str_remove(pattern = "Mixed") %>% 
  str_replace(pattern = "-",replacement = "_") %>% 
  str_remove(pattern = "_$")
# find mean values of all variables for each inoc site
chem <- 
  chem %>% 
  group_by(inoculum) %>% 
  summarize(mean_NO3 = mean(NO3,na.rm=TRUE),
            mean_NH4 = mean(NH4,na.rm=TRUE),
            mean_P = mean(P,na.rm=TRUE),
            mean_K = mean(K,na.rm=TRUE),
            mean_SO4_S = mean(`SO4-S`,na.rm=TRUE),
            mean_B = mean(B,na.rm=TRUE),
            mean_OM = mean(OM,na.rm=TRUE),
            mean_pH = mean(pH,na.rm=TRUE),
            mean_EC = mean(EC,na.rm=TRUE),
            mean_Zn = mean(Zn,na.rm=TRUE),
            mean_Mn = mean(Mn,na.rm=TRUE),
            mean_Cu = mean(Cu,na.rm=TRUE),
            mean_Fe = mean(Fe,na.rm=TRUE),
            mean_Ca = mean(Ca,na.rm=TRUE),
            mean_Mg = mean(Mg,na.rm=TRUE),
            mean_Na = mean(Na,na.rm=TRUE),
            mean_TKN = mean(TKN,na.rm=TRUE),
            mean_C_N_ratio = mean(`C:N`,na.rm=TRUE),
            mean_Sand = mean(Sand,na.rm=TRUE),
            mean_Silt = mean(Silt,na.rm=TRUE),
            mean_Clay = mean(Clay,na.rm=TRUE))


# join soil chemistry data with sample metadata
meta <- microbiome::meta(ps) %>% full_join(chem,by="inoculum") %>% sample_data()
sample_names(meta) <- meta$sample_name
ps@sam_data <- sample_data(meta)

# fix "burn frequency" variable to ordered factor
ps@sam_data$inoculum_burn_freq <- 
  ps@sam_data$inoculum_burn_freq %>% as.numeric() %>% ordered()
ps_full@sam_data$inoculum_burn_freq <- 
  ps_full@sam_data$inoculum_burn_freq %>% as.numeric() %>% ordered()

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./Output/ITS_ASV_reference_sequences.RDS")

seqs <- taxa_names(inoc)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./Output/ITS_Inoculum_ASV_reference_sequences.RDS")

saveRDS(ps, file = "./Output/ITS_clean_phyloseq_object.RDS")
saveRDS(inoc, file = "./Output/ITS_inoculum_samples_clean_phyloseq_object.RDS")

