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
#                     zahntools v 0.1.0
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
library(zahntools);packageVersion('zahntools')
'%ni%' <- Negate('%in%')

# seed
set.seed(666)

# DATA ####
ps <- readRDS("./Output/18S_ps_not-cleaned.RDS") # change to non-phylogeny stuff

# clean up wilting_scale
ps@sam_data$wilting_scale <- (1 / ps@sam_data$wilting_scale)


# Taxonomic ranks from PR2/Maarjam
colnames(tax_table(ps)) <- c("Domain","Supergroup","Division","Kingdom","Phylum","Subdivision","Order","Genus","Species")

# Remove Non-Fungi?
ps <- ps %>% subset_taxa(Kingdom == "Fungi")
ps@tax_table[,5] %>% unique

##### STOPPED HERE Feb 9, 4:20pm - Need to go back and export positive controls #####

# check positive control efficacy...
pos_ctls <- c("AMF-P-3","AMF-P-4","AMF-P-5","AMF-Pos1","AMF-Pos2","AMF-Pos3","AMF-Pos4","AMF-Pos5")
pos <- ps %>% 
  subset_samples(sample_names(ps) %in% pos_ctls)
y <- tax_table(pos)[which(taxa_sums(pos) > 0),5:9] %>% as.data.frame()


x <- as.data.frame(otu_table(pos))[,which(taxa_sums(pos) > 0)]
y$abund <- as.numeric(x[1,])

y$Species %>% 
  table

y %>% group_by(Species) %>% summarize(total=sum(abund)) %>% arrange(desc(total)) %>% mutate(relabund=total/sum(total)) %>% 
  mutate(Species=factor(Species,levels=pluck(.,"Species"))) %>% 
  ggplot(aes(x=Species,y=relabund)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=60,hjust=1)) +
  labs(y="Relative abundance",title="Fungal positive control sample")
ggsave("./Output/figs/18S_Fungal_positive_control_taxonomy.png",height = 4,width = 6)

pos %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill = "Genus") +
  labs(y="Relative abundance")
ggsave("./Output/figs/18S_Positive_Control_stacked-barchart.png", height = 6, width = 8,dpi=300)


# Remove postive control
ps <- ps %>% subset_samples(sample_names(ps) %ni% pos_ctls)
ps@sam_data$type <- ifelse(is.na(ps@sam_data$other_frompreviouscolumn),"Root","Inoculum")

# SEPARATE INOCULUM SAMPLES ####
# back-up full ps object
ps_full <- ps

# pull out inoculum reads
inoc <- subset_samples(ps, type == "Inoculum")
# remove empty taxa/samples
inoc <- subset_taxa(inoc, taxa_sums(inoc) > 0)
inoc <- subset_samples(inoc, sample_sums(inoc) > 0)

# Remove inoculum samples
ps <- subset_samples(ps, type == "Root")
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
saveRDS(site1.inoc.full,"Output/phyloseq_objects/18S_site1.inoc.full.RDS")
saveRDS(site2.inoc.full,"Output/phyloseq_objects/18S_site2.inoc.full.RDS")
saveRDS(site3.inoc.full,"Output/phyloseq_objects/18S_site3.inoc.full.RDS")
saveRDS(site4.inoc.full,"Output/phyloseq_objects/18S_site4.inoc.full.RDS")
saveRDS(site5.inoc.full,"Output/phyloseq_objects/18S_site5.inoc.full.RDS")
saveRDS(site6.inoc.full,"Output/phyloseq_objects/18S_site6.inoc.full.RDS")
saveRDS(site1.inoc,"Output/phyloseq_objects/18S_site1.inoc.RDS")
saveRDS(site2.inoc,"Output/phyloseq_objects/18S_site2.inoc.RDS")
saveRDS(site3.inoc,"Output/phyloseq_objects/18S_site3.inoc.RDS")
saveRDS(site4.inoc,"Output/phyloseq_objects/18S_site4.inoc.RDS")
saveRDS(site5.inoc,"Output/phyloseq_objects/18S_site5.inoc.RDS")
saveRDS(site6.inoc,"Output/phyloseq_objects/18S_site6.inoc.RDS")
# full objects
saveRDS(ps, file = "./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")
saveRDS(inoc, file = "./Output/phyloseq_objects/18S_inoculum_samples_clean_phyloseq_object.RDS")

