# -----------------------------------------------------------------------------#
# Validating taxonomy with positive controls
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
#                     zahntools v 0.1.0
# -----------------------------------------------------------------------------#

# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(adegenet); packageVersion("adegenet")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")
library(microbiome); packageVersion("microbiome")
library(zahntools); packageVersion("zahntools")

# Load data and subset to pos ctls
bact <- readRDS("./Output/16S_ps_not-cleaned_w_tree.RDS") 
bact <- bact %>% subset_samples(sample_names(bact) == "P1")
bact <- bact %>% subset_taxa(taxa_sums(bact) > 0)

its <- readRDS("./Output/ITS_ps_not-cleaned.RDS")
its <- its %>% subset_samples(sample_names(its) == "F-P-2")
its <- its %>% subset_taxa(taxa_sums(its) > 0)

amf <- readRDS("./Output/18S_ps_not-cleaned.RDS")
amf <- amf %>% subset_taxa(!is.na(amf@tax_table[,9]))
amf <- amf %>% subset_samples(sample_names(amf) %in% c("AMF-P-3","AMF-P-4","AMF-P-5","AMF-Pos1","AMF-Pos2","AMF-Pos3","AMF-Pos4","AMF-Pos5"))
amf <- amf %>% subset_taxa(taxa_sums(amf) > 0)

# clean up genus/species names  
bact@tax_table[,7] <- paste0(bact@tax_table[,6]," ",bact@tax_table[,7]) %>% str_replace("NA NA","Undetermined") %>% str_replace("NA","sp.")
its@tax_table[,6] <- its@tax_table[,6] %>% str_remove("g__")
its@tax_table[,6][is.na(its@tax_table[,6])] <- "Undetermined"
its@tax_table[,7] <- paste(str_remove(its@tax_table[,6],"g__"),str_remove(its@tax_table[,7],"s__")) %>% str_replace("NA NA","Undescribed") %>% str_replace("NA","sp.")
amf@tax_table[,9] <- amf@tax_table[,9] %>% str_replace_all("_"," ") %>% str_split("\\d+") %>% map_chr(1) %>% str_squish()



# Quick taxa plots
bact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill = "Species")+ scale_fill_viridis_d()
ggsave("./Output/figs/16S_Pos_Control_Barplot.png",width = 3,height = 4)

its %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill = "Genus") + scale_fill_viridis_d()
ggsave("./Output/figs/ITS_Pos_Control_Barplot.png",width = 3,height = 4)

amf %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill = "Species") + scale_fill_viridis_d()
ggsave("./Output/figs/18S_Pos_Control_Barplot.png",width = 8,height = 4)

# Write csv file of all positive control taxa from each amplicon
bact@sam_data
psmelt(bact) %>% 
  select(Sample,Species,OTU) %>% mutate(Amplicon="16S") %>% 
  full_join(psmelt(its) %>% 
              select(Sample,Species=Genus,OTU) %>% mutate(Amplicon="ITS")) %>%  
  full_join(psmelt(amf) %>% 
              select(Sample,Species,OTU) %>% mutate(Amplicon="18S")
  ) %>% 
  write_csv("./Output/Positive_Control_Table.csv")



