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

# Check AMF taxonomy assignments from both databases
amf_rdp <- readRDS("./Output/18S_postive_ctl_RDP_Taxonomy_from_dada2.RDS")
amf_euk <- readRDS("./Output/18S_RDP_positive_ctl_Taxonomy_from_dada2_Eukaryome.RDS")


# Load data and subset to pos ctls
bact <- readRDS("./Output/16S_ps_not-cleaned_w_tree.RDS") 
bact <- bact %>% subset_samples(sample_names(bact) == "P1")
bact <- bact %>% subset_taxa(taxa_sums(bact) > 0)

its <- readRDS("./Output/ITS_ps_not-cleaned.RDS")
its <- its %>% subset_samples(sample_names(its) == "F-P-2")
its <- its %>% subset_taxa(taxa_sums(its) > 0)

marjaam <- readRDS("./Output/18S_pos_ctl_marjaam_physeq.RDS") %>% subset_taxa(Class == "Fungi")
eukaryome <- readRDS("./Output/18S_pos_ctl_eukaryome_physeq.RDS") %>% subset_taxa(Kingdom == "k__Fungi")

eukaryome %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar(fill="Order") +
  scale_fill_viridis_d(option = 'turbo')


marjaam %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar(fill="Genus") +
  scale_fill_viridis_d(option = 'turbo')

p5_marjaam <- marjaam %>% subset_samples(sample_names(marjaam) == "AMF-P-5")
p5_marjaam <- p5_marjaam %>% subset_taxa(taxa_sums(p5_marjaam) > 0)
p5_eukaryome <- eukaryome %>% subset_samples(sample_names(eukaryome) == "AMF-P-5")
p5_eukaryome <- p5_eukaryome %>% subset_taxa(taxa_sums(p5_eukaryome) > 0)

taxa_names(p5_marjaam)


data.frame(
  marjaam_taxonomy = paste(p5_marjaam@tax_table[,5],p5_marjaam@tax_table[,6],p5_marjaam@tax_table[,7],
                           p5_marjaam@tax_table[,8],p5_marjaam@tax_table[,9]),
  eukaryome_taxonomy = paste(p5_eukaryome@tax_table[,2],p5_eukaryome@tax_table[,3],p5_eukaryome@tax_table[,4],
                             p5_eukaryome@tax_table[,5],p5_eukaryome@tax_table[,6]) %>% str_remove_all(".__"),
  top_blast_hit = c("Rhizophagus","Glomus","Pleosporaceae sp.","Pleosporaceae sp.","Rhizophagus or Glomus sp.","Pleosporaceae sp."),
  asv = taxa_names(p5_marjaam)
  
) %>% 
  write_csv("./Output/AMF_P5_taxonomy.csv")

ShortRead()
asvs <- DNAStringSet(taxa_names(p5_marjaam))
names(asvs) <- paste0("ASV_",1:6)


library(DECIPHER)
glomeros <- readFasta("./Taxonomy/glomeromcota.fasta")
names(glomeros) <- glomeros@id
c(glomeros, asvs)
glomeros_reads <- DNAStringSet(glomeros@sread)
names(glomeros_reads) <- glomeros@id

aln <- DECIPHER::AlignSeqs(c(RemoveGaps(glomeros_reads),RemoveGaps(asvs)),)


# tax_table(amf)[14,6] 
amf <- amf %>% subset_taxa(!is.na(amf@tax_table[,6]))
amf <- amf %>% subset_samples(sample_names(amf) %in% c("AMF-P-3","AMF-P-4","AMF-P-5","AMF-Pos1","AMF-Pos2","AMF-Pos3","AMF-Pos4","AMF-Pos5"))
sample_names(amf)

amf <- amf %>% subset_taxa(taxa_sums(amf) > 0)
plot_bar(amf %>% transform_sample_counts(function(x){x/sum(x)}),fill="Phylum")




amf <- amf %>% #subset_samples(sample_names(amf) == "AMF-Pos5") %>% 
  subset_taxa(Phylum == "p__Solanum")
otu_table(amf)


amf = amf %>% 
  subset_taxa(taxa_sums(amf) > 0)

plot_bar(amf,fill="Family")

tax_table(amf)[,7] %>% unique %>% unname

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



