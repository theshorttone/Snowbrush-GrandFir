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
#                     vegan v 2.6.4
#                     patchwork v 1.1.3
#                     ecodist v 2.1.3
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
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(ecodist); packageVersion("ecodist")

# options
theme_set(theme_minimal())
source("./R/palettes.R")

# functions

find_remaining_taxa <- function(z){
  
  # inoculum
  i <- 
    z %>% 
    subset_samples(community == "inoculum") %>% 
    subset_taxa(taxa_sums(z %>% 
                            subset_samples(community == "inoculum")) > 0) %>% 
    transform_sample_counts(function(x){x/sum(x)})
  i@phy_tree <- NULL
  
  # final
  f <- 
    z %>% 
    subset_samples(community == "final") %>% 
    subset_taxa(taxa_sums(z %>% 
                            subset_samples(community == "final")) > 0) #%>%
  # transform_sample_counts(function(x){x/sum(x)})
  f@phy_tree <- NULL
  
  # remaining taxa from inoculum
  r <- taxa_names(f)[which(taxa_names(f) %in% taxa_names(i))]
  
  fr <- 
    f %>% 
    subset_taxa(taxa_names(f) %in% r) #%>% 
  # transform_sample_counts(function(x){x/sum(x)})
  # ir <- 
  #   i %>% 
  #   subset_taxa(taxa_names(i) %in% r) %>% 
  #   transform_sample_counts(function(x){x/sum(x)})
  
  # merge initial and final subsets
  full <- merge_phyloseq(fr,i)
  return(full)
}

plot_remaining_taxa_bar <- function(x){
  x %>% 
    plot_bar(fill="Phylum") +
    facet_wrap(~community*drought,scales = 'free') +
    labs(y="Relative abundance") +
    theme(strip.text = element_text(face='bold',size=14),
          axis.title = element_text(face='bold',size=12))
}

plot_remaining_heatmap <- function(z){
  x <- z %>% 
    merge_samples("community") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    otu_table() %>% 
    as("matrix") 
  x[is.nan(x)] <- 0
  x <- x[,colSums(x) != 0]
  x <- x %>% as.data.frame()
  
  max_95 <- x %>% 
    mutate(group=row.names(.)) %>% 
    pivot_longer(-group) %>% 
    pluck("value") %>% 
    quantile(.95)
  
  x %>% 
    mutate(group=row.names(.)) %>% 
    pivot_longer(-group) %>% 
    mutate(newvalue = ifelse(value > max_95,max_95,value)) %>% 
    ggplot(aes(x=group,y=name,fill=newvalue)) +
    geom_tile() +
    theme(axis.text.y = element_blank(),axis.text.x = element_text(face='bold',size=12)) +
    scale_fill_viridis_c() +
    labs(fill="Relative\nabundance",x="",y="ASV")
  
}

MRM_remaining_taxa <- function(z){
  fin <- 
    z %>% 
    subset_samples(community == "final") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    otu_table() %>% 
    as("Matrix") %>% 
    t() %>% 
    dist()
  ini <- 
    z %>% 
    subset_samples(community == "inoculum") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    otu_table() %>% 
    as("Matrix") %>% 
    t() %>% 
    dist()
  
  return(ecodist::MRM(fin ~ ini))
}
# seed
set.seed(666)

# Load Data
site1.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site1.inoc.full.RDS")
site2.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site2.inoc.full.RDS")
site3.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site3.inoc.full.RDS")
site4.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site4.inoc.full.RDS")
site5.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site5.inoc.full.RDS")
site6.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site6.inoc.full.RDS")

# ORDINATIONS (DCA) ####
ord_1 <- 
  site1.inoc.full %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DCA",distance = "unifrac")
inoc_plot_1 <- plot_ordination(site1.inoc.full,ord_1,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_1,"./Output/figs/ITS_inoc_plot_1.RDS")

ord_2 <- 
  site2.inoc.full %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DCA",distance = "unifrac")
inoc_plot_2 <- plot_ordination(site2.inoc.full,ord_2,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_2,"./Output/figs/ITS_inoc_plot_2.RDS")

ord_3 <- 
  site3.inoc.full %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DCA",distance = "unifrac")
inoc_plot_3 <- plot_ordination(site3.inoc.full,ord_3,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_3,"./Output/figs/ITS_inoc_plot_3.RDS")

ord_4 <- 
  site4.inoc.full %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DCA",distance = "unifrac")
inoc_plot_4 <- plot_ordination(site4.inoc.full,ord_4,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_4,"./Output/figs/ITS_inoc_plot_4.RDS")

ord_5 <- 
  site5.inoc.full %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DCA",distance = "unifrac")
inoc_plot_5 <- plot_ordination(site5.inoc.full,ord_5,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_5,"./Output/figs/ITS_inoc_plot_5.RDS")

ord_6 <- 
  site6.inoc.full %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "DCA",distance = "unifrac")
inoc_plot_6 <- plot_ordination(site6.inoc.full,ord_6,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_6,"./Output/figs/ITS_inoc_plot_6.RDS")

# FIND OVERLAP BETWEEN INOCULUM AND FINAL ####
z <- site1.inoc.full
successful_1 <- find_remaining_taxa(site1.inoc.full)
z <- site2.inoc.full
successful_2 <- find_remaining_taxa(site2.inoc.full)
z <- site3.inoc.full
successful_3 <- find_remaining_taxa(site3.inoc.full)
z <- site4.inoc.full
successful_4 <- find_remaining_taxa(site4.inoc.full)
z <- site5.inoc.full
successful_5 <- find_remaining_taxa(site5.inoc.full)
z <- site6.inoc.full
successful_6 <- find_remaining_taxa(site6.inoc.full)

# PLOT HEATMAPS ####
# color scale maxed out at 95th percentile
p1 <- plot_remaining_heatmap(successful_1) + ggtitle("Inoculum 1") + theme(legend.position = 'none')
p2 <- plot_remaining_heatmap(successful_2) + ggtitle("Inoculum 2") + theme(legend.position = 'none')
p3 <- plot_remaining_heatmap(successful_3) + ggtitle("Inoculum 3") + theme(legend.position = 'none')
p4 <- plot_remaining_heatmap(successful_4) + ggtitle("Inoculum 4") + theme(legend.position = 'none')
p5 <- plot_remaining_heatmap(successful_5) + ggtitle("Inoculum 5") + theme(legend.position = 'none')
p6 <- plot_remaining_heatmap(successful_6) + ggtitle("Inoculum 6") + theme(legend.position = 'none')

(p1 + p2 + p3) / (p4 + p5 + p6) +plot_layout(guides = "collect") 
ggsave("./Output/figs/ITS_Inoculum_Taxa_Before_and_After.png",width = 8,height = 6)


# MRM ####
# does initial community predict final community (only taxa found initially in inoculum)?
mrm_1 <- MRM_remaining_taxa(successful_1)
mrm_2 <- MRM_remaining_taxa(successful_2)
mrm_3 <- MRM_remaining_taxa(successful_3)
mrm_4 <- MRM_remaining_taxa(successful_4)
mrm_5 <- MRM_remaining_taxa(successful_5)
mrm_6 <- MRM_remaining_taxa(successful_6)

MRM_df <- 
  cbind(inoculum=c("Inoc_1","Inoc_2","Inoc_3","Inoc_4","Inoc_5","Inoc_6"),
        rbind(
          unlist(mrm_1),
          unlist(mrm_2),
          unlist(mrm_3),
          unlist(mrm_4),
          unlist(mrm_5),
          unlist(mrm_6)) 
  ) %>% 
  as.data.frame() %>% 
  mutate(coef1=as.numeric(coef1),
         coef2=as.numeric(coef2),
         coef3=as.numeric(coef3),
         coef4=as.numeric(coef4),
         r.squared.R2=as.numeric(r.squared.R2),
         r.squared.pval=as.numeric(r.squared.pval),
         F.test.F=as.numeric(F.test.F),
         F.test.F.pval=as.numeric(F.test.F.pval)) %>% 
  mutate(across(where(is.numeric),function(x){round(x,3)}))
saveRDS(MRM_df,"./Output/ITS_MRM_stats_table.RDS")



# BARPLOTS OF REMAINING TAXA ####
plot_remaining_taxa_bar(successful_1)
plot_remaining_taxa_bar(successful_2)
plot_remaining_taxa_bar(successful_3)
plot_remaining_taxa_bar(successful_4)
plot_remaining_taxa_bar(successful_5)
plot_remaining_taxa_bar(successful_6)


