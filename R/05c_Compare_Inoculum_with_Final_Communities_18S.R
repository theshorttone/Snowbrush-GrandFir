# -----------------------------------------------------------------------------#
# Determining ASV transplantation success (ITS)
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.46.0
#                     ShortRead v 1.60.0
#                     Biostrings v 2.70.1
#                     janitor v 2.2.0
#                     microbiome v 1.24.0
#                     vegan v 2.6.4
#                     patchwork v 1.1.3
#                     ecodist v 2.1.3
#                     corncob v 0.4.1
#                     patchwork v 1.1.3
# -----------------------------------------------------------------------------#

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
library(zahntools); packageVersion('zahntools')#github: gzahn/zahntools
# seed
set.seed(666)

# options and variables
theme_set(theme_minimal())
source("./R/palettes.R")
'%ni%' <- Negate('%in%')

# FUNCTIONS ####

condense_ps_to_species <- function(x){
  # remove bad taxa assignments
  x <- x %>% subset_taxa(!is.na(Phylum))
  # build data frame
  sp <- x@tax_table[,7] %>% as.character()
  gn <- x@tax_table[,6] %>% as.character()
  or <- x@tax_table[,4] %>% as.character()
  sb <- x@tax_table[,3] %>% as.character()
  ph <- x@tax_table[,2] %>% as.character()
  
  condensed_taxonomy <- 
    data.frame(ph,sb,or,gn,sp) %>% 
    mutate(spp = case_when(is.na(or) & is.na(gn) & is.na(sp) ~ paste0(sb," sp."),
                           !is.na(or) & is.na(gn) & is.na(sp) ~ paste0(or," sp."),
                           !is.na(or) & !is.na(gn) & is.na(sp) ~ paste0(gn," sp."),
                           !is.na(or) & !is.na(gn) & !is.na(sp) ~ sp))
  
  x@tax_table[,7] <- str_replace_all(condensed_taxonomy$spp,"_"," ")
  x <- tax_glom(x,"Species")
  return(x)
}


find_remaining_taxa <- function(z){
  
  # inoculum
  i <- 
    z %>% 
    subset_samples(community == "inoculum") %>% 
    subset_taxa(taxa_sums(z %>% 
                            subset_samples(community == "inoculum")) > 0) %>% 
    # should we merge the inoculum samples to give an easier overview?
    merge_samples("community",fun = 'sum') %>% 
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
    subset_taxa(taxa_names(f) %in% r) 
  
  fr <- 
    fr %>% 
    # subset_taxa(taxa_sums(fr) > 0) %>% 
    transform_sample_counts(function(x){x/sum(x)})
  
  # return list of ps objects
  output <- 
    list(inoculum_ra = i,
    final_ra = fr)
  
  # merge initial and final subsets
  # full <- merge_phyloseq(fr,i)
  return(output)
}


build_remaining_dfs <- function(z, inoc.site, only.final.taxa = FALSE){
  
  inoc_site <- paste0("Inoculum ",inoc.site)
  
  inoc <- 
    z[["inoculum_ra"]] %>% 
    otu_table() %>% 
    as("matrix") %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(taxon_name = corncob::otu_to_taxonomy(taxa_names(z[["inoculum_ra"]]), z[["inoculum_ra"]])) %>% 
    mutate(rel_abund = inoculum,
           sample = "inoculum",
           rel_abund_transformed = ifelse(rel_abund > quantile(rel_abund,.95,na.rm=TRUE),quantile(rel_abund,.95,na.rm=TRUE),rel_abund),
           inoc_source = inoc_site) %>% 
    select(taxon_name,sample,rel_abund,rel_abund_transformed,inoc_source)
  
  drought <- z[["final_ra"]]@sam_data$drought
  host <- z[["final_ra"]]@sam_data$species
  sample <- sample_names(z[["final_ra"]])
  
  meta <- data.frame(drought,host,sample) %>% 
    mutate(host_sciname = case_when(host == "GrandFir" ~ "A. grandis",
                                    host == "Snowbrush" ~ "C. velutinus"))
  
  
  taxa_names(successful_1[[1]])
  
  final <- 
    z[["final_ra"]] %>% 
    otu_table() %>% 
    as("matrix") %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(taxon_name = corncob::otu_to_taxonomy(row.names(.),z[["final_ra"]])) %>% 
    pivot_longer(-taxon_name,names_to = "sample",values_to = "rel_abund") %>% 
    mutate(rel_abund_transformed = ifelse(rel_abund > quantile(rel_abund,.95,na.rm=TRUE),quantile(rel_abund,.95,na.rm=TRUE),rel_abund),
           inoc_source = inoc_site)
  
  final <- full_join(final,meta)
  
  if(only.final.taxa){
    inoc <- 
      inoc %>% 
      dplyr::filter(taxon_name %in% unique(final$taxon_name))
  }
  
  return(list(inoculum = inoc,
              pot = final))
  
}

plot_remaining_taxa_bar <- function(x){
  x %>% 
    plot_bar(fill="Phylum") +
    facet_wrap(~community*drought,scales = 'free') +
    labs(y="Relative abundance") +
    theme(strip.text = element_text(face='bold',size=14),
          axis.title = element_text(face='bold',size=12))
}

plot_remaining_heatmap <- function(z, viridis.option='mako'){
  inoc_plot <- 
    z$inoculum %>% 
    ggplot(aes(x=sample,y=taxon_name,fill=rel_abund_transformed)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_viridis_c(option = viridis.option) +
    labs(x="",
         y="ASVs")  +
    facet_wrap(~ inoc_source) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(face='bold',size=14),
          legend.position = 'none',
          strip.text = element_text(face='bold',size=10),
          axis.ticks = element_blank())
  
  
  final_plot <- 
    z$pot %>% 
    ggplot(aes(x=sample,y=taxon_name,fill=rel_abund_transformed)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_viridis_c(option = viridis.option) +
    labs(x="",y="") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle=90,hjust=1,vjust=.5),
          legend.position = 'none',
          axis.title = element_text(face='bold',size=12)) +
    facet_wrap(~host_sciname*drought,scales = 'free_x',nrow = 1) +
    theme(strip.text = element_text(face='bold.italic',size=10),
          axis.ticks = element_blank())
  
  return(list(inoc = inoc_plot,
              final = final_plot))  
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

find_dropped_taxa <- function(ps_list){
  dropped_seqs <- taxa_names(ps_list$inoculum_ra)[taxa_names(ps_list$inoculum_ra) %ni% taxa_names(ps_list$final_ra)]
  if(length(dropped_seqs) == 0){
    dropped_taxa <- dropped_seqs
  } else {
    dropped_taxa <- corncob::otu_to_taxonomy(dropped_seqs,ps_list$inoculum_ra)
  }
    return(dropped_taxa)
}

run_MRM_inoc_final <- 
  function(ps){
    inoc <- 
      ps$inoculum_ra %>% 
      subset_taxa(taxa_names(ps$inoculum_ra) %in% taxa_names(ps$final_ra))
    inoc@sam_data$community <- "inoculum"
    
    final <- ps$final_ra %>% 
      subset_taxa(taxa_names(ps$final_ra) %in% taxa_names(ps$inoculum_ra))
    
    # build distance matrices for both communities
    sample_sums(inoc)
    final_inoc_dist <- final@otu_table %>% as('matrix') %>% t() %>% vegan::vegdist(na.rm = TRUE)
    initial_inoc_dist <- inoc@otu_table %>% as('matrix') %>% t() %>% vegan::vegdist(na.rm = TRUE)
    
    # save ecodist tables
    MRM_results <- 
      ecodist::MRM(final_inoc_dist ~ initial_inoc_dist) %>% 
      as.data.frame() 
    
    return(MRM_results)
  }

plot_sterile_heatmap <- function(comparison){
  sterile %>% 
    tax_glom("Species") %>% 
    transform_sample_counts(function(x){x/sum(x)}) %>% 
    psmelt() %>% 
    mutate(taxon_name = otu_to_taxonomy(OTU,sterile),
           sample = Sample,
           rel_abund = Abundance,
           rel_abund_transformed = ifelse(rel_abund > quantile(rel_abund,.95,na.rm=TRUE),quantile(rel_abund,.95,na.rm=TRUE),rel_abund),
           inoc_source = inoculum) %>% 
    select(all_of(remaining_1_p$inoculum %>% names)) %>% 
    dplyr::filter(taxon_name %in% comparison[["inoculum"]][["taxon_name"]]) %>% 
    ggplot(aes(x=sample,y=taxon_name,fill=rel_abund_transformed)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'mako') +
    labs(x="",
         y="ASVs",
         fill="Relative abundance")  +
    facet_wrap(~ inoc_source) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=90,hjust=1,vjust=.5),
          axis.title.x = element_text(face='bold',size=14),
          legend.position = 'bottom',
          strip.text = element_text(face='bold',size=10),
          axis.ticks = element_blank())
}


# LOAD DATS ####
## Load Data ####
site1.inoc.full <- readRDS("Output/phyloseq_objects/18S_site1.inoc.full.RDS") %>% condense_ps_to_species()
site2.inoc.full <- readRDS("Output/phyloseq_objects/18S_site2.inoc.full.RDS") %>% condense_ps_to_species()
site3.inoc.full <- readRDS("Output/phyloseq_objects/18S_site3.inoc.full.RDS") %>% condense_ps_to_species()
site4.inoc.full <- readRDS("Output/phyloseq_objects/18S_site4.inoc.full.RDS") %>% condense_ps_to_species()
site5.inoc.full <- readRDS("Output/phyloseq_objects/18S_site5.inoc.full.RDS") %>% condense_ps_to_species()
site6.inoc.full <- readRDS("Output/phyloseq_objects/18S_site6.inoc.full.RDS") %>% condense_ps_to_species()
ps <- readRDS("Output/phyloseq_objects/18S_clean_phyloseq_object.RDS") %>% condense_ps_to_species()



# ORDINATIONS (DCA) ####
ord_1 <- 
  site1.inoc.full %>% 
  subset_samples(sample_sums(site1.inoc.full) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  ordinate(method = "NMDS",distance = "bray")
inoc_plot_1 <- plot_ordination(site1.inoc.full,ord_1,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_1,"./Output/figs/18S_inoc_plot_1.RDS")

ord_2 <- 
  site2.inoc.full %>% 
  subset_samples(sample_sums(site2.inoc.full) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",distance = "bray")
inoc_plot_2 <- plot_ordination(site2.inoc.full,ord_2,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_2,"./Output/figs/18S_inoc_plot_2.RDS")

ord_3 <- 
  site3.inoc.full %>% 
  subset_samples(sample_sums(site3.inoc.full) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",distance = "bray")
inoc_plot_3 <- plot_ordination(site3.inoc.full,ord_3,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_3,"./Output/figs/18S_inoc_plot_3.RDS")

ord_4 <- 
  site4.inoc.full %>% 
  subset_samples(sample_sums(site4.inoc.full) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",distance = "bray")
inoc_plot_4 <- plot_ordination(site4.inoc.full,ord_4,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_4,"./Output/figs/18S_inoc_plot_4.RDS")

ord_5 <- 
  site5.inoc.full %>% 
  subset_samples(sample_sums(site5.inoc.full) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",distance = "bray")
inoc_plot_5 <- plot_ordination(site5.inoc.full,ord_5,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_5,"./Output/figs/18S_inoc_plot_5.RDS")

ord_6 <- 
  site6.inoc.full %>% 
  subset_samples(sample_sums(site6.inoc.full) > 0) %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",distance = "bray")
inoc_plot_6 <- plot_ordination(site6.inoc.full,ord_6,color = "community") + scale_color_manual(values=pal.discrete[c(5,2)])
saveRDS(inoc_plot_6,"./Output/figs/18S_inoc_plot_6.RDS")


# FIND OVERLAP BETWEEN INOCULUM AND FINAL ####
# function is not pulling argument from function call for some reason. Whatever, I'm tired and this works:
z=site1.inoc.full
i <- 
  z %>% 
  subset_samples(community == "inoculum") %>% 
  subset_taxa(taxa_sums(z %>% 
                          subset_samples(community == "inoculum")) > 0) %>% 
  merge_samples("community",fun = 'sum') %>% 
  transform_sample_counts(function(x){x/sum(x)})
i@phy_tree <- NULL
f <- 
  z %>% 
  subset_samples(community == "final") %>% 
  subset_taxa(taxa_sums(z %>% 
                          subset_samples(community == "final")) > 0)
f@phy_tree <- NULL
r <- taxa_names(f)[which(taxa_names(f) %in% taxa_names(i))]
fr <- 
  f %>% 
  subset_taxa(taxa_names(f) %in% r) 
fr <- 
  fr %>% 
  transform_sample_counts(function(x){x/sum(x)})
output <- 
  list(inoculum_ra = i,
       final_ra = fr)

z=site1.inoc.full; successful_1 <- find_remaining_taxa(z)
z=site2.inoc.full; successful_2 <- find_remaining_taxa(z)
z=site3.inoc.full; successful_3 <- find_remaining_taxa(z)
z=site4.inoc.full; successful_4 <- find_remaining_taxa(z)
z=site5.inoc.full; successful_5 <- find_remaining_taxa(z)
z=site6.inoc.full; successful_6 <- find_remaining_taxa(z)

## build data frames for each inoculum set ####
remaining_1 <- build_remaining_dfs(successful_1,inoc.site = "1")
remaining_2 <- build_remaining_dfs(successful_2,inoc.site = "2")
remaining_3 <- build_remaining_dfs(successful_3,inoc.site = "3")
remaining_4 <- build_remaining_dfs(successful_4,inoc.site = "4")
remaining_5 <- build_remaining_dfs(successful_5,inoc.site = "5")
remaining_6 <- build_remaining_dfs(successful_6,inoc.site = "6")

successful_2$final_ra@tax_table

remaining_1_p <- build_remaining_dfs(successful_1,inoc.site = "1",only.final.taxa = TRUE)
remaining_2_p <- build_remaining_dfs(successful_2,inoc.site = "2",only.final.taxa = TRUE)
remaining_3_p <- build_remaining_dfs(successful_3,inoc.site = "3",only.final.taxa = TRUE)
remaining_4_p <- build_remaining_dfs(successful_4,inoc.site = "4",only.final.taxa = TRUE)
remaining_5_p <- build_remaining_dfs(successful_5,inoc.site = "5",only.final.taxa = TRUE)
remaining_6_p <- build_remaining_dfs(successful_6,inoc.site = "6",only.final.taxa = TRUE)

# COMPARE W/ STERILE INOC POTS ####
# remove taxa also found in "sterile inoculum" pots
sterile <- ps %>% 
  subset_samples(inoculum_site == "Sterile")
sterile <- sterile %>% 
  subset_taxa(taxa_sums(sterile) > 0)




remaining_1_p$inoculum$taxon_name %in% corncob::otu_to_taxonomy( taxa_names(sterile),sterile)
remaining_2_p$inoculum$taxon_name %in% corncob::otu_to_taxonomy( taxa_names(sterile),sterile)
remaining_3_p$inoculum$taxon_name %in% corncob::otu_to_taxonomy( taxa_names(sterile),sterile)
remaining_4_p$inoculum$taxon_name %in% corncob::otu_to_taxonomy( taxa_names(sterile),sterile)
remaining_5_p$inoculum$taxon_name %in% corncob::otu_to_taxonomy( taxa_names(sterile),sterile)
remaining_6_p$inoculum$taxon_name %in% corncob::otu_to_taxonomy( taxa_names(sterile),sterile)


# PLOT HEATMAPS ####

# use versions where only taxa shown for inoculum are also present in the final pot community
heatmaps_1 <- plot_remaining_heatmap(remaining_1_p)
heatmaps_2 <- plot_remaining_heatmap(remaining_2_p)
heatmaps_3 <- plot_remaining_heatmap(remaining_3_p)
heatmaps_4 <- plot_remaining_heatmap(remaining_4_p)
heatmaps_5 <- plot_remaining_heatmap(remaining_5_p)
heatmaps_6 <- plot_remaining_heatmap(remaining_6_p)

# include heatmaps of taxa also found in pots with sterile inoculum
p1 <- heatmaps_1$inoc + plot_sterile_heatmap(comparison = remaining_1_p) + heatmaps_1$final + 
  plot_layout(widths = c((3/4), 3, 3),guides = 'keep') & theme(legend.position = 'bottom')
p2 <- heatmaps_2$inoc + plot_sterile_heatmap(comparison = remaining_2_p) + heatmaps_2$final  + 
  plot_layout(widths = c((3/4), 3, 3),guides = 'keep') & theme(legend.position = 'bottom')
p3 <- heatmaps_3$inoc + plot_sterile_heatmap(comparison = remaining_3_p) + heatmaps_3$final  + 
  plot_layout(widths = c((3/4), 3, 3),guides = 'keep') & theme(legend.position = 'bottom')
p4 <- heatmaps_4$inoc + plot_sterile_heatmap(comparison = remaining_4_p) + heatmaps_4$final  + 
  plot_layout(widths = c((3/4), 3, 3),guides = 'keep') & theme(legend.position = 'bottom')
p5 <- heatmaps_5$inoc + plot_sterile_heatmap(comparison = remaining_5_p) + heatmaps_5$final  + 
  plot_layout(widths = c((3/4), 3, 3),guides = 'keep') & theme(legend.position = 'bottom')
p6 <- heatmaps_6$inoc + plot_sterile_heatmap(comparison = remaining_6_p) + heatmaps_6$final  + 
  plot_layout(widths = c((3/4), 3, 3),guides = 'keep') & theme(legend.position = 'bottom')

# save heatmaps
# p1; ggsave("./Output/figs/18S_inoculum_taxa_heatmaps_inoc-1.png",dpi=300,height = 8,width = 12)
# p2; ggsave("./Output/figs/18S_inoculum_taxa_heatmaps_inoc-2.png",dpi=300,height = 8,width = 12)
# p3; ggsave("./Output/figs/18S_inoculum_taxa_heatmaps_inoc-3.png",dpi=300,height = 8,width = 12)
# p4; ggsave("./Output/figs/18S_inoculum_taxa_heatmaps_inoc-4.png",dpi=300,height = 8,width = 12)
# p5; ggsave("./Output/figs/18S_inoculum_taxa_heatmaps_inoc-5.png",dpi=300,height = 8,width = 12)
# p6; ggsave("./Output/figs/18S_inoculum_taxa_heatmaps_inoc-6.png",dpi=300,height = 8,width = 12)



# FIND TAXA THAT DIDN'T SURVIVE ####
# Taxa that are in the inoculum, but not in the final pot community

dropped_1 <- find_dropped_taxa(successful_1)
dropped_2 <- find_dropped_taxa(successful_2)
dropped_3 <- find_dropped_taxa(successful_3)
dropped_4 <- find_dropped_taxa(successful_4)
dropped_5 <- find_dropped_taxa(successful_5)
dropped_6 <- find_dropped_taxa(successful_6)
# list of fungi found in all inoculum sources that did not persist in any pots
unsuccessful_transplant_taxa <- 
  dropped_1[dropped_1 %in% c(dropped_2,dropped_3,dropped_4,dropped_5,dropped_6)]
# list of all unique taxa that didn't transplant successfully
unsuccessful_asvs <- 
  list(dropped_1,dropped_2,dropped_3,dropped_4,dropped_5,dropped_6) %>% 
  map(names) %>% 
  map(unique) %>% 
  unlist() %>% 
  unique()

## Taxonomic breakdown of unsuccessful vs successful taxa ####

# load up full fungal phyloseq (inoculum only)
fung <- merge_phyloseq(site1.inoc.full,site2.inoc.full,site3.inoc.full,site4.inoc.full,site5.inoc.full,site6.inoc.full) %>% 
  subset_samples(community == "inoculum")
fung@sam_data$inoculum_site <- 
  fung@sam_data$sample_name %>% str_split("W") %>% map_chr(2)


fung_successful <- fung %>% 
  subset_taxa(taxa_names(fung) %ni% unsuccessful_asvs)
sample_names(fung_successful) <- paste0(sample_names(fung_successful),"_S")


fung_unsuccessful <- fung %>% 
  subset_taxa(taxa_names(fung) %in% unsuccessful_asvs)
sample_names(fung_unsuccessful) <- paste0(sample_names(fung_unsuccessful),"_U")

fung <- merge_phyloseq(fung_successful,fung_unsuccessful)
fung@sam_data$success <- ifelse(sample_names(fung) %>% grepl(pattern = "_U$"),
                                FALSE,TRUE)

fung@sam_data$mergevar <- paste0(fung@sam_data$other_frompreviouscolumn,"_", fung@sam_data$success)
fung_merged <- fung %>% 
  merge_samples("mergevar") %>% 
  transform_sample_counts(function(x){x/sum(x)})
# repair metadata
fung_merged@sam_data$success <- sample_names(fung_merged) %>% str_split("_") %>% map_chr(2)
fung_merged@sam_data$site <- sample_names(fung_merged) %>% str_split("_") %>% map_chr(1)
genera <- fung_merged@tax_table[,6] %>% 
  unique %>% 
  unname
genera
rank_names(fung_merged)
ord <- fung_merged@tax_table[,4] %>% 
  unique %>% 
  unname
ord



melted <- 
  fung_merged %>% 
  psmelt() %>% 
  mutate(success=ifelse(success == "FALSE","Unsuccessful","Successful")) %>%
  # clean up order names...
  mutate(Order = Order %>% str_remove("o__"),
         Order = case_when(is.na(Order) ~ "Undetermined",
                           TRUE ~ Order),
         Family = case_when(is.na(Family) | Family == "unclassified" ~ "Undetermined",
                            TRUE ~ Family %>% str_remove("f__")))
melted <- 
melted %>% 
  mutate(Family = case_when(Family == "unclassified" ~ "Undetermined",
                            TRUE ~ Family))
  
melted$Family %>% unique

  ord_order <- 
  melted %>% 
  group_by(Order) %>% 
  summarize(sum_ra = sum(Abundance,na.rm=TRUE)) %>% 
  arrange(desc(sum_ra)) %>% 
  pluck("Order")

  fam_order <- 
    melted %>%
    dplyr::filter(Phylum == "p__Glomeromycota") %>% 
    group_by(Family) %>% 
    summarize(sum_ra = sum(Abundance,na.rm=TRUE)) %>% 
    arrange(desc(sum_ra)) %>% 
    pluck("Family") %>% 
    str_remove("f__")
  
melted <- melted %>% 
  mutate(Order = factor(Order, levels = ord_order),
         Family = factor(Family, levels = fam_order))
  
melted <- melted %>% dplyr::filter(Phylum == "p__Glomeromycota" )

p <- 
  melted %>% 
  ggplot(aes(x=site,y=Abundance,fill=Family)) +
  geom_col(color='white',linewidth=.1) +
  coord_flip() +
  facet_wrap(~success) +
  scale_fill_viridis_d(option = 'turbo') +
  theme_minimal() +
  theme(strip.text = element_text(face='bold',size = 18),
        axis.text.y = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=18),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5),
        legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        # legend.position = 'right',
        legend.box.margin = margin(0, 10, 0, 0)) +
  labs(x="",y="Relative abundance")
p
ggsave("./Output/figs/18S_successful_vs_unsuccessful_taxa_breakdown.png",width = 12,height = 8)
ggsave("./Output/figs/manuscript_versions/18S_successful_vs_unsuccessful_taxa_breakdown.tiff",width = 12,height = 8,dpi=500)
saveRDS(p,"./Output/figs/18S_successful_vs_unsuccessful_taxa_breakdown.RDS")
p <- readRDS("./Output/figs/18S_successful_vs_unsuccessful_taxa_breakdown.RDS")


melted %>% 
  ggplot(aes(x=success,fill=site)) +
  geom_histogram(alpha=.5,stat = 'count') +
  facet_wrap(~Family)

successful_AMF <- 
  melted %>% 
  dplyr::select(Abundance,success,site,Kingdom,Phylum,Class,Order,Family,Genus,Species) %>% 
  unique.data.frame %>% 
  dplyr::filter(success == "Successful" & Abundance > 0) %>% 
  dplyr::select(site,Kingdom,Phylum,Class,Order,Family,Genus,Species) %>% 
  lapply(function(x){str_remove(x,".__") %>% str_remove(".  ")}) %>% 
  as.data.frame()
write_csv(successful_AMF,"./Output/18S_Successful_dataframe.csv")

unsuccessful_AMF <- 
  melted %>% 
  dplyr::select(Abundance,success,site,Kingdom,Phylum,Class,Order,Family,Genus,Species) %>% 
  unique.data.frame %>% 
  dplyr::filter(success == "Unsuccessful" & Abundance > 0) %>% 
  dplyr::select(site,Kingdom,Phylum,Class,Order,Family,Genus,Species) %>% 
  lapply(function(x){str_remove(x,".__") %>% str_remove(".  ")}) %>% 
  as.data.frame
write_csv(unsuccessful_AMF,"./Output/18S_Unsuccessful_dataframe.csv")

p_b <- readRDS("./Output/figs/16S_successful_vs_unsuccessful_taxa_breakdown_withlines.RDS") +
  labs(x="Bacteria\n",y="\n\n") +
  theme(legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"))
p_f <- readRDS("./Output/figs/ITS_successful_vs_unsuccessful_taxa_breakdown.RDS") +
  labs(x="Fungi\n",y="\n\n")
p_a <- p + labs(x="AMF\n",y="Relative abundance")

p_f / p_b / p_a
ggsave("./Output/figs/manuscript_versions/Successful_vs_Unsuccessful_Inoculum_Taxa_Final.png",width = 12,height = 10)
ggsave("./Output/figs/manuscript_versions/Successful_vs_Unsuccessful_Inoculum_Taxa_Final.tiff",width = 12,height = 10,dpi=500)

fung_merged %>% 
  psmelt() %>% 
  mutate(success=ifelse(success == "FALSE","Unsuccessful","Successful")) %>% 
  dplyr::filter(Subdivision == "Glomeromycotina") %>% 
  mutate(Genus = case_when(is.na(Genus) ~ "Uncertain",
                           TRUE ~ Genus)) %>% 
  ggplot(aes(x=site,y=Abundance,fill=Genus)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~success) +
  scale_fill_viridis_d() +
  theme(strip.text = element_text(face='bold',size = 12),
        axis.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=12)) +
  labs(x="",y="Relative abundance",title = "Taxonomy of successfully vs \nunsuccessfully transplanted Glomeromycotina ASVs")
ggsave("./Output/figs/18S_successful_vs_unsuccessful_taxa_breakdown.png",height = 6, width = 10)

fung_merged %>% 
  psmelt() %>% 
  mutate(success=ifelse(success == "FALSE","Unsuccessful","Successful")) %>% 
  # dplyr::filter(Subdivision == "Glomeromycotina") %>% 
  mutate(Genus = case_when(is.na(Genus) ~ "Uncertain",
                           TRUE ~ Genus),
         Sample = str_remove(Sample,"_TRUE")) %>% 
  dplyr::filter(success == "Successful") %>% 
  select(Species,Sample,Abundance,OTU) %>% 
  dplyr::filter(Abundance > 0) %>% 
  group_by(Species,Sample,OTU) %>% 
  summarize(Abundance = mean(Abundance)) %>% 
  arrange(Species,Sample,desc(Abundance)) %>% 
  select(Species,Sample,OTU) %>% 
  saveRDS("./Output/18S_Successfully_Transplanted_Taxa.RDS")




data.frame(
  Taxonomy = corncob::otu_to_taxonomy(unsuccessful_asvs,fung),
  ASV = names(corncob::otu_to_taxonomy(unsuccessful_asvs,fung))
) %>% saveRDS("./Output/18S_unsuccessful_ASV_transplants.RDS")




# MRM ####
# does initial community predict final community (only taxa found initially in inoculum)?

ps <- successful_1; mrm_1 <- run_MRM_inoc_final(ps)
ps <- successful_2; mrm_2 <- run_MRM_inoc_final(ps)
# ps <- successful_3; mrm_3 <- run_MRM_inoc_final(ps)
ps <- successful_4; mrm_4 <- run_MRM_inoc_final(ps)
# ps <- successful_5; mrm_5 <- run_MRM_inoc_final(ps)
# ps <- successful_6; mrm_6 <- run_MRM_inoc_final(ps)

# Pull MRM results tables into single table
MRM_df <- 
  cbind(inoculum=c("Inoc_1","Inoc_2","Inoc_4"),
        rbind(
          unlist(mrm_1 %>% filter(row.names(.) != "Int")),
          unlist(mrm_2 %>% filter(row.names(.) != "Int")),
          # unlist(mrm_3 %>% filter(row.names(.) != "Int")),
          unlist(mrm_4 %>% filter(row.names(.) != "Int"))
          # unlist(mrm_5 %>% filter(row.names(.) != "Int")),
          # unlist(mrm_6 %>% filter(row.names(.) != "Int"))) 
  )) %>% 
  as.data.frame() %>% 
  mutate(coef.final_inoc_dist=as.numeric(coef.final_inoc_dist),
         coef.pval=as.numeric(coef.pval),
         r.squared=as.numeric(r.squared),
         F.test=as.numeric(F.test)) %>% 
  mutate(across(where(is.numeric),function(x){round(x,3)}))
saveRDS(MRM_df,"./Output/18S_MRM_stats_table_inoc_vs_final.RDS")


# CHECK SUCCESSFUL AGAINST STERILE ####
# Check whether "successful ASVs" are also found in sterile control samples

# Get full list of successful ASVs that made it into root samples
successful_taxa <- fung_successful %>% subset_taxa(taxa_sums(fung_successful) > 0) %>% taxa_names

# Get full list of ASVs from Sterile Control root samples
x <- readRDS("Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")
x <- x %>% subset_samples(inoculum_site == "Sterile")
sterile_taxa <- x %>% subset_taxa(taxa_sums(x) > 0) %>% taxa_names

# Find matches
data.frame(successful_asv = successful_taxa,
           taxonomy = corncob::otu_to_taxonomy(successful_taxa,fung),
           found_in_sterile = successful_taxa %in% sterile_taxa) %>% 
  write_csv("./Output/18S_successful_taxa_comparison_w_sterile.csv")

# Remove "successful" transplants that were also found in sterile controls
sterile_contams <- 
  data.frame(successful_asv = successful_taxa,
             taxonomy = corncob::otu_to_taxonomy(successful_taxa,fung),
             found_in_sterile = successful_taxa %in% sterile_taxa) %>% 
  dplyr::filter(found_in_sterile)


# find environmental or other predictors of successful taxa
# do any of these taxa  associate with drought, fire, block? what?
successful_taxa <- 
  data.frame(successful_asv = successful_taxa,
             taxonomy = corncob::otu_to_taxonomy(successful_taxa,fung),
             found_in_sterile = successful_taxa %in% sterile_taxa) %>% 
  dplyr::filter(!found_in_sterile) %>% 
  pluck("successful_asv")
saveRDS(successful_taxa,"./Output/Successful_ASVs_18S.RDS")


fung_merged %>% 
  psmelt() %>% 
  mutate(success=ifelse(success == "FALSE","Unsuccessful","Successful")) %>% 
  dplyr::filter(OTU %ni% sterile_contams$successful_asv) %>% 
  ggplot(aes(x=site,y=Abundance,fill=Order)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~success) +
  scale_fill_viridis_d() +
  theme(strip.text = element_text(face='bold',size = 12),
        axis.text = element_text(face='bold',size=12),
        axis.title = element_text(face='bold',size=12)) +
  labs(x="",y="Relative abundance",title = "Taxonomy of successfully vs \nunsuccessfully transplanted ASVs")
ggsave("./Output/figs/18S_successful_vs_unsuccessful_taxa_breakdown.png",width = 12,height = 8)

# inoc diversity vs final diversity

get_inoc_diversity <- function(x){
  x <- get(x)
  inoc_id <- na.omit(x@sam_data$inoculum_site %>% unique) %>% as.character()
  y <- x %>% 
    subset_samples(!grepl("-[A,B,C]-",x@sam_data$sample_name)) %>% 
    estimate_richness(measures = c("shannon")) %>% 
    mutate(inoculum_site=inoc_id)
    return(y)
}

inoc_shannon <- map(ls(pattern = "inoc.full"),get_inoc_diversity) %>% 
  purrr::reduce(full_join) %>% 
  group_by(inoculum_site) %>% 
  summarize(inoc_med = median(Shannon,na.rm=TRUE),
            inoc_mean = mean(Shannon,na.rm=TRUE),
            inoc_sd = sd(Shannon,na.rm=TRUE))

final <- readRDS("Output/phyloseq_objects/18S_clean_phyloseq_object.RDS")
final_shannon <- estimate_richness(final,measures = "shannon")
final_meta <- microbiome::meta(final)
final_meta$shannon <- final_shannon$Shannon
final_meta %>% 
  dplyr::select(inoculum_site,shannon) %>% 
  full_join(inoc_shannon) %>% unique.data.frame() %>% 
  filter(inoculum_site != "Sterile") %>% 
  mutate(inoc_min = inoc_mean - inoc_sd,inoc_max=inoc_mean + inoc_sd,
         inoc_min = ifelse(inoc_min < 0,0,inoc_min)) %>% 
  lmerTest::lmer(data=.,formula=inoc_mean ~ shannon + (1|inoculum_site)) %>% 
  summary
  # ggplot(aes(x=shannon,y=inoc_mean)) + 
  # geom_point() +
  # facet_wrap(~inoculum_site)
