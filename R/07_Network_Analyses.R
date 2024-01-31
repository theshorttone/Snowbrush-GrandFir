# -----------------------------------------------------------------------------#
# Network Analyses for both bacteria and fungi
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
#                     patchwork v 1.1.3
#                     igraph v 1.5.1
#                     SpiecEasi v 1.1.3
#                     zahntools v 0.1.0
#                     parallel v 4.3.2
#                     hubfindr v 0.1.0
# -----------------------------------------------------------------------------#

# SETUP ####
# devtools::install_github("zdk123/SpiecEasi")
# packages 
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(corncob); packageVersion("corncob")
library(patchwork); packageVersion("patchwork")
library(igraph); packageVersion("igraph")
library(SpiecEasi); packageVersion("SpiecEasi")
library(zahntools); packageVersion("zahntools") # devtools::install_github("gzahn/zahntools")
library(parallel); packageVersion("parallel")
library(hubfindr); packageVersion("hubfindr") # devtools::install_github("gzahn/hubfindr")


# Options
options(scipen=999)
options(warn=0)
set.seed(666)
theme_set(theme_bw() + theme(
  strip.text = element_text(face='bold',size=14),
  axis.title = element_text(face='bold',size=14),
  legend.title = element_text(face='bold',size=14),
  axis.text = element_text(face='bold',size=12)
))

source("./R/palettes.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]


# functions
'%ni%' <- Negate('%in%')

find_ig_subset_attr <- function(ps.subset,ig.full){
  
  # run checks and tests for function
  stopifnot(class(ps.subset) == "phyloseq")
  stopifnot(class(ig.full) == "igraph")
  
  # get present taxa in physeq subset
  present_vertices <- which(taxa_sums(ps.subset) > 0)

  # build subgraph from present taxa only
  ps.subset %>% ntaxa()
  ig.subset <- igraph::subgraph(graph = ig.full,vids = present_vertices)
  # calculate various ig params
  # mean_alpha_centrality <- mean(igraph::alpha_centrality(ig.subset,tol=0),na.rm = TRUE)
  clique_num <- igraph::clique_num(ig.subset)
  mean_betweenness <- mean(igraph::betweenness(ig.subset),na.rm = TRUE)
  mean_closeness <- mean(igraph::closeness(ig.subset),na.rm = TRUE)
  mean_coreness <- mean(igraph::coreness(ig.subset),na.rm = TRUE)
  deg_dist <- igraph::degree_distribution(ig.subset)
  global_effic <- igraph::global_efficiency(ig.subset)
  n_vertices <- igraph::vcount(ig.subset)
  n_edges <- igraph::ecount(ig.subset)
  mean_dist <- igraph::mean_distance(ig.subset)
  similarity_matrix <- igraph::similarity(ig.subset,method='jaccard')
  clustering_coeficient <- igraph::transitivity(ig.subset)
  max_degree <- max(igraph::degree(ig.subset),na.rm=TRUE)
  
  out <- list(ig=ig.subset,
              n_vertices=n_vertices,
              n_edges=n_edges,
              mean_dist=mean_dist,
              # mean_alpha_centrality=mean_alpha_centrality,
              clique_num=clique_num,
              mean_betweenness=mean_betweenness,
              mean_closeness=mean_closeness,
              mean_coreness=mean_coreness,
              deg_dist=deg_dist,
              global_effic=global_effic,
              similarity_matrix=similarity_matrix,
              clustering_coeficient=clustering_coeficient,
              max_degree=max_degree)
  
  
  return(out)
}

# Data
bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object_w_guilds.RDS")
bact <- bact %>% subset_taxa(taxa_sums(bact) > 0)

fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object_w_guilds.RDS")
fung <- fung %>% subset_taxa(taxa_sums(fung) > 0)

# join the two for 'full' microbiome (16S + ITS2)
full <- merge_phyloseq(bact,fung)
sample_names(full)

# Barcharts of bact and fung ####

# merge samples by inoc source for plotting
bact_merged <- bact %>% merge_samples("inoculum_site")
fung_merged <- fung %>% merge_samples("inoculum_site")
# repair metadata
bact_merged@sam_data$inoculum_site <- row.names(bact_merged@sam_data)
fung_merged@sam_data$inoculum_site <- row.names(fung_merged@sam_data)

# transform to relabund and plot
phylum_barplot_by_inoc_bact <- 
bact_merged %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill = "Phylum") +
  theme_minimal() +
  scale_fill_viridis_d(end=.8) +
  labs(x="Inoculum source",y="Relative abundance",title = "Bacteria")
phylum_barplot_by_inoc_fung <- 
fung_merged %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill = "Phylum") +
  theme_minimal() +
  scale_fill_viridis_d(end=.8,begin=.2) +
  labs(x="Inoculum source",y="Relative abundance",title = "Fungi")

saveRDS(phylum_barplot_by_inoc_bact,"./Output/figs/phylum_barplot_by_inoc_bact.RDS")
saveRDS(phylum_barplot_by_inoc_fung,"./Output/figs/phylum_barplot_by_inoc_fung.RDS")





# FULL NETWORK ANALYSES ####

# build presence/absence version for co-occurrence networks
full_pa <- transform_sample_counts(full,function(x){ifelse(x>0,1,0)})

# build relative abundance version
full_ra <- transform_sample_counts(full,function(x){x/sum(x)})




## SpiecEasi ####
se.params <- list(rep.num=20, ncores=(parallel::detectCores()-1))

se.mb.fung <- SpiecEasi::spiec.easi(data = fung,
                                       method='mb',
                                       sel.criterion = "bstars",
                                       pulsar.params=se.params)
saveRDS(se.mb.fung,"./Output/ITS_SpiecEasi_out.RDS")
se.mb.fung <- readRDS("./Output/ITS_SpiecEasi_out.RDS")

# run on bacteria grouped by genus
bact_genus <- tax_glom(bact,"Genus")
se.mb.bact <- SpiecEasi::spiec.easi(data = bact_genus,
                                       method='mb',
                                       sel.criterion = "bstars",
                                       pulsar.params=se.params)
saveRDS(se.mb.bact,"./Output/16S_SpiecEasi_genus_out.RDS")
se.mb.bact <- readRDS("./Output/16S_SpiecEasi_genus_out.RDS")

# get best model and build igraph (include vertex names)
fung_igraph <- adj2igraph(getRefit(se.mb.fung),vertex.attr = list(name=taxa_names(fung)))
bact_igraph <- adj2igraph(getRefit(se.mb.bact),vertex.attr = list(name=taxa_names(bact_genus)))

# find putative hub taxa & plot for entire study
fung_hub_taxa <- hubfindr::find_hubs(graph = fung_igraph,physeq = fung)
bact_hub_taxa <- hubfindr::find_hubs(graph = bact_igraph,physeq = bact)

plot_hubs(fung_igraph)
plot_hubs(bact_igraph)


# SUBSET NETWORKS ####
## Build subsets ####

### subset phyloseqs ####
# fungi
fung_sb_1 <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '1')
fung_sb_2 <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '2')
fung_sb_3 <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '3')
fung_sb_4 <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '4')
fung_sb_5 <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '5')
fung_sb_6 <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '6')
fung_sb_s <- fung %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == 'Sterile')

fung_gf_1 <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '1')
fung_gf_2 <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '2')
fung_gf_3 <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '3')
fung_gf_4 <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '4')
fung_gf_5 <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '5')
fung_gf_6 <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '6')
fung_gf_s <- fung %>% 
  subset_samples(species == "GrandFir" & inoculum_site == 'Sterile')

# bacteria
bact_sb_1 <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '1')
bact_sb_2 <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '2')
bact_sb_3 <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '3')
bact_sb_4 <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '4')
bact_sb_5 <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '5')
bact_sb_6 <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == '6')
bact_sb_s <- bact_genus %>% 
  subset_samples(species == "Snowbrush" & inoculum_site == 'Sterile')

bact_gf_1 <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '1')
bact_gf_2 <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '2')
bact_gf_3 <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '3')
bact_gf_4 <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '4')
bact_gf_5 <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '5')
bact_gf_6 <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == '6')
bact_gf_s <- bact_genus %>% 
  subset_samples(species == "GrandFir" & inoculum_site == 'Sterile')

## Calculate attributes ####
# function that can handle pulling out these stats
# for fungal subsets
fung_snowbrush_graph_stats <- list()
for(i in ls(pattern = "fung_sb_")){
  stat_list <- find_ig_subset_attr(ps.subset = get(i),ig.full = fung_igraph)
  fung_snowbrush_graph_stats[[i]] <- stat_list
}
fung_grandfir_graph_stats <- list()
for(i in ls(pattern = "fung_gf_")){
  stat_list <- find_ig_subset_attr(ps.subset = get(i),ig.full = fung_igraph)
  fung_grandfir_graph_stats[[i]] <- stat_list
}

# for bacterial subsets
bact_snowbrush_graph_stats <- list()
for(i in ls(pattern = "bact_sb_")){
  stat_list <- find_ig_subset_attr(ps.subset = get(i),ig.full = bact_igraph)
  bact_snowbrush_graph_stats[[i]] <- stat_list
}
bact_grandfir_graph_stats <- list()
for(i in ls(pattern = "bact_gf_")){
  stat_list <- find_ig_subset_attr(ps.subset = get(i),ig.full = bact_igraph)
  bact_grandfir_graph_stats[[i]] <- stat_list
}

# add all to single object
full_subset_graphs <- 
list(fungi=list(grandfir=fung_grandfir_graph_stats,snowbrush=fung_snowbrush_graph_stats),
     bacteria=list(grandfir=bact_grandfir_graph_stats,snowbrush=bact_snowbrush_graph_stats))


# deeply nested list... access as follows:

full_subset_graphs$fungi$grandfir$fung_gf_1$mean_dist

## Extract attributes ####

### Degree Distributions ####
# isolate degree_distributions From all subsets
f_gf <- 
data.frame(domain = "fungi",
           host = "grandfir",
           deg_dist = I(
             full_subset_graphs$fungi$grandfir %>% 
               map('deg_dist'))
           )
f_sb <- 
  data.frame(domain = "fungi",
             host = "snowbrush",
             deg_dist = I(
               full_subset_graphs$fungi$snowbrush %>% 
                 map('deg_dist'))
  )
b_gf <- 
  data.frame(domain = "bacteria",
             host = "grandfir",
             deg_dist = I(
               full_subset_graphs$bacteria$grandfir %>% 
                 map('deg_dist'))
  )
b_sb <- 
  data.frame(domain = "bacteria",
             host = "snowbrush",
             deg_dist = I(
               full_subset_graphs$bacteria$snowbrush %>% 
                 map('deg_dist'))
  )
# stick them together
deg_dist_df <- 
f_gf %>% 
  bind_rows(f_sb) %>% 
  bind_rows(b_gf) %>% 
  bind_rows(b_sb)
deg_dist_df$group <- row.names(deg_dist_df)
deg_dist_df$inoc <- deg_dist_df$group %>% str_split("_") %>% map_chr(3)

# make deg_dist lengths match
m <-
  deg_dist_df$deg_dist %>% 
  map_dbl(length) %>% 
  max
deg_dist_df$deg_dist <- lapply(deg_dist_df$deg_dist, "length<-",m)  

deg_dist_df <- 
deg_dist_df$deg_dist %>% 
  as.data.frame() %>% 
  pivot_longer(everything(),names_to = "group") %>% 
  full_join(deg_dist_df)
deg_dist_df <- 
  deg_dist_df %>% 
  mutate(fire_freq = case_when(inoc %in% 1:2 ~ "0",
                               inoc %in% 3:4 ~ "1",
                               inoc %in% 5:6 ~ "3",
                               inoc == "s" ~ "sterile"))


deg_dist_df %>% 
  ggplot(aes(x=value,color=inoc)) +
  geom_density() +
  facet_wrap(~domain*host)

# model
glm(data = deg_dist_df,
    formula = value ~ inoc * domain * host) %>% 
  summary()
  

### Max degree ####
max_degree_df <- 
c(
  full_subset_graphs$fungi$grandfir %>% 
    map_dbl(pluck("max_degree")),
  full_subset_graphs$fungi$snowbrush %>% 
    map_dbl(pluck("max_degree")),
  full_subset_graphs$bacteria$grandfir %>% 
    map_dbl(pluck("max_degree")),
  full_subset_graphs$bacteria$snowbrush %>% 
    map_dbl(pluck("max_degree"))
) %>% 
  as.data.frame()
max_degree_df$group <- row.names(max_degree_df)
names(max_degree_df)[1] <- "max_degree"

graph_atributes_df <- 
full_join(deg_dist_df,max_degree_df) %>% 
  select(-value) %>% 
  unique.data.frame()


### N vertices ####
n_vertices_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("n_vertices")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("n_vertices")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("n_vertices")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("n_vertices"))
  ) %>% 
  as.data.frame()
n_vertices_df$group <- row.names(n_vertices_df)
names(n_vertices_df)[1] <- "n_vertices"

graph_atributes_df <- 
  full_join(graph_atributes_df,n_vertices_df)
  
### N_edges ####
n_edges_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("n_edges")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("n_edges")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("n_edges")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("n_edges"))
  ) %>% 
  as.data.frame()
n_edges_df$group <- row.names(n_edges_df)
names(n_edges_df)[1] <- "n_edges"

graph_atributes_df <- 
  full_join(graph_atributes_df,n_edges_df)

### Mean_dist ####
mean_dist_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("mean_dist")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("mean_dist")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("mean_dist")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("mean_dist"))
  ) %>% 
  as.data.frame()
mean_dist_df$group <- row.names(mean_dist_df)
names(mean_dist_df)[1] <- "mean_dist"

graph_atributes_df <- 
  full_join(graph_atributes_df,mean_dist_df)

### Clique number ####
clique_num_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("clique_num")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("clique_num")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("clique_num")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("clique_num"))
  ) %>% 
  as.data.frame()
clique_num_df$group <- row.names(clique_num_df)
names(clique_num_df)[1] <- "clique_num"

graph_atributes_df <- 
  full_join(graph_atributes_df,clique_num_df)

### Mean betweenness ####
mean_betweenness_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("mean_betweenness")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("mean_betweenness")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("mean_betweenness")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("mean_betweenness"))
  ) %>% 
  as.data.frame()
mean_betweenness_df$group <- row.names(mean_betweenness_df)
names(mean_betweenness_df)[1] <- "mean_betweenness"

graph_atributes_df <- 
  full_join(graph_atributes_df,mean_betweenness_df)


### Mean closeness ####
mean_closeness_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("mean_closeness")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("mean_closeness")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("mean_closeness")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("mean_closeness"))
  ) %>% 
  as.data.frame()
mean_closeness_df$group <- row.names(mean_closeness_df)
names(mean_closeness_df)[1] <- "mean_closeness"

graph_atributes_df <- 
  full_join(graph_atributes_df,mean_closeness_df)

### Mean coreness ####
mean_coreness_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("mean_coreness")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("mean_coreness")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("mean_coreness")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("mean_coreness"))
  ) %>% 
  as.data.frame()
mean_coreness_df$group <- row.names(mean_coreness_df)
names(mean_coreness_df)[1] <- "mean_coreness"

graph_atributes_df <- 
  full_join(graph_atributes_df,mean_coreness_df)

### Global efficiency ####
global_effic_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("global_effic")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("global_effic")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("global_effic")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("global_effic"))
  ) %>% 
  as.data.frame()
global_effic_df$group <- row.names(global_effic_df)
names(global_effic_df)[1] <- "global_effic"

graph_atributes_df <- 
  full_join(graph_atributes_df,global_effic_df)

### Clustering coeficient ####
clustering_coeficient_df <- 
  c(
    full_subset_graphs$fungi$grandfir %>% 
      map_dbl(pluck("clustering_coeficient")),
    full_subset_graphs$fungi$snowbrush %>% 
      map_dbl(pluck("clustering_coeficient")),
    full_subset_graphs$bacteria$grandfir %>% 
      map_dbl(pluck("clustering_coeficient")),
    full_subset_graphs$bacteria$snowbrush %>% 
      map_dbl(pluck("clustering_coeficient"))
  ) %>% 
  as.data.frame()
clustering_coeficient_df$group <- row.names(clustering_coeficient_df)
names(clustering_coeficient_df)[1] <- "clustering_coeficient"

graph_atributes_df <- 
  full_join(graph_atributes_df,clustering_coeficient_df)

# Save data frame
saveRDS(graph_atributes_df,"./Output/network_attributes_grouped.RDS")



# PLOT ATTRIBUTES ####

# betweenness
graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=mean_betweenness,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)

# clique number 
graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=clique_num,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)
# max degree
graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=max_degree,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)
# n edges
graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=n_edges,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)

# SINGLE-SAMPLE SUBSETS ####
# Build list objects with attributes and subgraphs for each sample

## Fungal subsets ####
fungal_subsets <- list()
for(i in sample_names(fung)){
  x <- subset_samples(fung,sample_names(fung) == i)
  fungal_subsets[[i]] <- find_ig_subset_attr(ps.subset = x, ig.full = fung_igraph)
}

## Bacterial subsets ####
bacterial_subsets <- list()
for(i in sample_names(bact_genus)){
  x <- subset_samples(bact_genus,sample_names(bact_genus) == i)
  bacterial_subsets[[i]] <- find_ig_subset_attr(ps.subset = x, ig.full = bact_igraph)
}


## Build data frames ####
fungal_network_attributes_df <- 
data.frame(sample_name = names(fungal_subsets),
           n_vertices = fungal_subsets %>% map_dbl(pluck("n_vertices")),
           n_edges = fungal_subsets %>% map_dbl(pluck("n_edges")),
           mean_dist = fungal_subsets %>% map_dbl(pluck("mean_dist")),
           clique_num = fungal_subsets %>% map_dbl(pluck("clique_num")),
           mean_betweenness = fungal_subsets %>% map_dbl(pluck("mean_betweenness")),
           mean_closeness = fungal_subsets %>% map_dbl(pluck("mean_closeness")),
           mean_coreness = fungal_subsets %>% map_dbl(pluck("mean_coreness")),
           global_effic = fungal_subsets %>% map_dbl(pluck("global_effic")),
           clustering_coeficient = fungal_subsets %>% map_dbl(pluck("clustering_coeficient")),
           max_degree = fungal_subsets %>% map_dbl(pluck("max_degree"))
           ) %>% 
  full_join(fung@sam_data)

bacterial_network_attributes_df <- 
  data.frame(sample_name = names(bacterial_subsets),
             n_vertices = bacterial_subsets %>% map_dbl(pluck("n_vertices")),
             n_edges = bacterial_subsets %>% map_dbl(pluck("n_edges")),
             mean_dist = bacterial_subsets %>% map_dbl(pluck("mean_dist")),
             clique_num = bacterial_subsets %>% map_dbl(pluck("clique_num")),
             mean_betweenness = bacterial_subsets %>% map_dbl(pluck("mean_betweenness")),
             mean_closeness = bacterial_subsets %>% map_dbl(pluck("mean_closeness")),
             mean_coreness = bacterial_subsets %>% map_dbl(pluck("mean_coreness")),
             global_effic = bacterial_subsets %>% map_dbl(pluck("global_effic")),
             clustering_coeficient = bacterial_subsets %>% map_dbl(pluck("clustering_coeficient")),
             max_degree = bacterial_subsets %>% map_dbl(pluck("max_degree"))
  ) %>% 
  full_join(bact_genus@sam_data)

full_network_attributes_df <- 
bacterial_network_attributes_df  %>%
  mutate(block = as.character(block)) %>% 
  mutate(inoculum_burn_freq = ordered(inoculum_burn_freq,levels = c("0","1","3"))) %>% 
  bind_rows(fungal_network_attributes_df) %>% 
  mutate(kingdom = case_when(amplicon == "16S" ~ "Bacteria",
                             amplicon == "ITS2" ~ "Fungi"))


## Plot attributes ####
net_attrib_names <- c("n_vertices","n_edges","mean_dist","clique_num",
                      "mean_betweenness","mean_closeness","mean_coreness",
                      "global_effic","clustering_coeficient","max_degree")

attribute_plots <- list()
for(i in net_attrib_names){
  
  p <- full_network_attributes_df %>% 
    ggplot(aes(x=inoculum_site,
               y= .data[[i]],
               fill=species)) +
    geom_boxplot() +
    facet_wrap(~kingdom) +
    scale_fill_manual(values = host_colors) +
    theme(axis.text.x = element_text(angle = 90,hjust = 1),
          axis.title.x = element_blank())
attribute_plots[[i]] <- p  
}

### plot all together ####
attribute_multiplot <- 
attribute_plots$n_vertices + attribute_plots$n_edges + attribute_plots$mean_dist + attribute_plots$clique_num + attribute_plots$mean_betweenness + 
  attribute_plots$mean_closeness + attribute_plots$mean_coreness + attribute_plots$global_effic + attribute_plots$clustering_coeficient + attribute_plots$max_degree +
  patchwork::plot_layout(guides = 'collect',nrow = 5) +
  patchwork::plot_annotation(title = "Community network properties")

# save plot for later (Supplementary Info)
saveRDS(attribute_multiplot,"./Output/figs/Network_Attributes_Plot_by_host-inoc-kingdom.RDS")

### Plot against plant health ####

# build composite plant health score
full_network_attributes_df_long <- 
  full_network_attributes_df %>% 
  mutate(across(c("wilting_scale","bud_number","leaf_number",
                  "leaf_length","height","shoot_dm","final_root_dm"),
                scale01)) %>% # scale/center all indicators
  pivot_longer(c("wilting_scale","bud_number","leaf_number",
                 "leaf_length","height","shoot_dm","final_root_dm"),
               names_to="indicator",values_to = "plant_health_indicator")
  

attribute_plots <- list()
for(i in net_attrib_names){
  
p <- 
  full_network_attributes_df_long %>% 
  ggplot(aes(x=.data[[i]],
             y=plant_health_indicator,
             color=species)) +
  geom_point(alpha=.3) +
  facet_wrap(~kingdom,scales = 'free_x') +
  scale_color_manual(values = host_colors) +
  geom_smooth(method='lm',se=FALSE) +
  labs(y="Plant health score")

attribute_plots[[i]] <- p  

}

### plot all together ####
attribute_multiplot <- 
  attribute_plots$n_vertices + attribute_plots$n_edges + attribute_plots$mean_dist + attribute_plots$clique_num + attribute_plots$mean_betweenness + 
  attribute_plots$mean_closeness + attribute_plots$mean_coreness + attribute_plots$global_effic + attribute_plots$clustering_coeficient + attribute_plots$max_degree +
  patchwork::plot_layout(guides = 'collect',nrow = 2) +
  patchwork::plot_annotation(title = "Community network properties")

# save plot for later (Supplementary Info)
saveRDS(attribute_multiplot,"./Output/figs/Network_Attributes_vs_Plant_Health.RDS")



## Model ####
mod_full <- 
  glm(data = full_network_attributes_df_long,
      formula = plant_health_indicator ~ species * kingdom *
        n_vertices + n_edges + mean_dist + mean_betweenness + mean_closeness + global_effic + clustering_coeficient + max_degree)
mod_full %>% summary

saveRDS(mod_full,"./Output/Model_Network_Attributes_vs_Plant_Health.RDS")


# PLOT NETWORKS ####

# isolate network subsets and add to data frame as list-column

full_network_attributes_df$igraph <- 
c(fungal_subsets %>% 
  map(pluck("ig")),
  bacterial_subsets %>% 
    map(pluck("ig")))

# Save all plots to external image files

# remove empty rows
full_network_attributes_df[full_network_attributes_df$sample_name == "F-R-081",'igraph']
dev.off()
for(i in 1:nrow(full_network_attributes_df)){
  ig <- full_network_attributes_df$igraph[[i]]
  p.title <- paste0(full_network_attributes_df$species[i], " - ",
                    full_network_attributes_df$kingdom[i], " - ",
                    "Inoc: ",full_network_attributes_df$inoculum_site[i], " - ",
                    full_network_attributes_df$sample_name[i])
  
  
  p.title_2 <- p.title %>% str_replace_all(" - ","_") %>% str_replace_all(": ","-")
  path <- file.path("./Output/figs/igraphs",paste0(p.title_2,".png"))
  
  if(length(igraph::V(ig)) < 10){ # some networks are essentially empty (F-R-081)
    next
  }
  png(path)
  plot_hubs(ig)
  title(p.title)
  dev.off()
}
  




# HUB TAXA ####
bact_hub_taxa
fung_hub_taxa

# Look at (relative) abundance of hub taxa vs plant health

# STOPPED HERE Jan 28 ###################################################


# Run Specieasi for each marker and inoculum source combination
# moved to separate script for clarity
source("./R/SpiecEasi_by_inoc_sources.R")
# CROSS-DOMAIN NETWORKS ####
# These are run in the above script as well


# identify RDS object files generated by above script
igraph_files <- list.files("./Output", full.names = TRUE, pattern = "igraph.*out.RDS")
se
se <- readRDS("./Output/ITS_SpiecEasi_4_out.RDS")
fung_igraph <- adj2igraph(getRefit(se),vertex.attr = )
readRDS("Output/ITS_SpiecEasi_1_out.RDS") %>% plot


















###########################################################

# NETWORK ANALYSIS ####
set.seed(666)
fung_network <- ps %>% 
  make_network(keep.isolates = TRUE)

# quick plot
plot_net(ps,distance = 'jaccard',color = 'inoculum_site',rescale = TRUE,maxdist = .7,point_label = "inoculum_burn_freq")

fung_network
igraph_ITS_1

data.frame(inoculum_site = c(1,2,3,4,5,6,"Sterile"),
           mean_degree = c(igraph::degree(igraph_ITS_1) %>% mean(),
                           igraph::degree(igraph_ITS_2) %>% mean(),
                           igraph::degree(igraph_ITS_3) %>% mean(),
                           igraph::degree(igraph_ITS_4) %>% mean(),
                           igraph::degree(igraph_ITS_5) %>% mean(),
                           igraph::degree(igraph_ITS_6) %>% mean(),
                           igraph::degree(igraph_ITS_Sterile) %>% mean()),
           hub_score = c(igraph::hub_score(igraph_ITS_1)$value,
                         igraph::hub_score(igraph_ITS_2)$value,
                         igraph::hub_score(igraph_ITS_3)$value,
                         igraph::hub_score(igraph_ITS_4)$value,
                         igraph::hub_score(igraph_ITS_5)$value,
                         igraph::hub_score(igraph_ITS_6)$value,
                         igraph::hub_score(igraph_ITS_Sterile)$value)) %>% 
  arrange((mean_degree)) %>% 
  right_join(fung %>% microbiome::meta(), by="inoculum_site") %>% 
  ggplot(aes(x=mean_degree,y=leaf_number)) +
  geom_smooth(method='lm',color = "#6e4618",fill="#c9c9c7") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(y="Leaf number",x="Network complexity") 
ggsave("./Output/figs/new_leaf_number_plot.png",height = 4,width = 12)








igraph::alpha_centrality(igraph_ITS_1)
igraph::assortativity_degree(fung_network,directed = FALSE)
deg <- igraph::degree(fung_network)
tmax <- igraph::centr_degree_tmax((fung_network),loops = FALSE)
igraph::centralize(deg, tmax)
igraph::cliques(fung_network)
igraph::assortativity(fung_network,
                      types1 = ps@sam_data$inoculum_site %>% factor %>% as.numeric)
igraph::largest_cliques(fung_network)


# find hub taxa (in inoc 4)

fung4 <- 
fung %>% 
  subset_samples(inoculum_site == "4") 
fung4 %>% 
  subset_taxa(taxa_sums(fung4) > 1)
igraph::betweenness(igraph_ITS_4)[igraph::betweenness(igraph_ITS_4) > 1000]
tax_table(fung4)[which(igraph::betweenness(igraph_ITS_4) > 1000)] %>% unname %>% 
  as.data.frame() %>% 
  reduce(paste)
tax_table(fung4)[which(igraph::degree(igraph_ITS_4) > 4)]



fung %>% ntaxa()
# grab response of each plant variable
leaf_mod <- fung %>% 
  microbiome::meta() %>% 
  select(leaf_number,bud_number,height,shoot_dm,final_root_dm,inoculum_site) %>% 
  glm(data=., formula = leaf_number ~ inoculum_site)
coef(leaf_mod)

log(coef(leaf_mod)[1] + coef(leaf_mod)[2]) - log(coef(leaf_mod)[1])

beepr::beep(sound = 4)

fung %>% 
  microbiome::meta() %>% 
  select(leaf_number,bud_number,height,shoot_dm,final_root_dm,inoculum_site) %>% 
  mutate(response_coef = case_when(inoculum_site == "Sterile" ~ coef(leaf_mod)[7],
                                   inoculum_site == "6" ~ coef(leaf_mod)[6],
                                   inoculum_site == "5" ~ coef(leaf_mod)[5],
                                   inoculum_site == "4" ~ coef(leaf_mod)[4],
                                   inoculum_site == "3" ~ coef(leaf_mod)[3],
                                   inoculum_site == "2" ~ coef(leaf_mod)[2],
                                   inoculum_site == "1" ~ median(leaf_number))) %>%
  mutate(scaled_leaf_number = scale(leaf_number),
         scaled_bud_number = scale(bud_number),
         scaled_height = scale(height)) %>% 
  pivot_longer(starts_with("scaled_"),names_to = "measure",values_to = "response",names_prefix = "scaled_",values_transform = as.numeric) %>% 
  # mutate(z_leaf_number = zscore(leaf_number) %>% log,
  #        z_bud_number = zscore(bud_number) %>% log,
  #        z_height = zscore(height) %>% log) %>% 
  # pivot_longer(starts_with("z_"),names_to = "measure",values_to = "log_zscore",names_prefix = "z_") %>% 
  ggplot(aes(x=factor(inoculum_site,levels=c("Sterile","3","2","5","6","1","4")),
             y=response,
             fill = measure)) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_boxplot() + # fill="#6e4618"
  # stat_summary(fun='mean',fill="#6e4618",geom = 'bar') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face='bold',size=18),panel.grid = element_blank(),
        legend.position = 'bottom',legend.text = element_text(size=12,face='bold'),
        legend.title = element_text(size=18,face='bold')) +
  labs(y="Plant response (scaled/centered)",fill="Measure") +
  lims(y=c(-1.25,2)) +
  scale_fill_manual(values = c("#6e4618","#97b031","#86808c"),labels = c("Bud number","Height","Leaf number"))

ggsave("./Output/figs/new_leaf_number_plot3.png",height = 4,width = 12)
