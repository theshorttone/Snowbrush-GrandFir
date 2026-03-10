library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(corncob); packageVersion("corncob")
library(patchwork); packageVersion("patchwork")
library(igraph); packageVersion("igraph")
library(SpiecEasi); packageVersion("SpiecEasi") # >devtools::install_github("zdk123/SpiecEasi")
library(zahntools); packageVersion("zahntools") # devtools::install_github("gzahn/zahntools")
library(parallel); packageVersion("parallel")
library(hubfindr); packageVersion("hubfindr") # devtools::install_github("gzahn/hubfindr")
library(ggthemes); packageVersion("ggthemes")

source("./R/palettes.R")

host_colors <- pal.discrete[c(7,10)] 

graph_atributes_df <- readRDS("./Output/network_attributes_grouped.RDS")

se.mb.fung <- readRDS("./Output/ITS_SpiecEasi_out.RDS")

fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS")

all_plants_fungi_stats <- readRDS("./Output/all_plants_fungi_network_stats.RDS")

fung_ig <- all_plants_fungi_stats %>% map('ig')



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

fung_igraph <- adj2igraph(getRefit(se.mb.fung),vertex.attr = list(name=taxa_names(fung)))



# PLOT ATTRIBUTES ####

# betweenness
p1 <- graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=mean_betweenness,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)

# clique number 
p2 <- graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=clique_num,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)
# max degree
p3 <- graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=max_degree,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)
# n edges
p4 <- graph_atributes_df %>% 
  ggplot(aes(x=inoc,y=n_edges,fill=host)) +
  geom_col(position = 'dodge') +
  facet_wrap(~domain) +
  scale_fill_manual(values = host_colors)

(p1 + p2) / (p3 + p4) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave("./Output/figs/network_attributes_borplots_by_host.png",width = 8,height = 8)

# SINGLE-SAMPLE SUBSETS ####
# Build list objects with attributes and subgraphs for each sample

## Fungal subsets ####
fungal_subsets <- list()
for(i in sample_names(fung)){
  x <- subset_samples(fung,sample_names(fung) == i)
  fungal_subsets[[i]] <- find_ig_subset_attr(ps.subset = x, ig.full = fung_igraph)
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
fungal_network_attributes_df$block <- as.character(fungal_network_attributes_df$block)
fungal_network_attributes_df$inoculum_burn_freq <- ordered(fungal_network_attributes_df$inoculum_burn_freq,levels = c("0","1","3"))


atributes <- c("n_edges","mean_dist","clique_num",
               "mean_betweenness","mean_closeness","mean_coreness",
               "global_effic","clustering_coeficient","max_degree") 

for (a in atributes){
  p <- fungal_network_attributes_df %>%
    ggplot(aes(x = n_vertices, y = .data[[a]])) +
    geom_point() +
    geom_smooth(aes(color = host), method = "lm", se = TRUE) +
    labs(title = a, y = a)
  
  print(p)
  #ggsave(paste0("./Output/figs/network_attributes_borplots_by_host_",a,".png"),width = 8,height = 8)
}

df_long <- fungal_network_attributes_df %>%
  pivot_longer(
    cols = all_of(atributes),
    names_to = "attribute",
    values_to = "value"
  )

p <- ggplot(df_long, aes(x = n_vertices, y = value)) +
  geom_point(alpha = 0.7) +
  geom_smooth(aes(color = host), method = "lm", se = TRUE) +
  scale_color_manual(values = host_colors) +
  facet_wrap(~ attribute, scales = "free_y") +
  labs(x = "n_vertices", y = NULL)+
  theme_few()

p

