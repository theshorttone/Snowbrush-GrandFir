# -----------------------------------------------------------------------------#
# SpiecEasi plots for both bacteria and fungi by inoculum source
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
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

plot_hubs <- function(graph,bigpoint=10,littlepoint=3){
  am.coord <- layout.auto(graph)
  art <- articulation_points(graph)
  plot(graph,
       layout = am.coord,
       vertex.size=(scale01(abs(igraph::authority_score(graph)$vector)) * 10)+3,
       vertex.label=NA,main=paste0(marker,"_inoc_",inoc))
}

scale01 <- function(x){
  if(max(x) == 0){return(x)}
  if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))}
}

# Data
bact <- readRDS("./Output/16S_clean_phyloseq_object.RDS") %>% 
  subset_samples(species == "GrandFir") 
bact <- bact %>% subset_taxa(taxa_sums(bact) >= 100)

fung <- readRDS("./Output/ITS_clean_phyloseq_object.RDS") %>% 
  subset_samples(species == "GrandFir")
fung <- fung %>% subset_taxa(taxa_sums(fung) > 0)


# Run for fungi
for(i in unique(fung@sam_data$inoculum_site)){
  
  ps <- fung
  inoc <- i
  marker="ITS"
  ncores <- 20
  
  set.seed(666)
  
  # set pulsar parameters for SpiecEasi
  se.params <- list(rep.num=20, ncores=ncores, seed=666)
  
  # subset to a given inoculum source
  ps_sub <- ps %>% 
    subset_samples(inoculum_site == inoc)
  # remove empty ASVs
  ps_sub <- ps_sub %>% 
    subset_taxa(taxa_sums(ps_sub) > 0)
  ps_sub
  # run spieceasi
  se <- SpiecEasi::spiec.easi(data = ps_sub,
                              method='mb',
                              sel.criterion = "bstars",
                              pulsar.params=se.params,
                              verbose=TRUE)
  # build file name
  fn <- paste0("./Output/",marker,"_SpiecEasi_",inoc,"_out.RDS")
  # save spieceasi object for later recall
  saveRDS(se, fn)
  
  # get best model and build igraph
  se_igraph <- adj2igraph(getRefit(se), vertex.attr = list(name=NA))
  # save that, as well
  fn2 <- paste0("./Output/",marker,"_igraph_",inoc,"_out.RDS")
  saveRDS(se_igraph, fn2)
  assign(paste0("igraph_",marker,"_",inoc),se_igraph,envir = .GlobalEnv)
  
  # Plot with igraph
  ## set size of vertex proportional to sum relabund
  vsize    <- transform_sample_counts(ps_sub,function(x){x/sum(x)}) %>% taxa_sums() + 3
  am.coord <- layout.auto(se_igraph)
  png(filename = paste0("./Output/figs/igraph_",marker,"_",inoc,".png"),width = 4,height = 4,res = 200,units = "in")
  plot_hubs(se_igraph)
  dev.off()
}


# look at connectivity degree distributions for fungi
degree_list <- 
  list(igraph_ITS_1,igraph_ITS_2,igraph_ITS_3,igraph_ITS_4,igraph_ITS_5,igraph_ITS_6,igraph_ITS_Sterile) %>% 
  map(degree_distribution)
degree_df <- lapply(degree_list, "length<-", max(lengths(degree_list))) %>% 
  as.data.frame()
names(degree_df) <- c("Inoc_1","Inoc_2","Inoc_3","Inoc_4","Inoc_5","Inoc_6","Inoc_Sterile")

# convert NA to 0 for regression plotting
degree_df2 <- degree_df
degree_df2[is.na(degree_df2)] <- 0 
c(unlist(degree_df2)) %>% length

# make better network figure for proposal

df <- 
data.frame(
  inoc_source=fung@sam_data$inoculum_site,
  leaf_number=fung@sam_data$leaf_number,
  bud_count=fung@sam_data$bud_number,
  height=fung@sam_data$height,
  shoot_mass=fung@sam_data$shoot_dm,
  root_mass=fung@sam_data$final_root_dm
)
zscore <- function(x){(x-mean(x,na.rm=TRUE)/sd(x,na.rm=TRUE))}

df_long <- 
bind_cols(apply(df[,-1],2,zscore),inoc_source=df$inoc_source) %>% 
  pivot_longer(-inoc_source,
               names_to = "measure",
               values_to = "zscore")


deglist <- 
  c("igraph_ITS_1",'igraph_ITS_2','igraph_ITS_3','igraph_ITS_4','igraph_ITS_5','igraph_ITS_6',"igraph_ITS_Sterile")
x <- c()
for(i in deglist){
  x[i] <- igraph::degree((get(i)),loops = FALSE) %>% mean(na.rm=FALSE)
}

x
df_long %>% 
  mutate(connectivity = 
           case_when(inoc_source == "1" ~ x[1],
                     inoc_source == "2" ~ x[2],
                     inoc_source == "3" ~ x[3],
                     inoc_source == "4" ~ x[4],
                     inoc_source == "5" ~ x[5],
                     inoc_source == "6" ~ x[6],
                     inoc_source == "Sterile" ~ x[7])) %>% 
  dplyr::filter(measure == "leaf_number") %>% 
  mutate(measure=case_when(measure=="height" ~ "Height",
                           measure=="leaf_number" ~ "Leaf count",
                           measure=="root_mass" ~ "Root mass",
                           measure=="shoot_mass" ~ "Shoot mass")) %>% 
  ggplot(aes(x=connectivity,y=zscore)) +
  # geom_point(alpha=.5) +
  geom_smooth(method='lm',color='black',fill='lightgray') +
  labs(x="Inoculum network connectivity",y="Leaf number") +
  # facet_wrap(~measure,scales='free') +
  theme_minimal() +
  theme(strip.text = element_text(face='bold',size=16),
        axis.text.x = element_blank(),
        axis.title = element_text(face='bold',size = 16),panel.grid = element_blank()
        )
ggsave("Output/figs/fungal_network_connectivity_vs_leaf_number_noaxis.jpg",width = 6,height = 6,dpi=300)


(
fplot <- 
degree_df %>% 
  mutate(degree=1:nrow(.)) %>% 
  pivot_longer(starts_with("Inoc")) %>% 
  ggplot(aes(x=degree,y=value,color=name)) +
  geom_point() +
  geom_path() +
  theme_minimal() +
  labs(x="Degree of connectivity",y="Frequency",title = "Fungal community connectivity",color="Inoculum source")
)
saveRDS(fplot,"./Output/figs/fungal_network_connectivity.RDS")

# Run for bacteria
for(i in unique(bact@sam_data$inoculum_site)){
  
  ps <- bact
  inoc <- i
  marker="V6V8"
  ncores <- 24
  
  set.seed(666)
  
  # set pulsar parameters for SpiecEasi
  se.params <- list(rep.num=20, ncores=ncores, seed=666)
  
  # subset to a given inoculum source
  ps_sub <- ps %>% 
    subset_samples(inoculum_site == inoc)
  # remove empty ASVs
  ps_sub <- ps_sub %>% 
    subset_taxa(taxa_sums(ps_sub) > 50)
  ps_sub

  # run spieceasi
  se <- SpiecEasi::spiec.easi(data = ps_sub,
                              method='mb',
                              sel.criterion = "bstars",
                              pulsar.params=se.params,
                              verbose=TRUE)
  # build file name
  fn <- paste0("./Output/",marker,"_SpiecEasi_",inoc,"_out.RDS")
  # save spieceasi object for later recall
  saveRDS(se, fn)

  # get best model and build igraph
  se_igraph <- adj2igraph(getRefit(se), vertex.attr = list(name=NA))
  # save that, as well
  fn2 <- paste0("./Output/",marker,"_igraph_",inoc,"_out.RDS")
  saveRDS(se_igraph, fn2)
  assign(paste0("igraph_",marker,"_",inoc),se_igraph,envir = .GlobalEnv)
  layout_
  # Plot with igraph
  ## set size of vertex proportional to sum relabund
  vsize    <- transform_sample_counts(ps_sub,function(x){x/sum(x)}) %>% taxa_sums() + 3
  am.coord <- layout.auto(se_igraph)
  png(filename = paste0("./Output/figs/igraph_",marker,"_",inoc,".png"),width = 4,height = 4,res = 200,units = "in")
  plot(se_igraph, layout=am.coord, vertex.size=vsize, vertex.label=NA,main=paste0(marker,"_inoc_",inoc))
  dev.off()
}

# look at connectivity degree distributions for bacteria
degree_list <- 
  list(igraph_V6V8_1,igraph_V6V8_2,igraph_V6V8_3,igraph_V6V8_4,igraph_V6V8_5,igraph_V6V8_6,igraph_V6V8_Sterile) %>% 
  map(degree_distribution)
degree_df <- lapply(degree_list, "length<-", max(lengths(degree_list))) %>% 
  as.data.frame()
names(degree_df) <- c("Inoc_1","Inoc_2","Inoc_3","Inoc_4","Inoc_5","Inoc_6","Inoc_Sterile")
(
  fplot <- 
    degree_df %>% 
    mutate(degree=1:nrow(.)) %>% 
    pivot_longer(starts_with("Inoc")) %>% 
    ggplot(aes(x=degree,y=value,color=name)) +
    geom_point() +
    geom_path() +
    theme_minimal() +
    labs(x="Degree of connectivity",y="Frequency",title = "Bacterial community connectivity",color="Inoculum source")
)
saveRDS(fplot,"./Output/figs/bacterial_network_connectivity.RDS")



# Run for cross-domain
for(i in unique(bact@sam_data$inoculum_site)){
  
  ps1 <- bact
  ps2 <- fung
  inoc <- i
  marker="cross-domain"
  ncores <- 24
  
  set.seed(666)
  
  # set pulsar parameters for SpiecEasi
  se.params <- list(rep.num=20, ncores=ncores, seed=666)
  
  # subset to a given inoculum source
  ps_sub1 <- ps1 %>% 
    subset_samples(inoculum_site == inoc)
  ps_sub2 <- ps2 %>% 
    subset_samples(inoculum_site == inoc)
  
  # remove empty ASVs
  ps_sub1 <- ps_sub1 %>% 
    subset_taxa(taxa_sums(ps_sub1) > 50)
  ps_sub2 <- ps_sub2 %>% 
    subset_taxa(taxa_sums(ps_sub2) > 0)
  # make sample names match
  sample_names(ps_sub2) <- sample_names(ps_sub2) %>% str_remove("F-")
  ps_sub2 %>% sample_names()
  ps_sub1 %>% sample_names()
  
  # run spieceasi
  se <- SpiecEasi::spiec.easi(list(ps_sub1,ps_sub2),
                              method='mb',
                              sel.criterion = "bstars",
                              pulsar.params=se.params,
                              verbose=TRUE)
  # build file name
  fn <- paste0("./Output/",marker,"_SpiecEasi_",inoc,"_out.RDS")
  # save spieceasi object for later recall
  saveRDS(se, fn)
  
  # get best model and build igraph
  se_igraph <- adj2igraph(getRefit(se), vertex.attr = list(name=NA))
  # save that, as well
  fn2 <- paste0("./Output/",marker,"_igraph_",inoc,"_out.RDS")
  saveRDS(se_igraph, fn2)
  assign(paste0("igraph_",marker,"_",inoc),se_igraph,envir = .GlobalEnv)
  
  # Plot with igraph
  ## set size of vertex proportional to sum relabund
  hmp216S
  ps <- merge_phyloseq(ps_sub1,ps_sub2)
  taxa_sums(ps) %>% summary
  vsize    <- transform_sample_counts(ps,
                                      function(x){x/sum(x)}) %>% taxa_sums() %>% sqrt() + 4
  
  vsize %>% summary
  am.coord <- layout_with_mds(se_igraph)
  grouping <- c(rep('blue',ntaxa(ps_sub1)), rep('orange',ntaxa(ps_sub2)))
  png(filename = paste0("./Output/figs/igraph_",marker,"_",inoc,".png"),width = 4,height = 4,res = 200,units = "in")
  plot(se_igraph, layout=am.coord, vertex.size=vsize, 
       vertex.label=NA, main=paste0(marker," inoc_source: ",inoc),
       vertex.color=grouping)
  dev.off()
}

# look at connectivity degree distributions for both domains
degree_list <- 
  list(`igraph_cross-domain_1`,`igraph_cross-domain_2`,`igraph_cross-domain_3`,
       `igraph_cross-domain_4`,`igraph_cross-domain_5`,`igraph_cross-domain_6`,`igraph_cross-domain_Sterile`) %>% 
  map(degree_distribution)
degree_df <- lapply(degree_list, "length<-", max(lengths(degree_list))) %>% 
  as.data.frame()
names(degree_df) <- c("Inoc_1","Inoc_2","Inoc_3","Inoc_4","Inoc_5","Inoc_6","Inoc_Sterile")
(
  fplot <- 
    degree_df %>% 
    mutate(degree=1:nrow(.)) %>% 
    pivot_longer(starts_with("Inoc")) %>% 
    ggplot(aes(x=degree,y=value,color=name)) +
    geom_point() +
    geom_path() +
    theme_minimal() +
    labs(x="Degree of connectivity",y="Frequency",title = "Full community connectivity",color="Inoculum source")
)
saveRDS(fplot,"./Output/figs/cross-domain_network_connectivity.RDS")


print.igraph(`igraph_cross-domain_1`)


plot.igraph(igraph_V6V8_1)
plot.igraph(igraph_V6V8_2)
plot.igraph(igraph_V6V8_3)
plot.igraph(igraph_V6V8_4)
plot.igraph(igraph_V6V8_5)
plot.igraph(igraph_V6V8_6)
plot.igraph(igraph_V6V8_Sterile)
