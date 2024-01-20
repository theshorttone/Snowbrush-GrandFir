# inoc="1"
# ps=fung
# marker="ITS"
# ncores=20
make_specieasi_by_inoc <- 
function(ps=fung,inoc="1",marker="ITS",ncores=(parallel::detectCores()-1)){

  set.seed(666)
  
  # set pulsar parameters for SpiecEasi
  se.params <- list(rep.num=20, ncores=ncores)
  
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
                                      pulsar.params=se.params)
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
  plot(se_igraph, layout=am.coord, vertex.size=vsize, vertex.label=NA,main=paste0(marker,"_inoc_",inoc))
  dev.off()
  
  
}


