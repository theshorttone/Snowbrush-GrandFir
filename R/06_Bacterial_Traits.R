# -----------------------------------------------------------------------------#
# Testing trait enrichment between inoculum sources 
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.42.0
#                     patchwork v 1.1.2
#                     BacDive v 0.8.0
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
#  install.packages("BacDive", repos="http://R-Forge.R-project.org")

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")

# Random seed
set.seed(666)


# Data
bact <- readRDS("./Output/16S_clean_phyloseq_object.RDS")

# initialize bacdive connection (add your userid and password)
# BacDive URL: https://bacdive.dsmz.de/
bacdive <- BacDive::open_bacdive(username = Sys.getenv("BACDIVE_USER"),
                                 password = Sys.getenv("BACDIVE_PW"))

# get list of unique genera
genus_list <- bact@tax_table[,6] %>% 
  table()
genus_list <- genus_list %>% 
  as.data.frame() %>% 
  pluck(".") %>% 
  levels() 
names(genus_list) <- genus_list

# BacDive API ####

# run in for-loop to build giant database of all traits for our taxa
bact_trait_db <- list()

for(i in genus_list){
  x <- BacDive::request(object = bacdive,
                        query = genus_list[i],
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  # if no results found, return NA
  if(length(y) == 0){
    bact_trait_db[[i]] <- NA
  } else {
    # otherwise, get info on that id
    z <- BacDive::fetch(bacdive,y)
    bact_trait_db[[i]] <- z$results  
  }
}

# save output from previous slow step
saveRDS(bact_trait_db,"./Output/16S_Bacterial_Trait_Database.RDS")
# reload point, for convenience
bact_trait_db <- readRDS("./Output/16S_Bacterial_Trait_Database.RDS")
        
# EXPLORE ####
# This database is a highly nested and confusing pile of information
# access looks something like
bact_trait_db$Rhizobium$`132155`$`Physiology and metabolism`
# DB_OBJECT   Genus   AccessionID   Trait Aspect

bact_trait_db$Rhizobium$`132155`$`Physiology and metabolism` %>% 
  str


# Best way night be a nested for-loop with lots of tests for missing data :(

genus_physiology <- list()
for(genus in names(bact_trait_db)){
  
  all_strains <- bact_trait_db[[genus]]
  
  
  
  for(strain in names(all_strains)){
    
    strain_x <- all_strains[[strain]]
    physiology_x <- strain_x[["Physiology and metabolism"]]
    
    
    
    if(length(physiology_x) < 1){
      genus_physiology[[genus]] <- NA
    }
  }  
}

bact_trait_db$Rhizobium %>% pluck(1) %>% names

