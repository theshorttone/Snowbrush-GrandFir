# ------------------------------------------------------------------------------------------------#
# Extract ITS1 region for a given SRA project
# Apply this script to each SRA Accession directory on FWD reads only
# Author: Geoffrey Zahn
# Requirements: itsxpress (available in miniconda3)
# ------------------------------------------------------------------------------------------------#

# load packages
library(tidyverse)
library(parallel)

# get file paths (fwd only)
datapath <- "./Data/Raw/ITS"
files <- list.files(datapath,full.names = TRUE,
           pattern = "R1_001.fastq.gz") %>% 
  grep(pattern="F-R-",x=.,value=TRUE)

nthreads <- parallel::detectCores() - 1

fwd %>% 


for(fwd in files){
  out <- fwd %>% 
    str_replace("001.fastq.gz","ITS2.fastq.gz")
  arguments <- paste0("--fastq ",fwd, " --region ITS2 --taxa Fungi --outfile ", out, " --threads ",nthreads)
  
  # run itsxpress
  system2("itsxpress",args = arguments)
}
