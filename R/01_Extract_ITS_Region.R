# -----------------------------------------------------------------------------#
# Remove adaptors using cutadapt
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     parallel v 4.3.2
#                     itsxpress v 2.0.1
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
