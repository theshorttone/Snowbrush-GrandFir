# download data

library(tidyverse)
df <- read_csv("./Data/Sample_Metadata_Full.csv")

accessions <- df$sra_accession
filenames_f <- file.path("./Data/Raw",df$filename_r1)
filenames_r <- file.path("./Data/Raw",df$filename_r2)

file.path("./Data/Raw",filenames_f)

# run download in a for-loop for all accessions
for(i in seq_along(accessions)){
  x <- paste0("--split-files -o ",filenames_f[i]," ",filenames_r[i]," --skip-technical ",accessions[i])
  system2("fasterq-dump",args = x)
}