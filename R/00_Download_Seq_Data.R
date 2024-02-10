
#------------------------------------------------------------------#
# Download raw sequence data from SRA Project PRJNA1064977         #
#                                                                  #
# Requires that sra-toolkit be installed and in your system path   #
# (https://github.com/ncbi/sra-tools)                              #
#------------------------------------------------------------------#


# Setup ####
library(tidyverse)
df <- read_csv("./Data/Sample_Metadata_Full.csv")
fung <- df %>% 
  dplyr::filter(amplicon == "ITS2")
bact <- df %>% 
  dplyr::filter(amplicon == "16S")
amf <- df %>% 
  dplyr::filter(amplicon == "18S")

# Download fungal amplicons ####
accessions <- fung$sra_accession
filenames_f <- file.path(str_remove(fung$filename_r1,".gz"))
filenames_r <- file.path(str_remove(fung$filename_r2,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")
dl_filenames_2 <- paste0(accessions,"_2.fastq")


# run download in a for-loop for all accessions
# compress and rename and move files
if(!dir.exists("./Data/Raw/ITS")){
  dir.create("./Data/Raw/ITS",recursive = TRUE)
}

for(i in seq_along(accessions)){
  system2("fasterq-dump",args = accessions[i])
  system2("gzip", args = c(dl_filenames_1[i],dl_filenames_2[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),paste0("./Data/Raw/ITS/",fung$filename_r1[i]))
  file.rename(paste0(dl_filenames_2[i],".gz"),paste0("./Data/Raw/ITS/",fung$filename_r2[i]))
}



# Download bacterial amplicons ####
accessions <- bact$sra_accession
filenames_f <- file.path(str_remove(bact$filename_r1,".gz"))
filenames_r <- file.path(str_remove(bact$filename_r2,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")
dl_filenames_2 <- paste0(accessions,"_2.fastq")


# run download in a for-loop for all accessions
# compress and rename and move files
if(!dir.exists("./Data/Raw/16S")){
  dir.create("./Data/Raw/16S",recursive = TRUE)
}

for(i in seq_along(accessions)){
  system2("fasterq-dump",args = accessions[i])
  system2("gzip", args = c(dl_filenames_1[i],dl_filenames_2[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),paste0("./Data/Raw/16S/",bact$filename_r1[i]))
  file.rename(paste0(dl_filenames_2[i],".gz"),paste0("./Data/Raw/16S/",bact$filename_r2[i]))
}


# Download amf amplicons ####
accessions <- amf$sra_accession
filenames_f <- file.path(str_remove(fung$filename_r1,".gz"))
filenames_r <- file.path(str_remove(fung$filename_r2,".gz"))
dl_filenames_1 <- paste0(accessions,"_1.fastq")
dl_filenames_2 <- paste0(accessions,"_2.fastq")


# run download in a for-loop for all accessions
# compress and rename and move files
if(!dir.exists("./Data/Raw/18S")){
  dir.create("./Data/Raw/18S",recursive = TRUE)
}

for(i in seq_along(accessions)){
  system2("fasterq-dump",args = accessions[i])
  system2("gzip", args = c(dl_filenames_1[i],dl_filenames_2[i]))
  file.rename(paste0(dl_filenames_1[i],".gz"),paste0("./Data/Raw/18S/",fung$filename_r1[i]))
  file.rename(paste0(dl_filenames_2[i],".gz"),paste0("./Data/Raw/18S/",fung$filename_r2[i]))
}
