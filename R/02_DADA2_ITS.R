# -----------------------------------------------------------------------------#
# 16S DADA2 pipeline
# Processing raw amplicon reads
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     dada2 v 1.24.0
#                     decontam v 1.16.0
#                     phyloseq v 1.40.0
#                     Biostrings v 2.64.0
#                     patchwork v 1.1.1
#                     readxl v 1.4.1
#                     janitor::clean_names() v 2.1.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(readxl); packageVersion("readxl")
library(parallel); packageVersion("parallel")

set.seed(123) # "random" seed for reproducibility
nthreads <- parallel::detectCores() - 1

# PARSE FILE PATHS ####

# File parsing - 
path <- "./Data/Raw/ITS" # CHANGE to the directory containing your adaptor-free demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads

# using only fwd reads for fungal ITS2
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_001.ITS2.fastq.gz")) # make pattern match your FWD reads
sample.names <- basename(fns) %>% str_split("_") %>% map_chr(1)

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[1]) + ggtitle("Example forward reads")

# display and save the plots
p1
ggsave("./Output/ITS_unfiltered_quality_plots.png",dpi=300,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))

# this is the actual qualit control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     # maxEE=c(3), # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     # truncLen = c(200), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=nthreads) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./Output/ITS_trackreads.RDS")
# had to truncate at 150 BP because of bad problem with quality at 3' ends

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[1]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[1])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4
ggsave("./Output/ITS_filtered_quality_comparison.png",dpi=300,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows

saveRDS(errF,"./Output/ITS_errF.RDS")
errF <- readRDS("./Output/ITS_errF.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./Output/ITS_error_model.png",dpi=300,height = 6,width = 6)

out

###########################################################
# Run DADA2 algorithm

derep_f <- derepFastq(filts_f)
dada_f <- dada(derep_f,err = errF,multithread=TRUE,selfConsist = TRUE,pool = "pseudo")
saveRDS(dada_f,"Output/ITS_dadaFs.RDS")


# Make sequence table
seqtab <- makeSequenceTable(dada_f)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./Output/ITS_seqtab.nochim.RDS")
seqtab.nochim <- readRDS("./Output/ITS_seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
# out <- readRDS("./Output/ITS_trackreads.RDS")
out <- out[as.data.frame(out)$reads.out > 0,]
# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_f, getN),  rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF",  "nonchim")
rownames(track) <- sample.names %>% str_split("_") %>% map_chr(1)
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)
track %>% 
  ggplot(aes(x=total.loss.proportion)) + geom_density()
write.csv(track, file = "./Output/ITS_read_counts_at_each_step.csv", row.names = TRUE)

# IMPORT METADATA ####

# IMPORT METADATA ####

# import and clean
meta <- read_csv("./Data/Sample_Metadata_Full.csv") %>% 
  janitor::clean_names() %>% 
  dplyr::filter(amplicon == "ITS2") # just bacteria samples for this
# subset to match seq table sample names
meta <- meta[meta$sample_name %in% (sample.names %>% str_split("_") %>% map_chr(1)), ]
row.names(seqtab.nochim) <- row.names(seqtab.nochim) %>% str_split("_") %>% map_chr(1)
row.names(meta) <- meta$sample_name

negative_ctls <- row.names(seqtab.nochim[which(!row.names(seqtab.nochim) %in% row.names(meta)),])
row.names(seqtab.nochim)
row.names(meta)
seqtab.nochim <- seqtab.nochim[!row.names(seqtab.nochim) %in% negative_ctls,]

# reorder metadata to match seqtab
df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                 sample_name=row.names(seqtab.nochim))
df2 <- left_join(meta,df,by="sample_name")
row.names(df2) <- df2$sample_name
row.names(meta) <- meta$sample_name
meta <- meta[row.names(seqtab.nochim),]
row.names(meta) <- meta$sample_name

identical(row.names(meta),row.names(seqtab.nochim))


# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./Output/ITS_seqtab.nochim.clean.RDS")

# ASSIGN TAXONOMY ####

# Use RDP training set for 16S
taxa <- assignTaxonomy(seqtab.nochim, "./Taxonomy/sh_general_release_dynamic_all_18.07.2023_dev.fasta.gz", multithread=nthreads,tryRC = TRUE,verbose = TRUE)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./Output/ITS_RDP_Taxonomy_from_dada2.RDS")

# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
sample_names(met) <- meta$sample_name

ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./Output/ITS_ps_not-cleaned.RDS")
ps
