# -----------------------------------------------------------------------------#
# 18S DADA2 pipeline
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
#                     ShortRead v 1.54.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(readxl); packageVersion("readxl")
library(tidyverse); packageVersion("tidyverse")
library(parallel); packageVersion("parallel")
library(ShortRead);packageVersion("ShortRead")
nthreads <- parallel::detectCores() - 1
set.seed(123) # "random" seed for reproducibility

condense_ps_to_species <- function(x){
  # remove bad taxa assignments
  # x <- x %>% subset_taxa(!is.na(Phylum))
  # build data frame
  sp <- x@tax_table[,9] %>% as.character()
  gn <- x@tax_table[,8] %>% as.character()
  or <- x@tax_table[,7] %>% as.character()
  sb <- x@tax_table[,6] %>% as.character()
  ph <- x@tax_table[,5] %>% as.character()
  
  condensed_taxonomy <- 
    data.frame(ph,sb,or,gn,sp) %>% 
    mutate(spp = case_when(is.na(or) & is.na(gn) & is.na(sp) ~ paste0(sb," sp."),
                           !is.na(or) & is.na(gn) & is.na(sp) ~ paste0(or," sp."),
                           !is.na(or) & !is.na(gn) & is.na(sp) ~ paste0(gn," sp."),
                           !is.na(or) & !is.na(gn) & !is.na(sp) ~ sp))
  
  x@tax_table[,9] <- str_replace_all(condensed_taxonomy$spp,"_"," ")
  x <- tax_glom(x,"Species")
  return(x)
}

# PARSE FILE PATHS ####
meta <- read_csv("./Data/Sample_Metadata_Full.csv") %>% 
  dplyr::filter(amplicon == "18S")

# File parsing - 
path <- "./Data/Raw/18S" # CHANGE to the directory containing your adaptor-free demultiplexed fastq files when using your own data
filtpath <- "./Data/Filtered/18S/Stephanie" # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads

full_file_list <- list.files(path,full.names = TRUE,pattern = "^1")
fns <- grep("_R1_001",full_file_list,value = TRUE)
rns <- grep("_R2_001",full_file_list,value = TRUE)
sample.names <- basename(fns) %>% str_split("_") %>% map_chr(1) %>% paste0("AMF-",.)

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[1:2]) + ggtitle("Example forward reads")
p2 <- plotQualityProfile(rns[1:2]) + ggtitle("Example reverse reads")

# display and save the plots
p1 / p2

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path("./Data/Filtered/18S/Stephanie", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path("./Data/Filtered/18S/Stephanie", paste0(sample.names, "_REV_filt.fastq.gz"))

# this is the actual qualit control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, rns, filts_r, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=c(2,2), # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     # truncLen = c(300,275), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=nthreads) # On Windows set multithread=FALSE

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(rns[1:2]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_r[1:2])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4

# LEARN ERROR RATES ####

# learn errors
errF <- learnErrors(filts_f, multithread=nthreads, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows
errR <- learnErrors(filts_r, multithread=nthreads, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows


# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
plotErrors(errR, nominalQ=FALSE)

###########################################################
# IF DATA WON'T FIT IN MEMORY
# Can do derep and dada one sample at a time in a for-loop
names(filts_f) <- sample.names
names(filts_r) <- sample.names
dada_f <- vector("list",length(sample.names))
names(dada_f) <- sample.names
dada_r <- vector("list",length(sample.names))
names(dada_r) <- sample.names
mergers <- vector("list",length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names){
  cat("Processing: ", sam, "\n")
  derep_f <- derepFastq(filts_f[[sam]])
  derep_r <- derepFastq(filts_r[[sam]])
  dada_f[[sam]] <- dada(derep_f,err = errF,multithread=TRUE)
  dada_r[[sam]] <- dada(derep_r,err = errR,multithread=TRUE)
  class(dada_f)
  merger <- mergePairs(dada_f[[sam]],derep_f,dada_r[[sam]],derep_r,trimOverhang = TRUE,minOverlap = 10,maxMismatch = 1)
  mergers[[sam]] <- merger
}

##############################################################


# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out <- out[as.data.frame(out)$reads.out > 0,]
# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dada_f, getN), sapply(dada_r, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names %>% str_split("_") %>% map_chr(1)
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)


# IMPORT METADATA ####

# import and clean
# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]


# ASSIGN TAXONOMY ####

# Use RDP training set for 18S
taxa <- assignTaxonomy(seqtab.nochim, 
                       "./Taxonomy/combined_pr2-maarjam_reference.fasta", 
                       multithread=nthreads,
                       taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
                       tryRC = TRUE,
                       verbose = TRUE)

# Save intermediate taxonomy file

# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- data.frame(sample.names,row.names = sample.names) %>% sample_data()
sample_names(met) <- row.names(met)

ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./Output/18S_ps_not-cleaned_Stephanie.RDS")
ps
rank_names(ps)


ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
plot_bar2(fill="Division")
ggsave("./Output/Stephanie)barplot_Division.jpg")

ps %>% 
  subset_taxa(Class == "Fungi") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Species")
ggsave("./Output/Stephanie)barplot_Division.jpg")

ps %>% 
  psmelt() %>% 
  write_csv("./Output/Stephanie_full_results.csv")
