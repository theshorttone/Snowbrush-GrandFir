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

# PARSE FILE PATHS ####
meta <- read_csv("./Data/Sample_Metadata_Full.csv") %>% 
  dplyr::filter(amplicon == "18S")

# File parsing - 
path <- "./Data/Raw/18S" # CHANGE to the directory containing your adaptor-free demultiplexed fastq files when using your own data
filtpath <- "./Data/Filtered/18S" # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads

full_file_list <- list.files(path,full.names = TRUE)[list.files(path) %in% meta$filename_r1 | list.files(path) %in% meta$filename_r2]
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
ggsave("./Output/figs/AMF_unfiltered_quality_plots.png",dpi=300,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path("./Data/Filtered/18S", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path("./Data/Filtered/18S", paste0(sample.names, "_REV_filt.fastq.gz"))

# this is the actual quality control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, rns, filts_r, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=c(2,2), # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     truncLen = c(300,275), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=nthreads) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./Output/18S_trackreads.RDS")
# out <- readRDS("./Output/18S_trackreads.RDS")

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(rns[1:2]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_r[1:2])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4
ggsave("./Output/figs/18S_filtered_quality_comparison.png",dpi=300,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
errF <- learnErrors(filts_f, multithread=nthreads, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows
errR <- learnErrors(filts_r, multithread=nthreads, MAX_CONSIST = 20,verbose = 1) # set multithread = FALSE on Windows

saveRDS(errF,"./Output/18S_errF.RDS")
saveRDS(errR,"./Output/18S_errR.RDS")
errF <- readRDS("./Output/18S_errF.RDS")
errR <- readRDS("./Output/18S_errR.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./Output/figs/18S_error_model.png",dpi=300,height = 6,width = 6)
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
saveRDS(seqtab.nochim,"./Output/18S_seqtab.nochim.RDS")
seqtab.nochim <- readRDS("./Output/18S_seqtab.nochim.RDS")

sum(seqtab.nochim)/sum(seqtab)
dim(seqtab.nochim)
# reassign "out" to remove any missing reads
out <- readRDS("./Output/18S_trackreads.RDS")
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

write.csv(track, file = "./Output/18S_read_counts_at_each_step.csv", row.names = TRUE)

# IMPORT METADATA ####

# import and clean
meta <- read_csv("./Data/Sample_Metadata_Full.csv") %>% 
  janitor::clean_names() %>% 
  dplyr::filter(amplicon == "18S") # just bacteria samples for this
meta$sample_name %>% duplicated
# subset to match seq table sample names
meta <- meta[meta$sample_name %in% (sample.names %>% str_split("_") %>% map_chr(1)), ]
duplicated(meta$sample_name)
row.names(seqtab.nochim) <- row.names(seqtab.nochim) %>% str_split("_") %>% map_chr(1)
row.names(meta) <- meta$sample_name

negative_ctls <- grep("AMF-N",meta$sample_name,value = TRUE)
positive_ctls <- grep("AMF-P",meta$sample_name,value = TRUE)

contams <- decontam::isContaminant(seqtab = seqtab.nochim,
                                   method = 'prevalence',
                                   neg = row.names(seqtab.nochim) %in% negative_ctls,
                                   normalize = TRUE)
table(contams$contaminant)
# remove contaminants
seqtab.nochim <- seqtab.nochim[,(which(contams$contaminant != TRUE))]
ncol(seqtab.nochim)

# remove negative and positive controls
# separate positive controls for later use
seqtab.posctl <- seqtab.nochim[row.names(seqtab.nochim) %in% positive_ctls,]
seqtab.nochim <- seqtab.nochim[!row.names(seqtab.nochim) %in% negative_ctls,]

# reorder metadata to match seqtab
df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                 sample_name=row.names(seqtab.nochim))
df2 <- left_join(meta,df,by="sample_name")

duplicated(df2$sample_name)
row.names(df2) <- df2$sample_name

duplicated(df2$sample_name)
row.names(meta) <- meta$sample_name
duplicated(row.names(seqtab.nochim))
meta <- meta[row.names(seqtab.nochim),]


row.names(meta) <- meta$sample_name
meta$sample_name[duplicated(meta$sample_name)]
duplicated(row.names(seqtab.nochim))

identical(row.names(meta),row.names(seqtab.nochim))

# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./Output/18S_seqtab.nochim.clean.RDS")

# ASSIGN TAXONOMY ####

# prepare PR2 + Marjaam combined database
maarjam <- readFasta("./Taxonomy/maarjam_database_SSU.fasta.gz")
pr2 <- readFasta("./Taxonomy/pr2_version_5.0.0_SSU_dada2.fasta.gz")

# fix maarjam headers to match pr2 format
maarjam_headers_pt1 <- 
maarjam@id %>% 
  as.character() %>% 
  sub(pattern = ".*?\\b\\| \\b", replacement = "",x=.) %>% 
  str_replace(" ",";") %>% 
  str_replace(" ",";")
maarjam_family <- 
  maarjam_headers_pt1 %>% 
  str_split(";") %>% 
  map_chr(1)

maarjam_genera <- 
maarjam_headers_pt1 %>% 
  str_split(";") %>% 
  map_chr(2)
maarjam_species <- 
  maarjam_headers_pt1 %>% 
  str_split(";") %>% 
  map_chr(3)
maarjam_species <- paste0(maarjam_genera," ",maarjam_species)

maarjam_headers <- 
paste0("Eukaryota;Obazoa;Opisthokonta;Fungi;Mucoromycota;Glomeromycotina;",
       maarjam_family,
       ";",
       maarjam_genera,
       ";",
       maarjam_species)
maarjam@id <- BStringSet(maarjam_headers)
longest_read <- pr2@sread@ranges %>% as.data.frame() %>% pluck('width') %>% max()
# write new fasta file that combines pr2 with maarjam databases
writeFasta(object = pr2,
           file = "./Taxonomy/combined_pr2-maarjam_reference.fasta",
           mode = 'w',
           width=longest_read)
writeFasta(object = maarjam,
           file = "./Taxonomy/combined_pr2-maarjam_reference.fasta",
           mode = 'a',
           width=longest_read)
beepr::beep(sound=8)

# check out taxaonomic assignments of just taxa found in positive controls...
# using both databases
seqtab.posctl <- seqtab.posctl[,which(colSums(seqtab.posctl) > 0)]

# Use RDP training set for 18S
taxa <- assignTaxonomy(seqtab.nochim, 
                       "./Taxonomy/combined_pr2-maarjam_reference.fasta", 
                       multithread=nthreads,
                       taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
                       tryRC = TRUE,
                       verbose = TRUE)
pos.taxa <- assignTaxonomy(seqtab.posctl, 
                           "./Taxonomy/combined_pr2-maarjam_reference.fasta", 
                           multithread=nthreads,
                           taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
                           tryRC = TRUE,
                           verbose = TRUE)
    
# Save intermediate taxonomy file
saveRDS(taxa, file = "./Output/18S_RDP_Taxonomy_from_dada2.RDS")
saveRDS(pos.taxa, file = "./Output/18S_postive_ctl_RDP_Taxonomy_from_dada2.RDS")

beepr::beep(sound=8)

# Use Eukaryome training set for 18S (reformatted for RDP Classifier to remove leading accession and replace | with ; )
# zcat Eukaryome_General_SSU_v1.8.fasta.gz | sed 's/>[^|]*|/>/' | sed 's/|/;/g' > Eukaryome_General_SSU_v1.8.reformat.fasta; gzip Eukaryome_General_SSU_v1.8.reformat.fasta
taxa2 <- assignTaxonomy(seqtab.nochim, 
                       "./Taxonomy/Eukaryome_General_SSU_v1.8.reformat.fasta.gz", 
                       multithread=nthreads,
                       taxLevels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                       tryRC = TRUE,
                       verbose = TRUE)

pos.taxa2 <- assignTaxonomy(seqtab.posctl, 
                        "./Taxonomy/Eukaryome_General_SSU_v1.8.reformat.fasta.gz", 
                        multithread=nthreads,
                        taxLevels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                        tryRC = TRUE,
                        verbose = TRUE)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./Output/18S_RDP_Taxonomy_from_dada2_Eukaryome.RDS")
saveRDS(pos.taxa2, file = "./Output/18S_RDP_positive_ctl_Taxonomy_from_dada2_Eukaryome.RDS")
beepr::beep(sound=8)
# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print2 <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print2) <- NULL
head(taxa.print2)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa2)
met <- sample_data(meta)
sample_names(met) <- meta$sample_name

pos.meta <- meta[grepl("^AMF-P",meta$sample_name),]
pos.met <- sample_data(pos.meta)
sample_names(pos.met) <- pos.meta$sample_name
pos.tax.marjaam <- tax_table(pos.taxa)
pos.tax.eukaryome <- tax_table(pos.taxa2)
pos.otu <- otu_table(seqtab.posctl,taxa_are_rows = FALSE)

pos.ps.marjaam <- phyloseq(pos.otu,pos.tax.marjaam,pos.met)
pos.ps.eukaryome <- phyloseq(pos.otu,pos.tax.eukaryome,pos.met)
saveRDS(pos.ps.marjaam,"./Output/18S_pos_ctl_marjaam_physeq.RDS")
saveRDS(pos.ps.eukaryome,"./Output/18S_pos_ctl_eukaryome_physeq.RDS")


ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./Output/18S_ps_not-cleaned.RDS")
saveRDS(ps,"./Output/18S_ps_not-cleaned_eukaryome.RDS")

# build marjaam version of 18S physeq
tax <- tax_table(taxa)
ps2 <- phyloseq(otu,met,tax)
saveRDS(ps2,"./Output/18S_ps_not-cleaned_marjaam.RDS")





