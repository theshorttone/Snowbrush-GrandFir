# -----------------------------------------------------------------------------#
# Building and adding a phylogeny to the cleaned phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.2
#                     phangorn v 2.10.0
#                     phyloseq v 1.40.0
#                     msa v 1.28.0
#                     ape v 5.6.2
#                     seqinr v 4.2.16
#                     DECIPHER v 2.24.0
#                     parallel v 4.3.2
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Perform multiple sequence alignment of all ASVs, build distance matrix,       # 
# construct and refine a phylogenetic tree, add the tree to the phyloseq object #
#           With larger data sets, this can be a long process...                #
# Further, proper phylogenetics is beyond the scope of this tutorial.           #
#################################################################################

# Packages and functions ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(phangorn); packageVersion("phangorn")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(DECIPHER); packageVersion("DECIPHER")
library(parallel); packageVersion("parallel")
# library(fastreeR); packageVersion("fastreeR")

# number of threads for parallel processing
nthreads <- parallel::detectCores() - 1

# Read in phyloseq object from first script output ####
ps <- readRDS("./Output/16S_ps_not-cleaned.RDS")


# simplify ASV names
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# create DNAStringSet object
seqs_StringSet <- DNAStringSet(seqs)
ShortRead::writeFasta(DNAStringSet(seqs),"./Output/16S_ASVs_seqs.fasta")

# Multiple sequence alignment  ####
decipher_alignment <- DECIPHER::AlignSeqs(seqs_StringSet, processors = parallel::detectCores() - 1,verbose = TRUE) # DECIPHER method
adj_decipher_alignment <- DECIPHER::AdjustAlignment(decipher_alignment,processors = parallel::detectCores() - 1)

saveRDS(adj_decipher_alignment,"./Output/16S_dna_alignment_DECIPHER.RDS")
adj_decipher_alignment <- readRDS("./Output/16S_dna_alignment_DECIPHER.RDS")
ShortRead::writeFasta(adj_decipher_alignment,"./Output/16S_ASVs_seqs_aligned.fasta")


# Convert to various formats
align_character <- as.character(decipher_alignment)

# phang.align <- as.phyDat(align_character, type = "DNA")
# dnabin.align <- as.DNAbin(adj_decipher_alignment)

# distance - maximum likelihood ####
dm <- DECIPHER::DistanceMatrix(adj_decipher_alignment)

#save
saveRDS(dm,"./Output/16S_ML_Distance.RDS")
dm <- readRDS("./Output/16S_ML_Distance.RDS")
dist_dm <- as.dist(dm)



# Initial maximum likelihood tree ####
# treeNJ <- fastreeR::dist2tree(inputDist = dist_dm)


# Build tree outside of R
# RAxML-NG ML tree
system2("raxml-ng",args = c("-msa ./Output/16S_ASVs_seqs_aligned.fasta",paste("--workers",nthreads),
                            "--seed 666","--all","--model GTR+G","--bs-trees autoMRE{1000}"))

# import tree
x <- ape::read.tree(file = "./Output/16S_raxml.tre")

# replace tip labels with ASV sequences for import into phyloseq
x$tip.label <- seqs[x$tip.label]

#save
saveRDS(x,"./Output/16S_RAxML_tree_w_labels.RDS")

# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(x))

# sanity check that tip labels are in right order
taxa_names(ps2)[1] == ps2@phy_tree$tip.label[1]

plot_tree(ps2,color = "Phylum")
ggsave("./Output/16S_Tree_Plot.png",width = 36,height = 8)
# Save updated phyloseq object with tree
saveRDS(ps2, "./Output/16S_ps_not-cleaned_w_tree.RDS")


