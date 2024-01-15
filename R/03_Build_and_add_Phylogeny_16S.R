# -----------------------------------------------------------------------------#
# Syringodium isoetifolium 16S phylogeny
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
#                     DECIPHER 2.24.0
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
# library(msa); packageVersion("msa")
library(ape); packageVersion("ape")
library(seqinr); packageVersion("seqinr")
library(DECIPHER); packageVersion("DECIPHER")

# Read in phyloseq object from first script output ####
ps <- readRDS("./Output/16S_ps_not-cleaned.RDS")


# simplify ASV names
seqs <- rownames(tax_table(ps))
names(seqs) <- paste0("ASV_",1:length(seqs)) # This propagates to the tip labels of the tree

# create DNAStringSet object
seqs_StringSet <- DNAStringSet(seqs)

# Multiple sequence alignment  ####
decipher_alignment <- DECIPHER::AlignSeqs(seqs_StringSet, processors = parallel::detectCores() - 1,verbose = TRUE) # DECIPHER method
adj_decipher_alignment <- DECIPHER::AdjustAlignment(decipher_alignment,processors = parallel::detectCores() - 1)

saveRDS(adj_decipher_alignment,"./Output/16S_dna_alignment_DECIPHER.RDS")


# Convert to various formats
phang.align <- as.phyDat(as.character(adj_decipher_alignment), type = "DNA")
dnabin.align <- as.DNAbin(adj_decipher_alignment)


# distance - maximum likelihood ####
dm <- dist.ml(phang.align,model = "JC69")


# dm.TN93 <- dist.dna(dnabin.align,model = "TN93")
# missing <- which(is.na(dm.TN93))
# possible_positions <- round(missing / length(labels(dm.TN93)),0)
# 
# labels(dm.TN93)[possible_positions]



#save
saveRDS(dm,"./Output/16S_ML_Distance.RDS")

# Initial neighbor-joining tree ####
treeNJ <- njs(dm) # Note, tip order != sequence order

# save progress
saveRDS(treeNJ, "./Output/16S_treeNJ.RDS")

# Estimate model parameters ####
fit <-  pml(treeNJ, data=as.phyDat(muscle_alignment))

#save
saveRDS(fit,"./Output/16S_fit_treeNJ_muscle.RDS")


fit$tree$tip.label <- seqs

# add tree to phyloseq object ####
ps2 <- phyloseq(tax_table(tax_table(ps)),
                otu_table(otu_table(ps)),
                sample_data(sample_data(ps)),
                phy_tree(fit$tree))


# Save updated phyloseq object with tree
saveRDS(ps2, "./Output/ps_not-cleaned_w_tree.RDS")


