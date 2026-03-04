# BUILD PHYLOQGENIES ####

# SETUP ####
set.seed(666)
#mylib <- "/uufs/chpc.utah.edu/common/home/u6033249/R/library-4.4"

.libPaths(c("/home/zachary.shortt/R/lib",
            "/data/lab/cheeke/R_libs", .libPaths()))

## packages ####
library(tidyverse)
library(phyloseq)
library(DECIPHER)
library(janitor)
library(phylogram) # on CRAN
# library(tidyverse)
# library(phyloseq)
# library(DECIPHER)
# library(phylogram)


## functions ####

## load physeq objects ####
#ps_its <- readRDS("./Output/ITS_clean_phyloseq_object.RDS)
ps_16s <- readRDS("./output/16S_clean_phyloseq_object.RDS")

## extract sequences ####

seqs_16s <- rownames(tax_table(ps_16s))
names(seqs_16s) <- paste0("ASV_",1:length(seqs_16s)) # This propagates to the tip labels of the tree
#seqs_its <- rownames(tax_table(ps_its))
#names(seqs_its) <- paste0("ASV_",1:length(seqs_its)) # This propagates to the tip labels of the tree



# ALIGNMENT ####
alignment_16s <- DECIPHER::AlignSeqs(DNAStringSet(seqs_16s),processors=NULL)
saveRDS(alignment_16s,"./output/16S_dna_alignment.RDS")

Biostrings::writeXStringSet(alignment_16s, "./output/16S_ASV_aligned.fasta")

#Build Tree with FastTree
cmd <- system("FastTree -nt -gtr -gamma ./output/16S_ASV_aligned.fasta > ./output/16S_fasttree.nwk")
status <- system(cmd)
stopifnot(status == 0, file.exists("./output/16S_fasttree.nwk"))

tree_16s <- ape::read.tree("./output/16S_fasttree.nwk")

# relabel tips from ASV_# -> phyloseq taxa_names (your original sequence strings)
map_16s <- setNames(taxa_names(ps_16s), names(seqs_16s))
tree_16s$tip.label <- map_16s[tree_16s$tip.label]

# sanity checks (critical for βNTI)
stopifnot(length(tree_16s$tip.label) == ntaxa(ps_16s))
stopifnot(setequal(tree_16s$tip.label, taxa_names(ps_16s)))
stopifnot(!is.null(tree_16s$edge.length))

# attach tree to phyloseq
ps_16s_w_tree <- merge_phyloseq(ps_16s, phy_tree(tree_16s))

# export
saveRDS(ps_16s_w_tree, "./output/16S_Physeq_cleaned_w_tree.RDS")
#saveRDS(ps_its_w_tree,"./output/ITS_Physeq_cleaned_w_tree.RDS")