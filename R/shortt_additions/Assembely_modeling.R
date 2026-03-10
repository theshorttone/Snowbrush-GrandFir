library(tidyverse)
library(phyloseq)

ps <- readRDS("./Output/16S_Physeq_cleaned_w_tree.RDS")

phylo <- ps
sort(sample_sums(ps))


ps <- readRDS("./Output/16S_Physeq_cleaned_w_tree.RDS")

depth_df <- data.frame(
  SampleID = sample_names(ps),
  Depth = sample_sums(ps)
)

summary(depth_df$Depth)
quantile(depth_df$Depth, probs = c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1))


ggplot(depth_df, aes(x = Depth)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  labs(title = "Library size distribution", x = "Reads per sample", y = "Number of samples")

otu_mat <- as(otu_table(ps), "matrix")

# make sure samples are rows for vegan
if (taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

rarecurve(
  otu_mat,
  step = 1000,
  cex = 0.4,
  label = FALSE,
  xlab = "Sequencing depth",
  ylab = "Observed ASVs"
)
set.seed(666)
ps <- readRDS("./Output/16S_Physeq_cleaned_w_tree.RDS")
ps_5000 <- prune_samples(sample_sums(ps) >= 5190, ps)
ps_5000_rare <- rarefy_even_depth(ps_5000, sample.size = 5190, rngseed = 1, replace = FALSE)
#####written September 5, 2025 as a null model tutorial####
#loosely based on: 
#Knelman, Joseph E., et al. "Multiple, compounding disturbances in a forest ecosystem: fire increases susceptibility of soil edaphic properties, bacterial community structure, and function to change with extreme precipitation event." Soil Systems 3.2 (2019): 40.

#NOTE: this code uses rarefaction and OTUs, but any normalization and species-grouping will work. The tree should be trimmed to the OTUs or ASVs that are in the final table.

#####set up####
rm(list=ls());graphics.off()#clear all
#change line below to directory where your files are
#load libraries
library(permute);library(gee);library(vegan);library(ape);library(picante);library(ecodist);library(picante);
#set rarefaction depth and dataset name
rare.depth = 5190 #this parameter should be adjusted to balance keeping samples and them having enough reads; lowest recommended is 5-10k, but it could be much higher (~50)

data.set.name = 'ps_5000_rare' #this is just a name that will be used in file names for output; it should match the name of the OTU table and tree files that are read in below

#####written September 5, 2025 as a null model tutorial####
#reps are set to 100 for testing; this should be much higher, recommended 999
no.reps = 100#number of reps for null models

#load functions, note these can be altered if needed
source("R/shortt_additions/assembly_model/raup_crick_abundance_one_comparison.R")
source("R/shortt_additions/assembly_model/beta.nti.R")

#####Optional Pre-processing####
set.seed(666)

ps <- readRDS("./Output/16S_Physeq_cleaned_w_tree.RDS")

ps_5000 <- prune_samples(sample_sums(ps) >= 5190, ps)

ps_5000_rare <- rarefy_even_depth(
  ps_5000,
  sample.size = 5190,
  rngseed = 1,
  replace = FALSE
)

#### Extract OTU table, trimmed tree, set number of samples ####
otu <- as(otu_table(ps_5000_rare), "matrix")

if (taxa_are_rows(ps_5000_rare)) {
  otu <- t(otu)
}

phylo <- phy_tree(ps_5000_rare)

phylo$tip.label <- gsub("'", "", phylo$tip.label)  # probably unnecessary, but OK if old code needed it

# reorder OTU columns to match tree tip labels
otu <- otu[, phylo$tip.label]

match.phylo.otu <- match.phylo.data(phylo, as.data.frame(t(otu)))

no.of.samples <- nrow(otu)

####### Step 1: calculate observed betaMNTD #####
#this can be computationally intensive
date(); beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T)); date();
#dim(beta.mntd.weighted);beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,paste(data.set.name,"_bMNTD_weighted.csv",sep=""),quote=F);

##### Step 2: calculate randomized betaMNTD #####
rep = 1:no.reps

for (i in rep){
  print(i)
  rand.weighted.beta.mntd = as.matrix(comdistnt(comm = t(match.phylo.otu$data), 
                                                dis = taxaShuffle(cophenetic(match.phylo.otu$phy)), 
                                                abundance.weighted = T)); ## randomized beta.mntd
  write.csv(rand.weighted.beta.mntd,paste(data.set.name ,"_bMNTD_weighted_rep_",rep[i],".csv",sep=""),quote=F);
  print(date());
  rm(rand.weighted.beta.mntd)
}

##### Step 3: calculate betaNTI (must have observed beta.mntd and all null model reps run) #####
beta.mntd.weighted = read.csv(paste(data.set.name,"_bMNTD_weighted.csv",sep=""),row.names=1);

#get all null model files
all.files = sort(as.numeric(gsub(".*_rep_|\\.csv", "", list.files("./", pattern = paste0(data.set.name, "_bMNTD_weighted_rep_\\d+\\.csv$")))))

## find beta.mntd reps that still need to be done; if there are any, need to do them before proceeding
reps.to.do = c(1:no.reps)[-which(1:no.reps %in% all.files)];reps.to.do#should be empty
#write.table(t(reps.to.do),"./reps.to.do.txt",sep=" ",quote=F,row.names=F,col.names=F)

beta.nti.weighted = beta.nti.calc.stegen(samp = otu, reps = all.files, 
                                         path.to.reps = paste('./',data.set.name,'_bMNTD_weighted_rep_',sep=""),
                                         beta.mntd.obs = beta.mntd.weighted);
write.csv(beta.nti.weighted,paste(data.set.name,'_bNTI_weighted.csv',sep=""),quote=F);

# optional plot to check distribution
# pdf(paste(data.set.name,'_bNTI_weighted.pdf',sep=""));
# hist(as.dist(beta.nti.weighted),xlim = c(min(c(-2,min(as.dist(beta.nti.weighted)))),max(c(2,max(as.dist(beta.nti.weighted))))),xlab="betaNTI - Weighted",cex.lab=1.3,main=""); abline(v=c(-2,2),col=2,lwd=2);
# dev.off();

##### Step 4: calculate Bray-Curtis ####
bray = as.matrix(distance(t(match.phylo.otu$data),method = 'bray-curtis'))
write.csv(bray,paste(data.set.name,"_bray_weighted.csv",sep=""),quote=F);

# optional quick check, should be a non-linearily increasing relationship
# identical(colnames(beta.mntd.weighted),colnames(bray.out))
# plot(as.dist(beta.mntd.weighted) ~ as.dist(bray.out))

#### Step 5 Raup-Crick Bray####
metric = 'RC'
abund = T #option to set to false and not have abundance-weighted RC (uncommon)
abund.for.names = 'weighted' # 'unweighted' (uncommon)

#generate null model output
for (i in 1:(no.of.samples - 1)){
  for (j in (i + 1):(no.of.samples)){
    raup.crick.out = raup_crick_abundance_one_comparison(null.one.use = i,null.two.use = j,otu, plot_names_in_col1=FALSE, classic_metric=FALSE, split_ties=TRUE, reps=999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE);
    print(raup.crick.out)
    print(c(i,j))
    write.csv(raup.crick.out,paste('RC_',abund.for.names,'_Diss_comm1_',i,'_comm2_',j,'.csv',sep=""),row.names=T,quote=F);
  }
}

all.RC <- grep("RC_weighted_Diss_comm1_", list.files(pattern = "RC_"), value = TRUE)

#create matrix of comparisons
raup.crick.out = matrix(NA, nrow = nrow(bray), ncol = ncol(bray), dimnames = dimnames(bray))

for (i in 1:length(all.RC)) {
  RC.temp = read.csv(paste(all.RC[i],sep=""),#reading in the file is more complicated because "T" and "F" were being converted to "TRUE" and "FALSE"; this fixes that
                     header = TRUE,check.names = FALSE,
                     stringsAsFactors = FALSE,colClasses = "character")
  row.names(RC.temp) <- RC.temp[[1]];RC.temp <- RC.temp[, -1, drop = FALSE]
  if (nrow(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  if (ncol(RC.temp) != 1 ) {print(c(i,"ERROR"))} 
  raup.crick.out[which(rownames(raup.crick.out) == rownames(RC.temp)),which(colnames(raup.crick.out) == colnames(RC.temp))] = RC.temp[1,1]
  print(RC.temp)
}

raup.crick.out = as.data.frame(as.matrix(as.dist(raup.crick.out)))
write.csv(raup.crick.out,paste(data.set.name,'_RC_weighted.csv', sep = ''),quote=F);