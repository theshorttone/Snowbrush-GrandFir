#!/usr/bin/env Rscript

rm(list = ls())
graphics.off()

.libPaths(c("/home/zachary.shortt/R/lib",
            "/data/lab/cheeke/R_libs",
            .libPaths()))

suppressPackageStartupMessages({
  library(phyloseq)
  library(ecodist)
  library(parallel)
})

set.seed(666)

source("./R/shortt_additions/assembly_model/raup_crick_abundance_one_comparison.R")

# ----------------------------
# Settings
# ----------------------------
data.set.name <- "ps_5000_rare"
input_rds <- "./Output/ps_5000_rare.RDS"

# Put RC outputs in their own folder
out_dir <- file.path("./Output/assembly_model", paste0(data.set.name, "_RC_files"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Number of RC null reps per comparison
rc_reps <- 999

# Number of cores from SLURM; fallback to 1 if running interactively
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(ncores) || ncores < 1) ncores <- 1

cat("Starting Raup-Crick Bray run\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Cores:", ncores, "\n")
cat("Output dir:", out_dir, "\n\n")

# ----------------------------
# Load phyloseq object and rebuild OTU matrix
# ----------------------------
ps_5000_rare <- readRDS(input_rds)

otu <- as(otu_table(ps_5000_rare), "matrix")
if (taxa_are_rows(ps_5000_rare)) {
  otu <- t(otu)
}

# sanity checks
stopifnot(nrow(otu) >= 2, ncol(otu) >= 2)
cat("OTU dimensions (samples x taxa):", nrow(otu), "x", ncol(otu), "\n")

# rebuild bray only to get final matrix dimensions/names
bray <- as.matrix(distance(otu, method = "bray-curtis"))
no.of.samples <- nrow(otu)
abund.for.names <- "weighted"

cat("Number of samples:", no.of.samples, "\n")
cat("Number of pairwise comparisons:", choose(no.of.samples, 2), "\n\n")

# ----------------------------
# Build pair table
# ----------------------------
pair_df <- expand.grid(
  i = seq_len(no.of.samples - 1),
  j = seq_len(no.of.samples),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

pair_df <- pair_df[pair_df$i < pair_df$j, , drop = FALSE]
row.names(pair_df) <- NULL

# ----------------------------
# One comparison
# ----------------------------
run_one_rc <- function(k) {
  i <- pair_df$i[k]
  j <- pair_df$j[k]
  
  out_file <- file.path(
    out_dir,
    paste0("RC_", abund.for.names, "_Diss_comm1_", i, "_comm2_", j, ".csv")
  )
  
  if (file.exists(out_file)) {
    return(paste("Skipping existing:", i, j))
  }
  
  rc <- raup_crick_abundance_one_comparison(
    null.one.use = i,
    null.two.use = j,
    otu,
    plot_names_in_col1 = FALSE,
    classic_metric = FALSE,
    split_ties = TRUE,
    reps = rc_reps,
    set_all_species_equal = FALSE,
    as.distance.matrix = TRUE,
    report_similarity = FALSE
  )
  
  write.csv(rc, out_file, row.names = TRUE, quote = FALSE)
  
  paste("Finished:", i, j)
}

# ----------------------------
# Run pairwise comparisons in parallel
# ----------------------------
cat("Launching parallel RC comparisons...\n")
cat("Time:", as.character(Sys.time()), "\n\n")

msgs <- mclapply(
  X = seq_len(nrow(pair_df)),
  FUN = run_one_rc,
  mc.cores = ncores,
  mc.preschedule = FALSE
)

# print a subset of messages to the log
cat("Example worker messages:\n")
print(utils::head(unlist(msgs), 20))
cat("\n")

# ----------------------------
# Rebuild full RC matrix from per-pair files
# ----------------------------
cat("Rebuilding final RC matrix...\n")
cat("Time:", as.character(Sys.time()), "\n\n")

all.RC <- grep(
  "^RC_weighted_Diss_comm1_.*_comm2_.*\\.csv$",
  list.files(out_dir),
  value = TRUE
)

cat("RC files found:", length(all.RC), "\n")
cat("Expected files:", choose(no.of.samples, 2), "\n\n")

if (length(all.RC) == 0) {
  stop("No RC files were found. Aborting.")
}

raup.crick.out <- matrix(
  NA,
  nrow = nrow(bray),
  ncol = ncol(bray),
  dimnames = dimnames(bray)
)

for (f in all.RC) {
  RC.temp <- read.csv(
    file.path(out_dir, f),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  
  row.names(RC.temp) <- RC.temp[[1]]
  RC.temp <- RC.temp[, -1, drop = FALSE]
  
  if (nrow(RC.temp) != 1 || ncol(RC.temp) != 1) {
    warning("Unexpected dimensions in file: ", f)
    next
  }
  
  rr <- which(rownames(raup.crick.out) == rownames(RC.temp))
  cc <- which(colnames(raup.crick.out) == colnames(RC.temp))
  
  if (length(rr) != 1 || length(cc) != 1) {
    warning("Could not match row/col names for file: ", f)
    next
  }
  
  raup.crick.out[rr, cc] <- RC.temp[1, 1]
}

# convert upper triangle vector into full dist-style object then matrix/data.frame
raup.crick.out <- as.data.frame(as.matrix(as.dist(raup.crick.out)))

final_out <- file.path(
  "./Output/assembly_model",
  paste0(data.set.name, "_RC_weighted.csv")
)

write.csv(raup.crick.out, final_out, quote = FALSE)

cat("Done.\n")
cat("Final matrix written to:", final_out, "\n")
cat("Finish time:", as.character(Sys.time()), "\n")