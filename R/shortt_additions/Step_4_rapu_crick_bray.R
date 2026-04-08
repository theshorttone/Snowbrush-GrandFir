#!/usr/bin/env Rscript

rm(list = ls())
graphics.off()

.libPaths(c("/home/zachary.shortt/R/lib",
            "/data/lab/cheeke/R_libs",
            .libPaths()))

suppressPackageStartupMessages({
  library(ecodist)
})

set.seed(666)

source("R/shortt_additions/assembly_model/Raup_Crick_Abundance_One_Comparison.R")

# ----------------------------
# Settings
# ----------------------------
data.set.name <- "snowbrush_april_2026"
otu_dir  <- "./outputs/assembly_model"
input_rds <- file.path(otu_dir, paste0(data.set.name, "_OTU_table.RDS"))
otu <- readRDS(input_rds)
out_dir <- file.path("./outputs/assembly_model/rc_files", paste0(data.set.name, "_RC_files"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rc_reps <- 999
abund.for.names <- "weighted"

cat("Starting Raup-Crick Bray run\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Working dir:", getwd(), "\n")
cat("Input exists:", file.exists(input_rds), "\n")
cat("Output dir:", out_dir, "\n\n")

# ----------------------------
# Load OTU matrix
# ----------------------------
otu <- readRDS(input_rds)

stopifnot(nrow(otu) >= 2, ncol(otu) >= 2)

cat("OTU dimensions (samples x taxa):", nrow(otu), "x", ncol(otu), "\n")
cat("OTU object size:", format(object.size(otu), units = "GB"), "\n")

# Build Bray only for final matrix dimensions/names
bray <- as.matrix(distance(otu, method = "bray-curtis"))
bray_nrow <- nrow(bray)
bray_ncol <- ncol(bray)
bray_dimnames <- dimnames(bray)

no.of.samples <- nrow(otu)

cat("Number of samples:", no.of.samples, "\n")
cat("Number of pairwise comparisons:", choose(no.of.samples, 2), "\n\n")

rm(bray)
gc()

# ----------------------------
# Build unique pair table
# ----------------------------
pair_df <- as.data.frame(t(combn(no.of.samples, 2)))
names(pair_df) <- c("i", "j")

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
  
  rm(rc)
  gc(verbose = FALSE)
  
  paste("Finished:", i, j)
}

# ----------------------------
# Run pairwise comparisons serially
# ----------------------------
msgs <- lapply(seq_len(nrow(pair_df)), run_one_rc)

# ----------------------------
# Rebuild full RC matrix from per-pair files
# ----------------------------
all.RC <- grep(
  "^RC_weighted_Diss_comm1_.*_comm2_.*\\.csv$",
  list.files(out_dir),
  value = TRUE
)

if (length(all.RC) == 0) {
  stop("No RC files were found. Aborting.")
}

raup.crick.out <- matrix(
  NA_real_,
  nrow = bray_nrow,
  ncol = bray_ncol,
  dimnames = bray_dimnames
)

for (f in all.RC) {
  RC.temp <- read.csv(
    file.path(out_dir, f),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
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
  
  val <- as.numeric(RC.temp[1, 1])
  raup.crick.out[rr, cc] <- val
  raup.crick.out[cc, rr] <- val
}

diag(raup.crick.out) <- 0

final_out <- file.path(
  "./Output/assembly_model",
  paste0(data.set.name, "_RC_weighted.csv")
)

write.csv(raup.crick.out, final_out, quote = FALSE)

cat("Done.\n")
cat("Final matrix written to:", final_out, "\n")
cat("Finish time:", as.character(Sys.time()), "\n")