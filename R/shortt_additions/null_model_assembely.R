#!/usr/bin/env Rscript

rm(list = ls())
graphics.off()

.libPaths(c("/home/zachary.shortt/R/lib",
            "/data/lab/cheeke/R_libs",
            .libPaths()))

library(picante)

args <- commandArgs(trailingOnly = TRUE)
start_rep <- as.integer(args[1])
end_rep   <- as.integer(args[2])

data.set.name <- "snowbrush_april_2026"
out_dir  <- "./outputs/assembly_model"
null_dir <- file.path(out_dir, "null_reps")

dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)

comm <- readRDS(file.path(out_dir, paste0(data.set.name, "_comm.rds")))
phylo_dist <- readRDS(file.path(out_dir, paste0(data.set.name, "_phylo_dist.rds")))

for (i in start_rep:end_rep) {
  
  outfile <- file.path(
    null_dir,
    sprintf("%s_bMNTD_weighted_rep_%03d.rds", data.set.name, i)
  )
  
  if (file.exists(outfile)) next
  
  set.seed(100000 + i)
  
  rand.weighted.beta.mntd <- as.matrix(
    comdistnt(
      comm = comm,
      dis = taxaShuffle(phylo_dist),
      abundance.weighted = TRUE
    )
  )
  
  saveRDS(rand.weighted.beta.mntd, outfile)
  rm(rand.weighted.beta.mntd)
  gc()
}