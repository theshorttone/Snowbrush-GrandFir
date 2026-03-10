library(tidyverse)
library(vegan)
library(phyloseq)

fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS") %>% 
  subset_samples(inoculum_site != "Sterile")

bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS") %>% 
  subset_samples(inoculum_site != "Sterile")


am <- readRDS("./Output/phyloseq_objects/18S_clean_phyloseq_object.RDS") %>% 
  subset_samples(inoculum_site != "Sterile")



a <- data.frame(sample_data(fung))

# permanova function:  ----------------------------------------------------
run_permanova_phyloseq <- function(
    ps,
    vars = c("inoculum_site", "species", "drought"),
    permutations = 9999,
    seed = 1,
    transform = c("relabund", "hellinger", "none"),
    distance = "bray"
) {
  transform <- match.arg(transform)
  
  if (!inherits(ps, "phyloseq")) stop("`ps` must be a phyloseq object.")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Need package: phyloseq")
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Need package: vegan")
  
  # --- community matrix (samples x taxa) ---
  Y <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) Y <- t(Y)
  
  # --- metadata ---
  meta <- data.frame(phyloseq::sample_data(ps), check.names = FALSE)
  
  # ensure vars exist
  missing_vars <- setdiff(vars, colnames(meta))
  if (length(missing_vars) > 0) {
    stop("Missing variables in sample_data(ps): ", paste(missing_vars, collapse = ", "))
  }
  
  n0 <- nrow(meta)  # starting sample count
  
  # drop samples with NA in any predictor
  keep <- stats::complete.cases(meta[, vars, drop = FALSE])
  removed_na <- sum(!keep)
  if (removed_na > 0) {
    message("Removed ", removed_na, " sample(s) due to NA in: ", paste(vars, collapse = ", "))
  }
  Y <- Y[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  
  # drop samples with 0 reads (avoids NaNs in relabund)
  rs <- rowSums(Y)
  removed_zero <- sum(rs == 0)
  if (removed_zero > 0) {
    message("Removed ", removed_zero, " sample(s) with 0 total reads after filtering.")
    Y <- Y[rs > 0, , drop = FALSE]
    meta <- meta[rs > 0, , drop = FALSE]
  }
  
  # final removal summary
  n1 <- nrow(meta)
  if (n1 == n0) {
    message("No samples removed (n = ", n1, ").")
  } else {
    message("Total samples kept: ", n1, " of ", n0, ".")
  }
  
  # drop taxa that are now all-zero (cleanup)
  Y <- Y[, colSums(Y) > 0, drop = FALSE]
  
  # coerce predictors to factors
  for (v in vars) meta[[v]] <- as.factor(meta[[v]])
  
  # --- transform ---
  Y_tr <- Y
  if (transform == "relabund") {
    Y_tr <- sweep(Y, 1, rowSums(Y), "/")
  } else if (transform == "hellinger") {
    Y_rel <- sweep(Y, 1, rowSums(Y), "/")
    Y_tr <- sqrt(Y_rel)
  }
  
  # --- distance ---
  d <- vegan::vegdist(Y_tr, method = distance)
  
  # --- PERMANOVA (sequential / order matters) ---
  set.seed(seed)
  form <- stats::as.formula(paste("d ~", paste(vars, collapse = " + ")))
  ad <- vegan::adonis2(form, data = meta, permutations = permutations, by = "terms")
  
  # --- betadisper for each predictor ---
  bd <- lapply(vars, function(v) {
    bdv <- vegan::betadisper(d, meta[[v]])
    list(
      betadisper = bdv,
      anova = stats::anova(bdv),
      permutest = vegan::permutest(bdv)
    )
  })
  names(bd) <- vars
  
  list(
    adonis2 = ad,
    betadisper = bd,
    kept_samples = rownames(meta),
    meta = meta
  )
}


# func 2 ------------------------------------------------------------------


res <- run_permanova_phyloseq(
  ps = fung,
  vars = c("inoculum_site", "species", "drought"),
  transform = "relabund",
  distance = "bray",
  permutations = 999,
)

res$adonis2
res$betadisper$species$permutest
res$betadisper$drought$permutest
res$betadisper$inoculum_site$permutest

res <- run_permanova_phyloseq(
  ps = bact,
  vars = c("inoculum_site", "species", "drought"),
  transform = "relabund",
  distance = "bray",
  permutations = 999,
)

res$adonis2
res$betadisper$species$permutest
res$betadisper$drought$permutest
res$betadisper$inoculum_site$permutest

res <- run_permanova_phyloseq(
  ps = am,
  vars = c("inoculum_site", "drought"),
  transform = "relabund",
  distance = "bray",
  permutations = 999,
)

res$adonis2
res$betadisper$drought$permutest
res$betadisper$inoculum_site$permutest



# -------------------------------------------------------------------------

