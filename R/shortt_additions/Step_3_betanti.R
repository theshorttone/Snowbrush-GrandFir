data.set.name <- "snowbrush_april_2026"

source("R/shortt_additions/assembly_model/beta.nti.R")

out_dir  <- "./outputs/assembly_model"
null_dir <- file.path(out_dir, "null_reps")

input_rds <- file.path(out_dir, paste0(data.set.name, "_OTU_table.RDS"))
otu <- readRDS(input_rds)

# observed bMNTD matrix
beta.mntd.weighted <- read.csv(
  file.path(out_dir, paste0(data.set.name, "_bMNTD_weighted.csv")),
  row.names = 1,
  check.names = FALSE
)

# get all null model rep numbers
all.files <- sort(as.integer(
  gsub(
    ".*_rep_|\\.csv", "",
    list.files(
      null_dir,
      pattern = paste0(data.set.name, "_bMNTD_weighted_rep_\\d+\\.csv$")
    )
  )
))

beta.nti.weighted <- beta.nti.calc.stegen(
  samp = otu,
  reps = all.files,
  path.to.reps = file.path(null_dir, paste0(data.set.name, "_bMNTD_weighted_rep_")),
  beta.mntd.obs = beta.mntd.weighted
)

write.csv(
  beta.nti.weighted,
  file.path(out_dir, paste0(data.set.name, "_bNTI_weighted.csv")),
  quote = FALSE
)