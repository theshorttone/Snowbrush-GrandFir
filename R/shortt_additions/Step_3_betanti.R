out_dir <- "./outputs/assembly_model"
null_dir <- file.path(out_dir, "null_reps")

beta.mntd.weighted <- read.csv(
  file.path(out_dir, paste0(data.set.name, "_bMNTD_weighted.csv")),
  row.names = 1
)

#get all null model files
all.files <- sort(as.numeric(
  gsub(
    ".*_rep_|\\.csv", "",
    list.files(
      null_dir,
      pattern = paste0(data.set.name, "_bMNTD_weighted_rep_\\d+\\.csv$")
    )
  )
))

## find beta.mntd reps that still need to be done; if there are any, need to do them before proceeding
reps.to.do <- c(1:no.reps)[-which(1:no.reps %in% all.files)]
reps.to.do

#write.table(t(reps.to.do),
#            file.path(out_dir, "reps.to.do.txt"),
#            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

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