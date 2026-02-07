library(randomForest)
library(tidyverse)
library(ranger)
library(phyloseq)
library(corncob)
library(vip)
library(openxlsx)



# import data ---------------------------------------------------------------
bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")
fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS")

# function  --------------------------------------------------------------



run_rf_growth <- function(ps,
                                     growth_predictors = c("shoot_dm", "leaf_number","final_root_dm"),
                                     gene,
                                     tax_rank          = "Genus",
                                     num_trees         = 999,
                                     seed              = 2,
                                     num_features_vip = 20,
                                     out_dir           = NULL) {
  
  # ---- packages ----
  pkgs <- c("phyloseq","janitor","corncob","ranger","tidyverse")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))
  
  # ---- checks ----
  if (!inherits(ps, "phyloseq")) stop("ps must be a phyloseq object.")
  meta0 <- as(phyloseq::sample_data(ps), "data.frame")
  
  missing_cols <- setdiff(growth_predictors, colnames(meta0))
  if (length(missing_cols)) stop("Missing columns in sample_data(ps): ", paste(missing_cols, collapse = ", "))
  
  # ---- tax glom (optional) ----
  ps_rank <- ps
  if (!is.null(tax_rank) && nzchar(tax_rank)) {
    ps_rank <- phyloseq::tax_glom(ps, taxrank = tax_rank)
  }
  
  # ---- RA table (samples x taxa) ----
  ps_ra <- phyloseq::transform_sample_counts(ps_rank, function(x) x / sum(x))
  otu_mat <- as(phyloseq::otu_table(ps_ra), "matrix")
  ra_table <- as.data.frame(otu_mat, check.names = FALSE)
  
  # ---- rename taxa columns to taxonomy labels ----
  taxa_print <- corncob::otu_to_taxonomy(phyloseq::taxa_names(ps_ra), ps_ra)
  colnames(ra_table) <- janitor::make_clean_names(taxa_print)
  
  # drop all-zero taxa columns
  ra_table <- ra_table[, colSums(ra_table) > 0, drop = FALSE]
  
  # ---- metadata (plain data.frame) ----
  meta <- data.frame(phyloseq::sample_data(ps_ra), check.names = FALSE, stringsAsFactors = FALSE)
  
  # ---- modeling base df: growth + taxa ----
  df <- meta %>%
    dplyr::select(dplyr::all_of(growth_predictors)) %>%
    dplyr::bind_cols(ra_table)
  
  # ---- fit RF per growth var ----
  top_taxa_dfs <- list()
  models <- list()
  vip_plots <- list()
  
  for (growth_var in growth_predictors) {
    
    df_1 <- df %>%
      dplyr::mutate(
        "{growth_var}" := as.numeric(as.character(.data[[growth_var]]))
      ) %>%
      dplyr::filter(!is.na(.data[[growth_var]])) %>% 
      select(-dplyr::all_of(growth_predictors[growth_predictors != growth_var]))
    
    set.seed(seed)
    rf_mod <- ranger::ranger(
      formula = stats::as.formula(paste0(growth_var, " ~ .")),
      data = df_1,
      importance = "permutation",
      num.trees = num_trees
    )
    
    models[[growth_var]] <- rf_mod
    
    p_vip <- vip::vip(rf_mod, num_features = num_features_vip) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold.italic")) +
      ggplot2::labs(
        title = "Top taxa predicting plant growth (RF regression)",
        subtitle = paste0("Response: ", growth_var)
      )
    
    vip_plots[[growth_var]] <- p_vip
    
    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      ggplot2::ggsave(
        filename = file.path(out_dir, paste0(gene, "/"),paste0(growth_var, ".png")),
        plot = p_vip,
        width = 10, height = 5, units = "in", dpi = 300
      )
    }
    
    top_20 <- vip::vi_model(rf_mod) %>%
      dplyr::arrange(dplyr::desc(Importance)) %>%
      dplyr::slice_head(n = 20)
    
    top_taxa_dfs[[growth_var]] <- top_20
    print (paste0("Completed for ", growth_var))

  }
  
  # ---- combine top taxa tables ----
  top_taxa_long <- dplyr::bind_rows(top_taxa_dfs, .id = "growth_metric") %>%
    dplyr::mutate(
      taxa = stringr::str_extract(Variable, "[^_]+_[^_]+$") # keep last two underscore chunks
    )
  
  return(list(
    df = df,
    ra_table = ra_table,
    models = models,
    vip_plots = vip_plots,
    top_taxa_dfs = top_taxa_dfs,
    top_taxa_long = top_taxa_long
  ))
}
# running rf for bact -----------------------------------------------------

res_bact <- run_rf_growth(
  ps = bact,
  gene = "16S",
  growth_predictors = c("shoot_dm", "final_root_dm", "leaf_number"),
  tax_rank = "Genus",
  out_dir = "./R/shortt_additions/figures/random_forest"
)

a <- res_bact$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

write.xlsx(exported_table, "./R/shortt_additions/figures/random_forest/16S/top_taxa_rf.xlsx")


# running rf for fung -----------------------------------------------------

res_fung <- run_rf_growth(
  ps = fung,
  gene = "ITS",
  growth_predictors = c("shoot_dm", "leaf_number","final_root_dm"),
  tax_rank = "Genus",
  out_dir = "./R/shortt_additions/figures/random_forest"
)

a <- res_fung$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

write.xlsx(exported_table, "./R/shortt_additions/figures/random_forest/ITS/top_taxa_rf.xlsx")


# subset by spec ----------------------------------------------------------
fung_sb <- subset_samples(fung, species == "Snowbrush")
fung_gr <- subset_samples(fung, species == 'GrandFir')
bact_sb <- subset_samples(bact, species == "Snowbrush")
bact_gr <- subset_samples(bact, species == "GrandFir")




#ITS snowbrush
sb_fung <- run_rf_growth(
  ps = fung_sb,
  gene = "ITS_snowbrush",
  growth_predictors = c("leaf_number"),
  tax_rank = "Genus",
  out_dir = "./R/shortt_additions/figures/random_forest"
)

a <- sb_fung$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

genus_vec <- a %>%
  filter(!grepl("g__", taxa)) %>% 
  pull(taxa) %>%
  sub("^g__", "", .) %>%
  unique()
saveRDS(genus_vec, "./R/shortt_additions/figures/random_forest/ITS_snowbrush/top_genera_ITSsb.rds")

write.xlsx(exported_table, "./R/shortt_additions/figures/random_forest/ITS_snowbrush/top_taxa_ITSsb.xlsx")


#ITS grandfir 
gr_fung <- run_rf_growth(
  ps = fung_gr,
  gene = "ITS_grandfir",
  growth_predictors = c("leaf_number"),
  tax_rank = "Genus",
  out_dir = "./R/shortt_additions/figures/random_forest"
)

a <- gr_fung$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

genus_vec <- a %>%
  filter(!grepl("g__", taxa)) %>% 
  pull(taxa) %>%
  sub("^g__", "", .) %>%
  unique()
saveRDS(genus_vec, "./R/shortt_additions/figures/random_forest/ITS_grandfir/top_genera_ITSgr.rds")

write.xlsx(exported_table, "./R/shortt_additions/figures/random_forest/ITS_grandfir/top_taxa_ITSgr.xlsx")

#16S snowbrush

sb_bact <- run_rf_growth(
  ps = bact_sb,
  gene = "16S_snowbrush",
  growth_predictors = c("leaf_number"),
  tax_rank = "Genus",
  out_dir = "./R/shortt_additions/figures/random_forest"
)

a <- sb_bact$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

genus_vec <- a %>%
  filter(!grepl("g__", taxa)) %>% 
  pull(taxa) %>%
  sub("^g__", "", .) %>%
  unique()
saveRDS(genus_vec, "./R/shortt_additions/figures/random_forest/16S_snowbrush/top_genera_16sb.rds")

write.xlsx(exported_table, "./R/shortt_additions/figures/random_forest/16S_snowbrush/top_taxa_16sb.xlsx")

#16S grandfir

gr_bact <- run_rf_growth(
  ps = bact_gr,
  gene = "16S_grandfir",
  growth_predictors = c("leaf_number"),
  tax_rank = "Genus",
  out_dir = "./R/shortt_additions/figures/random_forest"
)

a <- gr_bact$top_taxa_long
exported_table <- a %>%
  select(growth_metric, taxa, Importance)

genus_vec <- a %>%
  filter(!grepl("g__", taxa)) %>% 
  pull(taxa) %>%
  sub("^g__", "", .) %>%
  unique()
saveRDS(genus_vec, "./R/shortt_additions/figures/random_forest/16S_grandfir/top_genera_16gr.rds")
write.xlsx(exported_table, "./R/shortt_additions/figures/random_forest/16S_grandfir/top_taxa_16gf.xlsx")


# pdp --------------------------------------------------------------------

library(pdp)

ordered_tax <- sb_bact$top_taxa_long %>%  
  dplyr::arrange(dplyr::desc(Importance))

tax <- ordered_tax$Variable[4]
pdp::partial(
  object = sb_bact$models$leaf_number,
  pred.var = tax,
  train = sb_bact$df,
  plot = TRUE
)

a <- as.data.frame(sample_data(fung))
                   
#iml
#compare list to corncob after making catagorical maybe 3-5 but play around with it 
#see if there is a good cutoff for growth in terams of leaf number 
#can add burn freq into the random forest model too

#check out indic species 

#def should compare all these species to the network traits 
#individual traits have not been looked at with networks so lez do 
#git hub: statdivlab/radEmu but computationally huge 


#bact dive just sucks so don't use it 





