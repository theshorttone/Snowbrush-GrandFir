
library(corncob)
library(phyloseq)
library(tidyverse)
fung <- readRDS("./Output/phyloseq_objects/ITS_clean_phyloseq_object.RDS")
fung <- fung %>% tax_glom("Genus")


# functions ---------------------------------------------------------------

corncob_sig_to_rf_genus <- function(fit, ps) {
  # taxon IDs corncob returns (in your case: sequences)
  sig_taxa <- fit$significant_taxa
  if (is.null(sig_taxa)) sig_taxa <- names(fit$significant_models)
  
  # map taxon_id -> Genus
  genus_raw <- as.data.frame(tax_table(ps)) %>%
    rownames_to_column("taxon_id") %>%
    filter(taxon_id %in% sig_taxa) %>%
    pull(Genus)
  
  genus_rf <- genus_raw %>%
    tolower() %>%
    str_replace("^g__", "g_") %>%
    str_replace_all("[^a-z0-9]+", "_") %>%   # mimic janitor-ish style
    str_replace("_+$", "") %>%
    na.omit() %>%
    unique()
  
  genus_rf
}

# fungi -------------------------------------------------------------------


sam_df <- data.frame(sample_data(fung), check.names = FALSE)

sam_abies <- sam_df %>% filter(host == "Abies grandis") %>% 
  select("leaf_number") %>% 
  mutate(
    leaf_number = as.numeric(leaf_number)
  )

ggplot(sam_abies, aes(x = leaf_number))+
  geom_bar()

x <- sam_abies$leaf_number
mu <- mean(x, na.rm = TRUE)
sdx <- sd(x, na.rm = TRUE)
lo  <- mu - sdx
hi  <- mu + sdx

sam_cean <- sam_df %>% filter(host == "Ceanothus velutinus") %>% 
  select("leaf_number") %>% 
  mutate(
    leaf_number = as.numeric(leaf_number)
  )

ggplot(sam_cean, aes(x = leaf_number))+
  geom_bar()

x1 <- sam_cean$leaf_number
mu1 <- mean(x1, na.rm = TRUE)
sdx1 <- sd(x1, na.rm = TRUE)
lo1  <- mu1 - sdx1
hi1  <- mu1 + sdx1

# 1) create category + counts
sam_cat_df <- sam_df %>%
  mutate(
    grand_fir_growth = case_when(
      is.na(leaf_number)  ~ NA_character_,
      leaf_number < lo       ~ "low",
      leaf_number > hi       ~ "high",
      TRUE                ~ "medium"
    ),
    grand_fir_growth = factor(
      grand_fir_growth,
      levels = c("low", "medium", "high")
    ),
    ceanothus_growth = case_when(
      is.na(leaf_number)  ~ NA_character_,
      leaf_number < lo1      ~ "low",
      leaf_number > hi1      ~ "high",
      TRUE                ~ "medium"
    ),
    ceanothus_growth = factor(
      ceanothus_growth,
      levels = c("low", "medium", "high")
    )
  )
  

 sample_data(fung) <- sample_data(sam_cat_df)

ps_filt <- fung %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 10, .)

ps <- ps_filt

sample_data(ps)$inoculum_site <- factor(sample_data(ps)$inoculum_site)
sample_data(ps)$fire_freq <- factor(sample_data(ps)$fire_freq)

ps_ceanotus <- ps %>% subset_samples(host == "Ceanothus velutinus", na.rm = TRUE) %>% 
  subset_samples(!(ceanothus_growth == "medium")) %>%
  subset_samples(!is.na(ceanothus_growth))

ps_grandfir <- ps %>% subset_samples(host == "Abies grandis", na.rm = TRUE) %>% 
  subset_samples(!(grand_fir_growth == "medium")) %>%
  subset_samples(!is.na(grand_fir_growth))
  
a <- data.frame(sample_data(ps_ceanotus))

#change to the random forest model 
#use vIP package to interpret bc forest is confusing 
fitgr <- differentialTest(
  formula = ~ grand_fir_growth,
  phi.formula = ~ 1,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  test = "Wald",
  fdr_cutoff = 0.05,
  data = ps_gr_genus
)

fitcn <- differentialTest(
  formula = ~ ceanothus_growth,
  phi.formula = ~ 1,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  test = "Wald",
  fdr_cutoff = 0.05,
  data = ps_cn_genus
)

plot(fitgr)
plot(fitcn)

rf_ITS_cn <- readRDS("./R/shortt_additions/figures/random_forest/ITS_snowbrush/top_genera_ITSsb.rds")
rf_ITS_gr <- readRDS("./R/shortt_additions/figures/random_forest/ITS_grandfir/top_genera_ITSgr.rds")

corncob_gr_genus <- corncob_sig_to_rf_genus(fitgr, ps_grandfir)
corncob_cn_genus <- corncob_sig_to_rf_genus(fitcn, ps_ceanotus)

# now compare directly to your RF genus vectors
overlap_gr <- intersect(corncob_gr_genus, rf_ITS_gr)
overlap_cn <- intersect(corncob_cn_genus, rf_ITS_cn)

overlap_gr
overlap_cn

# bacteria ----------------------------------------------------------------

bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")
bact <- bact %>% tax_glom("Genus")
ASV_names <- otu_table(bact) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,bact,level = c("Phylum","Class","Order","Family","Genus"))



sam_df <- data.frame(sample_data(bact), check.names = FALSE)

sam_abies <- sam_df %>% filter(host == "Abies grandis") %>% 
  select("leaf_number") %>% 
  mutate(
    leaf_number = as.numeric(leaf_number)
  )

ggplot(sam_abies, aes(x = leaf_number))+
  geom_bar()

x <- sam_abies$leaf_number
mu <- mean(x, na.rm = TRUE)
sdx <- sd(x, na.rm = TRUE)
lo  <- mu - sdx
hi  <- mu + sdx

sam_cean <- sam_df %>% filter(host == "Ceanothus velutinus") %>% 
  select("leaf_number") %>% 
  mutate(
    leaf_number = as.numeric(leaf_number)
  )

ggplot(sam_cean, aes(x = leaf_number))+
  geom_bar()

x1 <- sam_cean$leaf_number
mu1 <- mean(x1, na.rm = TRUE)
sdx1 <- sd(x1, na.rm = TRUE)
lo1  <- mu1 - sdx1
hi1  <- mu1 + sdx1

# 1) create category + counts
sam_cat_df <- sam_df %>%
  mutate(
    grand_fir_growth = case_when(
      is.na(leaf_number)  ~ NA_character_,
      leaf_number < lo       ~ "low",
      leaf_number > hi       ~ "high",
      TRUE                ~ "medium"
    ),
    grand_fir_growth = factor(
      grand_fir_growth,
      levels = c("low", "medium", "high")
    ),
    ceanothus_growth = case_when(
      is.na(leaf_number)  ~ NA_character_,
      leaf_number < lo1      ~ "low",
      leaf_number > hi1      ~ "high",
      TRUE                ~ "medium"
    ),
    ceanothus_growth = factor(
      ceanothus_growth,
      levels = c("low", "medium", "high")
    )
  )

sample_data(bact) <- sample_data(sam_cat_df)

ps_filt <- bact %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  prune_taxa(taxa_sums(.) > 10, .)

ps <- bact

sample_data(ps)$inoculum_site <- factor(sample_data(ps)$inoculum_site)
sample_data(ps)$fire_freq <- factor(sample_data(ps)$fire_freq)

ps_ceanothus <- ps %>% subset_samples(host == "Ceanothus velutinus", na.rm = TRUE) %>% 
  subset_samples(!(ceanothus_growth == "medium")) %>%
  subset_samples(!is.na(ceanothus_growth))

ps_grandfir <- ps %>% subset_samples(host == "Abies grandis", na.rm = TRUE) %>% 
  subset_samples(!(grand_fir_growth == "medium")) %>%
  subset_samples(!is.na(grand_fir_growth))


a <- data.frame(sample_data(ps_ceanothus))

#corncob
fit_cean <- differentialTest(
  formula = ~ ceanothus_growth,
  phi.formula = ~ 1,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  test = "Wald",
  fdr_cutoff = 0.05,
  data = ps_ceanothus
)

fit_grand <- differentialTest(
  formula = ~ grand_fir_growth,
  phi.formula = ~ 1,
  formula_null = ~ 1,
  phi.formula_null = ~ 1,
  test = "Wald",
  fdr_cutoff = 0.05,
  data = ps_grandfir
)

plot(fit_cean)
plot(fit_grand)

rf_16S_cn <- readRDS("./R/shortt_additions/figures/random_forest/16S_snowbrush/top_genera_16sb.rds")
rf_16S_gr <- readRDS("./R/shortt_additions/figures/random_forest/16S_grandfir/top_genera_16gr.rds")

rf_16S_cn_g <- sub("^[^_]+_", "", rf_16S_cn)
rf_16S_gr_g <- sub("^[^_]+_", "", rf_16S_gr)

corncob_gr_genus <- corncob_sig_to_rf_genus(fit_grand, ps_grandfir)
corncob_cn_genus <- corncob_sig_to_rf_genus(fit_cean, ps_ceanothus)



# now compare directly to your RF genus vectors
overlap_gr_bact <- intersect(corncob_gr_genus, rf_16S_gr_g)
overlap_cn_bact <- intersect(corncob_cn_genus, rf_16S_cn_g)

overlap_gr
overlap_cn
overlap_gr_bact
overlap_cn_fung

# random_forest -----------------------------------------------------------

  



