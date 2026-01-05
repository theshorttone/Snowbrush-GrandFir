library(dplyr)
library(tibble)

site1.inoc.full <- readRDS("Output/phyloseq_objects/18S_site1.inoc.full.RDS") %>% condense_ps_to_species()
site2.inoc.full <- readRDS("Output/phyloseq_objects/18S_site2.inoc.full.RDS") %>% condense_ps_to_species()
site3.inoc.full <- readRDS("Output/phyloseq_objects/18S_site3.inoc.full.RDS") %>% condense_ps_to_species()
site4.inoc.full <- readRDS("Output/phyloseq_objects/18S_site4.inoc.full.RDS") %>% condense_ps_to_species()
site5.inoc.full <- readRDS("Output/phyloseq_objects/18S_site5.inoc.full.RDS") %>% condense_ps_to_species()
site6.inoc.full <- readRDS("Output/phyloseq_objects/18S_site6.inoc.full.RDS") %>% condense_ps_to_species()

inoc.full.am <- merge_phyloseq(site1.inoc.full,site2.inoc.full,site3.inoc.full,
                               site4.inoc.full,site5.inoc.full,site6.inoc.full)
inoc.full.am <- prune_samples(sample_sums(inoc.full.am) > 0, inoc.full.am)

df_samples <- sample_data(inoc.full.am) 
df_samples <- data.frame(df_samples)   

df_samples <- clean_func(df_samples)
sample_data(inoc.full.am) <- sample_data(df_samples)
## --- raw OTU table ---
otu_raw <- as(otu_table(inoc.full.am), "matrix")
if (taxa_are_rows(inoc.full.am)) otu_raw <- t(otu_raw)

## --- relative abundance (safe) ---
ps_rel <- prune_samples(sample_sums(inoc.full.am) > 0, inoc.full.am) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) x / sum(x))

otu_rel <- as(otu_table(ps_rel), "matrix")
if (taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)

## --- NMDS (use the ord you already ran if you want) ---
ord <- ordinate(ps_rel, method = "NMDS", distance = "bray", trymax = 100)

pts <- as.data.frame(vegan::scores(ord, display = "sites"))
colnames(pts)[1:2] <- c("NMDS1", "NMDS2")

## --- diagnostics table ---
qc_am <- tibble(
  sample = rownames(otu_raw),
  depth_raw = sample_sums(inoc.full.am),
  richness = rowSums(otu_raw > 0),
  top_taxon_prop = apply(otu_rel, 1, max, na.rm = TRUE)
) %>%
  left_join(
    pts %>%
      mutate(sample = rownames(.),
             nmds_r = sqrt(NMDS1^2 + NMDS2^2)),
    by = "sample"
  ) %>%
  arrange(desc(nmds_r))
qc_am <- qc_am %>%
  mutate(
    flag_low_depth   = depth_raw < 50,
    flag_low_rich    = richness < 3,
    flag_dominant    = top_taxon_prop > 0.9,
    flag_problematic = flag_low_depth | flag_low_rich
  )
table(qc_am$flag_problematic)

bad_samples <- qc_am %>%
  filter(flag_problematic) %>%
  pull(sample)

ps_am_clean <- prune_samples(
  !sample_names(inoc.full.am) %in% bad_samples,
  inoc.full.am
)
