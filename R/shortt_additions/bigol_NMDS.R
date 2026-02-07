library(tidyverse)
library(phyloseq)
library(zahntools)
library(vegan); packageVersion("vegan")
library(ggnewscale)
library(ggthemes)

theme_set(theme_minimal())
source("./R/palettes.R")
site_colors <- c(
  "Site 1" = "#4cbfe6", "Site 2" = "#2443f0",
  "Site 3" = "#f0c424", "Site 4" = "#db7e04",
  "Site 5" = "#cc6866", "Site 6" = "#9e0402"
)

# 3 mid-blend colors (averages) for ellipses
burn_colors <- c(
  "0" = "#3881F3",
  "1"   = "#E59D14",
  "3"   = "#A73634"
)


# plotting functions ------------------------------------------------------

#via Geooff:
condense_ps_to_species <- function(x){
  # remove bad taxa assignments
  x <- x %>% subset_taxa(!is.na(Phylum))
  # build data frame
  sp <- x@tax_table[,7] %>% as.character()
  gn <- x@tax_table[,6] %>% as.character()
  or <- x@tax_table[,4] %>% as.character()
  sb <- x@tax_table[,3] %>% as.character()
  ph <- x@tax_table[,2] %>% as.character()
  
  condensed_taxonomy <- 
    data.frame(ph,sb,or,gn,sp) %>% 
    mutate(spp = case_when(is.na(or) & is.na(gn) & is.na(sp) ~ paste0(sb," sp."),
                           !is.na(or) & is.na(gn) & is.na(sp) ~ paste0(or," sp."),
                           !is.na(or) & !is.na(gn) & is.na(sp) ~ paste0(gn," sp."),
                           !is.na(or) & !is.na(gn) & !is.na(sp) ~ sp))
  
  x@tax_table[,7] <- str_replace_all(condensed_taxonomy$spp,"_"," ")
  x <- tax_glom(x,"Species")
  return(x)
}

make_ord_plot <- function(ps,
                          method = "NMDS",
                          distance = "bray",
                          shape_var = "community",
                          point_color_var = "plot_nmds_color",
                          ellipse_color_var = "fire_freq",
                          ellipse_linetype_var = "community",
                          site_colors,
                          burn_colors,
                          ellipse_level = 0.8,
                          ellipse_type = "t",
                          trymax = 20) {
  
  
  # ---- ordination ----
  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  
  ord <- if (toupper(method) == "NMDS") {
    phyloseq::ordinate(ps_rel, method = "NMDS", distance = distance, trymax = trymax)} else {
    phyloseq::ordinate(ps_rel, method = method, distance = distance)
  }

  # ---- base plot (points) ----
  p <- phyloseq::plot_ordination(
    ps, ord,
    shape = shape_var,
    color = point_color_var,
    
  ) +
    scale_color_manual(values = site_colors,
                       name = "Inoculation site")+
    scale_shape_manual(
      values = c(16, 17, 11),
      name = "Community",
    )
  
  # ---- second color scale (ellipses) ----
  p <- p +
    new_scale_color() +
    stat_ellipse(
      aes(
        color = factor(.data[[ellipse_color_var]]),
        linetype = .data[[ellipse_linetype_var]],
        group = interaction(factor(.data[[ellipse_color_var]]), .data[[ellipse_linetype_var]])
      ),
      type = ellipse_type,
      level = ellipse_level,
      linewidth = 0.7
    ) +
    scale_color_manual(values = burn_colors, name = "Fire frequency") +
    scale_linetype_manual(values = c("solid", "longdash", "dotted"))+
    guides(
      color = guide_legend(order = 1),
      fill  = guide_legend(order = 2),
      linetype = "none"
    )+
    theme_few()
  
  
  return(list(ord = ord, plot = p))
}


# data cleaning function --------------------------------------------------


clean_func <- function (df){
  df_samples <- df %>%
    mutate(
      plot_nmds_color = case_when(
        !is.na(.data$inoculum) & .data$inoculum != "" ~ str_replace(.data$inoculum, "^.*_", ""),
        TRUE ~ recode(
          .data$other_frompreviouscolumn,
          "Site1" = "Site 1",
          "Site2" = "Site 2",
          "Site3" = "Site 3",
          "Site4" = "Site 4",
          "Site5" = "Site 5",
          "Site6" = "Site 6",
          .default = NA_character_
        )),
      plot_nmds_color = recode(plot_nmds_color,
                               "W1" = "Site 1",
                               "W2" = "Site 2",
                               "W3" = "Site 3",
                               "W4" = "Site 4",
                               "W5" = "Site 5",
                               "W6" = "Site 6",
                               .default = .data$plot_nmds_color
      ),
      fire_freq = if_else(
        .data$other_frompreviouscolumn %in% paste0("Site", 1:6),
        case_when(
          .data$other_frompreviouscolumn %in% c("Site1", "Site2") ~ 0,
          .data$other_frompreviouscolumn %in% c("Site3", "Site4") ~ 1,
          .data$other_frompreviouscolumn %in% c("Site5", "Site6") ~ 3
        ),
        as.numeric(.data$fire_freq)
      )
    )
}
# 16S ---------------------------------------------------------------------


site1.inoc.full <- readRDS("Output/phyloseq_objects/16S_site1.inoc.full.RDS")
site2.inoc.full <- readRDS("Output/phyloseq_objects/16S_site2.inoc.full.RDS")
site3.inoc.full <- readRDS("Output/phyloseq_objects/16S_site3.inoc.full.RDS")
site4.inoc.full <- readRDS("Output/phyloseq_objects/16S_site4.inoc.full.RDS")
site5.inoc.full <- readRDS("Output/phyloseq_objects/16S_site5.inoc.full.RDS")
site6.inoc.full <- readRDS("Output/phyloseq_objects/16S_site6.inoc.full.RDS")
ps <- readRDS("Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")

# merge all inoculum phyloseq objects
inoc.full <- merge_phyloseq(site1.inoc.full,site2.inoc.full,site3.inoc.full,
                           site4.inoc.full,site5.inoc.full,site6.inoc.full)

df_samples <- sample_data(inoc.full) 
df_samples <- data.frame(df_samples)   
  
df_samples <- clean_func(df_samples) %>% 
  mutate(
    community2 = dplyr::case_when(
      species == "Snowbrush" ~ "final snowbrush",
      species == "GrandFir"  ~ "final grandfir",
      TRUE ~ community
    )
  )

sample_data(inoc.full) <- sample_data(df_samples)


make_ord_plot(
  ps = inoc.full,
  method = "NMDS",
  distance = "bray",
  shape_var = "community2",
  point_color_var = "plot_nmds_color",
  ellipse_color_var = "fire_freq",
  ellipse_linetype_var = "community2",
  site_colors = site_colors,
  burn_colors = burn_colors,
  ellipse_level = 0.8,
  ellipse_type = "t",
  trymax = 20
)$plot +
  labs(
    title = "NMDS of 16S gene sequences Stress = 0.15",
    subtitle = " "
)
ggsave("R/shortt_additions/figures/bigol_16S_NMDS_inoculum.png",
       width = 6, height = 5)


# Fungal-----------------------------------------------------------------------
site1.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site1.inoc.full.RDS")
site2.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site2.inoc.full.RDS")
site3.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site3.inoc.full.RDS")
site4.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site4.inoc.full.RDS")
site5.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site5.inoc.full.RDS")
site6.inoc.full <- readRDS("Output/phyloseq_objects/ITS_site6.inoc.full.RDS")

inoc.full.fung <- merge_phyloseq(site1.inoc.full,site2.inoc.full,site3.inoc.full,
                            site4.inoc.full,site5.inoc.full,site6.inoc.full)

df_samples <- sample_data(inoc.full.fung) 
df_samples <- data.frame(df_samples)   

df_samples <- clean_func(df_samples) %>% 
  mutate(
    community2 = dplyr::case_when(
      species == "Snowbrush" ~ "final snowbrush",
      species == "GrandFir"  ~ "final grandfir",
      TRUE ~ community
    )
  )
  
sample_data(inoc.full.fung) <- sample_data(df_samples)

make_ord_plot(
  ps = inoc.full.fung,
  method = "NMDS",
  distance = "bray",
  shape_var = "community2",
  point_color_var = "plot_nmds_color",
  ellipse_color_var = "fire_freq",
  ellipse_linetype_var = "community2",  site_colors = site_colors,
  burn_colors = burn_colors,
  trymax = 20
)$plot +
  labs(
    title = "NMDS of ITS gene sequences but STRESS = 0.27 (bad)",
    subtitle = " "
  )
ggsave("R/shortt_additions/figures/bigol_ITS_DCA_inoculum.png",
       width = 6, height = 5)


# AM (18S)---------------------------------------------------------------------

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


am_res <- make_ord_plot(
  ps = ps_am_clean,
  method = "NMDS",
  distance = "bray",
  site_colors = site_colors,
  burn_colors = burn_colors,
  trymax = 20
)



names(am_res)      # should show e.g. "plot" and "ord" (depending on your function)
p_am   <- am_res$plot + labs(title = "AM of 18S gene sequences", subtitle = " ")
ord_am <- am_res$ord
ggsave("R/shortt_additions/figures/bigol_AM_NMDS_inoculum.png",
       width = 6, height = 5)
p_am


pts <- as.data.frame(ord_am$points)
colnames(pts)[1:2] <- c("NMDS1","NMDS2")

# Rank by distance from origin (big = likely problematic)
outliers_nmds <- pts %>%
  mutate(sample = rownames(.),
         r = sqrt(NMDS1^2 + NMDS2^2)) %>%
  arrange(desc(r))

outliers_nmds %>% head(15)

