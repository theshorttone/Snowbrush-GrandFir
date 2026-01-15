# -----------------------------------------------------------------------------#
# Pulling and processinf bacterial traits from BacDive database 
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.42.0
#                     patchwork v 1.1.2
#                     BacDive v 0.8.0
#                     lmerTest v 3.1.3
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
#  install.packages("BacDive", repos="http://R-Forge.R-project.org")

library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")
library(lmerTest); packageVersion("lmerTest")
source("./R/palettes.R")
source("./R/scale01.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]


clean_model_df <- function(x){
  broom.mixed::tidy(x) %>% 
    mutate(term=term %>% str_remove("indicator")) %>% 
    mutate(across(where(is.numeric),function(z){round(z,4)}))
}


# Random seed
set.seed(666)

readRenviron("./.Renviron")

# Data
bact <- readRDS("./Output/phyloseq_objects/16S_clean_phyloseq_object.RDS")

# reload point, for convenience
bact_trait_db <- readRDS("./Output/16S_Bacterial_Trait_Database.RDS")

pathogen_list <- readRDS("./Output/list_of_pathogenic_bacterial_genera.RDS")


pathogen_df <- readRDS("R/shortt_additions/outputs/pathogen_df.RDS")
mutualist_df <- readRDS("R/shortt_additions/outputs/mutualist_df.RDS")
guild_df <- readRDS("R/shortt_additions/outputs/guild_df.RDS")


# data manipulation --------------------------------------------------------
indicators <- c("wilting_scale","bud_number","leaf_number",
                "leaf_length","height","shoot_dm","final_root_dm")

guild_plot_df <- guild_df %>% 
  mutate(across(all_of(indicators), scale)) %>% 
  pivot_longer(all_of(indicators),
               names_to = "indicator",
               values_to = "value") %>% 
  mutate(
    indicator = indicator %>% str_replace_all("_"," ") %>% str_to_sentence(),
    
    proportion_mutualist = as.numeric(proportion_mutualist),
    proportion_pathogen  = as.numeric(proportion_pathogen),

    prop_mut_log10  = log10(proportion_mutualist + 1e-6),
    prop_path_log10 = log10(proportion_pathogen  + 1e-6),

    # nicer drought label if you want it
    Moisture = case_when(drought == "D"  ~ "Low",
                         drought == "ND" ~ "High",
                         TRUE ~ as.character(drought)),
    fire_freq = fct_na_value_to_level(factor(fire_freq), level = "control")
    
      
      )


# plotting  ---------------------------------------------------------------

plot_guild_effects <- function(df,
                               species_sel = c("GrandFir","Snowbrush"),
                               guild = c("mutualist","pathogen"),
                               color_by = c("drought","fire"),
                               log_x = TRUE) {
  
  species_sel <- match.arg(species_sel)
  guild    <- match.arg(guild)
  color_by <- match.arg(color_by)
  
  # choose x column
  x_col <- dplyr::case_when(
    guild == "mutualist"  & log_x ~ "prop_mut_log10",
    guild == "mutualist"  & !log_x ~ "proportion_mutualist",
    guild == "pathogen"   & log_x ~ "prop_path_log10",
    guild == "pathogen"   & !log_x ~ "proportion_pathogen",
    guild == "saprotroph" & log_x ~ "prop_sap_log10",
    guild == "saprotroph" & !log_x ~ "proportion_saprotroph"
  )
  
  df2 <- df %>% dplyr::filter(.data$species == species_sel)
  
  if (color_by == "drought") {
    col_aes <- rlang::sym("Moisture")   # or "drought" for D/ND
    color_scale <- scale_color_manual(values = drought_colors)
    col_title <- "Moisture"
  } else {
    df2 <- df2 %>% dplyr::mutate(fire_freq = ordered(fire_freq, levels = c("0","1","3")))
    col_aes <- rlang::sym("fire_freq")
    color_scale <- scale_color_manual(values = fire_colors)
    col_title <- "Fire frequency"
  }
  
  p <- ggplot(df2, aes(x = .data[[x_col]], y = value, color = !!col_aes)) +
    geom_point(alpha = .5) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~indicator, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 12)) +
    labs(
      title = species_sel,
      x = paste0("Proportion of ", guild, " bacteria", ifelse(log_x, " (log scale)", "")),
      y = "Scaled/Centered Value",
      color = col_title
    ) +
    color_scale
  
  # add log tick labels only when using prelogged x
  if (log_x) {
    p <- p + scale_x_continuous(
      breaks = log10(c(1e-6, 1e-4, 1e-2, 1e-1, 1)),
      labels = c("0", "0.0001", "0.01", "0.1", "1")
    )
  }
  
  p
}

p1 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="mutualist", color_by="drought", log_x=TRUE)
p2 <- plot_guild_effects(guild_plot_df, species_sel="GrandFir", guild="pathogen", color_by="drought", log_x=TRUE)
p3 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="mutualist", color_by="drought", log_x=TRUE)
p4 <- plot_guild_effects(guild_plot_df, species_sel="Snowbrush", guild="pathogen", color_by="drought", log_x=TRUE)

p1
p2
p3
p4

# combined growth ---------------------------------------------------------

growth_traits <- c("leaf_number","shoot_dm","final_root_dm")

eps <-  1e-6

guild_df_comp <- guild_df %>%
  dplyr::mutate(
    proportion_mutualist  = as.numeric(proportion_mutualist),
    proportion_pathogen   = as.numeric(proportion_pathogen),
    prop_mut_log10  = log10(proportion_mutualist  + eps),
    prop_path_log10 = log10(proportion_pathogen   + eps)
  ) %>%
  dplyr::mutate(across(all_of(growth_traits), scale01)) %>%
  dplyr::mutate(
    growth_index = rowMeans(dplyr::across(all_of(growth_traits)), na.rm = TRUE)
  )


# combined regressiong log ------------------------------------------------

plot_guild_growth_index_log <- function(df,
                                        species_sel = c("GrandFir","Snowbrush"),
                                        guild = c("mutualist","pathogen"),
                                        color_by = c("none","drought","fire_freq")) {
  
  species_sel <- match.arg(species_sel)
  guild       <- match.arg(guild)
  color_by    <- match.arg(color_by)
  
  if (guild == "mutualist") {
    lim = -6
  }else{
    lim = -1
  }
  
  if (lim == -6) {
    x_breaks <- c(-6, -4, -2, -1, 0)
    x_labels <- c("0", "0.0001", "0.01", "0.1", "1")
  } else { # lim == -1
    x_breaks <- c(-1, -0.5, 0)
    x_labels <- c("0.1", "0.3", "1")
  }
  
  df2 <- df %>%
    dplyr::filter(.data$species == species_sel) %>%
    dplyr::mutate(
      x_log = dplyr::case_when(
        guild == "mutualist"  ~ .data$prop_mut_log10,
        guild == "pathogen"   ~ .data$prop_path_log10
      ),
      Moisture = dplyr::case_when(
        drought == "D"  ~ "Low",
        drought == "ND" ~ "High",
        TRUE ~ as.character(drought)
      ),
      fire_freq = factor(fire_freq, levels = c("0","1","3"))
    )
  
  base <- ggplot2::ggplot(df2, ggplot2::aes(x = x_log, y = growth_index)) +
    ggplot2::geom_point(alpha = .5, size = 3) +
    ggplot2::geom_smooth(method = "lm", se = FALSE) +
    ggplot2::scale_x_continuous(
          limits = c(lim, 0),
          breaks = x_breaks,
          labels = x_labels,
        )+
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 16),
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = .5),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text  = ggplot2::element_text(face = "bold", size = 12)
    ) +
    ggplot2::labs(
      x = paste0("Proportion of ", guild, " bacteria (log scale)"),
      y = "Combined scaled growth",
      title = ifelse(species_sel == "GrandFir", "Grand Fir", "Snowbrush")
    )
  
  if (color_by == "none") {
    return(base)
  }
  
  if (color_by == "drought") {
    return(
      base +
        ggplot2::aes(color = Moisture) +
        ggplot2::scale_color_manual(values = drought_colors) +
        ggplot2::labs(color = "Moisture")
        
    )
  }
  # color_by == "fire_freq"
  base +
    ggplot2::aes(color = fire_freq) +
    ggplot2::scale_color_manual(values = fire_colors) +
    ggplot2::labs(color = "Fire frequency")
}



# Grand fir, pathogens, color by drought
p1 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "pathogen", "fire_freq")
p2 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "mutualist", "fire_freq")
p3 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "pathogen", "fire_freq")
p4 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "mutualist", "fire_freq")
p1
combined_growth <- (p1 + p2) / (p3 + p4)
combined_growth
ggsave("./R/shortt_additions/figures/16S_combined_growth_fire.png",
       plot = combined_growth,
       width = 10, height = 8, units = "in", dpi = 300)
p1 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "pathogen", "drought")
p2 <- plot_guild_growth_index_log(guild_df_comp, "GrandFir", "mutualist", "drought")
p3 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "pathogen", "drought")
p4 <- plot_guild_growth_index_log(guild_df_comp, "Snowbrush", "mutualist", "drought")
p1
combined_growth <- (p1 + p2) / (p3 + p4)
combined_growth
ggsave("./R/shortt_additions/figures/16S_combined_growth_drought.png",
       plot = combined_growth,
       width = 10, height = 8, units = "in", dpi = 300)

# results_drought <- make_results_df(guild_df_comp, interact = "drought")
# results_fire    <- make_results_df(guild_df2, interact = "fire_freq")

# bar graphs  -------------------------------------------------------------

# proportion mutualist bar graph ------------------------------------------

bar_plot_fire_df <- guild_df_comp %>% 
  group_by(species, fire_freq) %>%
  summarise(
    mean_mutualist = mean(proportion_mutualist, na.rm = TRUE),
    sd_mutualist   = sd(proportion_mutualist, na.rm = TRUE),
    n              = n(),
    se_mutualist   = sd_mutualist / sqrt(n),
    mean_pathogen = mean(proportion_pathogen, na.rm = TRUE),
    sd_pathogen   = sd(proportion_pathogen, na.rm = TRUE),
    se_pathogen   = sd_pathogen / sqrt(n)
  ) %>% 
  mutate(
    fire_freq = factor(fire_freq),
  )

bar_plot_drought_df <- guild_df_comp %>% 
  group_by(species, drought) %>%
  summarise(
    mean_mutualist = mean(proportion_mutualist, na.rm = TRUE),
    sd_mutualist   = sd(proportion_mutualist, na.rm = TRUE),
    n              = n(),
    se_mutualist   = sd_mutualist / sqrt(n),
    mean_pathogen = mean(proportion_pathogen, na.rm = TRUE),
    sd_pathogen   = sd(proportion_pathogen, na.rm = TRUE),
    se_pathogen   = sd_pathogen / sqrt(n)
  )

# plots -------------------------------------------------------------------

p1 <- ggplot(bar_plot_fire_df,
             aes(x = fire_freq, y = mean_mutualist, fill = fire_freq)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_mutualist - se_mutualist,
                    ymax = mean_mutualist + se_mutualist),
                width = .2,
                position = position_dodge(.9)) +
  facet_wrap(~species) +
  scale_fill_manual(values = fire_colors) +
  labs(
    x = "Fire frequency",
    y = "Mean proportion of mutualist bacteria"
  )+
  theme_few()

p2 <- ggplot(bar_plot_fire_df,
             aes(x = fire_freq, y = mean_pathogen, fill = fire_freq)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_pathogen - se_pathogen,
                    ymax = mean_pathogen + se_pathogen),
                width = .2,
                position = position_dodge(.9)) +
  facet_wrap(~species) +
  scale_fill_manual(values = fire_colors) +
  labs(
    x = "Fire frequency",
    y = "Mean proportion of pathogen bacteria"
  )+
  theme_few()

p3 <- ggplot(bar_plot_drought_df,
             aes(x = drought, y = mean_pathogen, fill = drought)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_pathogen - se_pathogen,
                    ymax = mean_pathogen + se_pathogen),
                width = .2,
                position = position_dodge(.9)) +
  facet_wrap(~species) +
  scale_fill_manual(values = drought_colors) +
  labs(
    x = "Fire frequency",
    y = "Mean proportion of pathogen bacteria"
  )+
  theme_few()

p4 <- ggplot(bar_plot_drought_df,
             aes(x = drought, y = mean_mutualist, fill = drought)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_mutualist - se_mutualist,
                    ymax = mean_mutualist + se_mutualist),
                width = .2,
                position = position_dodge(.9)) +
  facet_wrap(~species) +
  scale_fill_manual(values = drought_colors) +
  labs(
    x = "Fire frequency",
    y = "Mean proportion of mutualist bacteria"
  )+
  theme_few()

p4

combined_bar_plot_fire <- (p1 | p2)


combined_bar_plot_fire
ggsave(combined_bar_plot_fire,
       filename = "R/shortt_additions/figures/mean_bacteria_proportion_guilds_fire_bar.pdf",
       width = 10,
       height = 4,
       units = "in",
       dpi = 300)

combined_bar_plot_drought <- (p4 | p3)
combined_bar_plot_drought
ggsave(combined_bar_plot_drought,
       filename = "R/shortt_additions/figures/mean_proportion_bacteria_guilds_drought_bar.pdf",
       width = 10,
       height = 4,
       units = "in",
       dpi = 300)

