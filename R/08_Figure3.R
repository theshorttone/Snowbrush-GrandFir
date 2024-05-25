library(ggplot2,lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.3/")
packageVersion("ggplot2")
library(tidyverse)
library(patchwork)
library(janitor)
library(lmerTest)
library(broom)
# # functions ####
# scales_add_defaults <- function(scales, data, aesthetics, env) {
#   if (is.null(aesthetics)) return()
#   names(aesthetics) <- unlist(lapply(names(aesthetics), aes_to_scale))
#   
#   new_aesthetics <- setdiff(names(aesthetics), scales$input())
#   # No new aesthetics, so no new scales to add
#   if (is.null(new_aesthetics)) return()
#   
#   datacols <- lapply(aesthetics[new_aesthetics], rlang::eval_tidy, data = data)
#   datacols <- compact(datacols)
#   
#   for (aes in names(datacols)) {
#     scales$add(find_scale(aes, datacols[[aes]], env))
#   }
#   
# }
# aes_to_scale <- function(var) {
#   var[var %in% ggplot_global$x_aes] <- "x"
#   var[var %in% ggplot_global$y_aes] <- "y"
#   
#   var
# }
# ggplot_global <- new.env(parent = emptyenv())
# 
# # The current theme. Defined here only as placeholder, and defined properly
# # in file "theme-current.R". This setup avoids circular dependencies among
# # the various source files.
# ggplot_global$theme_current <- list()
# 
# # Element tree for the theme elements. Defined here only as placeholder, and
# # defined properly in file "theme-elements.r".
# ggplot_global$element_tree <- list()
# 
# # List of all aesthetics known to ggplot
# # (In the future, .all_aesthetics should be removed in favor
# # of direct assignment to ggplot_global$all_aesthetics, see below.)
# .all_aesthetics <- c(
#   "adj", "alpha", "angle", "bg", "cex", "col", "color",
#   "colour", "fg", "fill", "group", "hjust", "label", "linetype", "lower",
#   "lty", "lwd", "max", "middle", "min", "pch", "radius", "sample", "shape",
#   "size", "srt", "upper", "vjust", "weight", "width", "x", "xend", "xmax",
#   "xmin", "xintercept", "y", "yend", "ymax", "ymin", "yintercept", "z"
# )
# 
# ggplot_global$all_aesthetics <- .all_aesthetics
# 
# # Aesthetic aliases
# # (In the future, .base_to_ggplot should be removed in favor
# # of direct assignment to ggplot_global$base_to_ggplot, see below.)
# .base_to_ggplot <- c(
#   "col"   = "colour",
#   "color" = "colour",
#   "pch"   = "shape",
#   "cex"   = "size",
#   "lty"   = "linetype",
#   "lwd"   = "linewidth",
#   "srt"   = "angle",
#   "adj"   = "hjust",
#   "bg"    = "fill",
#   "fg"    = "colour",
#   "min"   = "ymin",
#   "max"   = "ymax"
# )
# 
# ggplot_global$base_to_ggplot <- .base_to_ggplot
# 
# # These two vectors must match in length and position of symmetrical aesthetics
# # xintercept2 is a filler to match to the intercept aesthetic in geom_abline
# ggplot_global$x_aes <- c("x", "xmin", "xmax", "xend", "xintercept",
#                          "xmin_final", "xmax_final", "xlower", "xmiddle", "xupper", "x0")
# 
# ggplot_global$y_aes <- c("y", "ymin", "ymax", "yend", "yintercept",
#                          "ymin_final", "ymax_final", "lower", "middle", "upper", "y0")


source("./R/palettes.R")
drought_colors <- pal.discrete[c(2,5)]
host_colors <- pal.discrete[c(7,10)] 
fire_colors <- pal.discrete[c(18,2,14)]
inoc_colors <- c("#4cbfe6","#2443f0",
                 "#f0c424","#db7e04",
                 "#cc6866","#9e0402",
                 "gray20")
# data ####
# load plant time series data
plants <- read_csv("./Data/MasterMeasurement_Sep_Sheet1.csv") %>% clean_names()

# load ordination plots
grandfir_bact_ord_plot <- readRDS("./Output/grandfir_bact_ord_plot.RDS") + theme(axis.title = element_blank()) + labs(y="Bacteria") + theme(axis.title.y = element_text(face='bold',size=14))
snowbrush_bact_ord_plot <- readRDS("./Output/snowbrush_bact_ord_plot.RDS") + theme(axis.title.y = element_text(face='bold',size=14),axis.title.x = element_blank())
grandfir_fungi_ord_plot <- readRDS("./Output/grandfir_fungi_ord_plot.RDS")  + theme(axis.title = element_blank()) + labs(y="Fungi") + theme(axis.title.y = element_text(face='bold',size=14))
snowbrush_fungi_ord_plot <- readRDS("./Output/snowbrush_fungi_ord_plot.RDS") + theme(axis.title.y = element_text(face='bold',size=14),axis.title.x = element_blank())
snowbrush_18S_ord_plot <- readRDS("./Output/snowbrush_18S_ord_plot.RDS") + theme(axis.title = element_blank())


# clean plant data
names(plants)
plants <- 
plants %>% 
  dplyr::select(all_of(c("unique_id","block","unique_id_2","block_2","seed_species","drought_trt","soil_inoculum","burn_frequency")),
                contains("leaf_num")) %>% 
  pivot_longer(contains("leaf_num"),names_to = "month",values_to = "leaf_number") %>% 
  mutate(month = month %>% str_remove("tru.*_leaf_num_") %>% str_to_sentence() %>% str_replace("June","Jun")) %>% 
  mutate(month = factor(month, levels = c("Nov","Dec","Jan","Feb","Mar","May","Jun","Aug"))) %>% 
  # mutate(plantgrowth = scale(leaf_number)[,1]) %>% 
  mutate(host = case_when(seed_species == "GrandFir" ~ "Abies grandis",
                          seed_species == "Snowbrush" ~ "Ceanothus velutinus"))

scaled_plantgrowth <- 
plants %>% 
  group_by(host) %>% 
  summarize(plantgrowth = scale(leaf_number)) %>% 
  pluck("plantgrowth")
plants$scaled_plantgrowth <- c(scaled_plantgrowth) %>% as.numeric()

plant_means <- 
plants %>% 
  group_by(host,drought_trt,month) %>% 
  summarize(plantgrowth = mean(scale(leaf_number)[,1]))

plants$inoc_source <-   
ifelse(grepl("_W",plants$soil_inoculum),plants$soil_inoculum,"_Sterile") %>% 
  str_split("_") %>% 
  map_chr(2) %>% 
  str_remove("W")


# plot plant growth over time
plants %>% 
  dplyr::filter(!is.na(host)) %>% 
  ggplot(aes(x=month,y=scaled_plantgrowth,color=drought_trt)) +
  geom_jitter(size=3,width = .2,alpha=.25) +
  geom_boxplot(alpha=.2,size=1.5,fill='black') +
  facet_wrap(~host,ncol = 2) +
  scale_color_manual(values = drought_colors,labels=c("Low","High")) +
  labs(color="Moisture",x="\nTimepoint",y="Scaled plant growth") +
  theme_bw() +
  theme(strip.text = element_text(face='bold.italic',size=18),
        strip.background = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5))
ggsave("./Output/figs/plantgrowth_time_series_host.png",dpi=300,height = 6,width = 10)


plants %>% 
  dplyr::filter(!is.na(host)) %>% 
  ggplot(aes(x=month,y=scaled_plantgrowth,color=drought_trt)) +
  geom_jitter(size=3,width = .2,alpha=.25) +
  geom_boxplot(alpha=.2,size=1.5,fill='black') +
  facet_wrap(~inoc_source) +
  scale_color_manual(values = drought_colors,labels=c("Low","High")) +
  labs(color="Moisture",x="\nTimepoint",y="Scaled plant growth") +
  theme_bw() +
  theme(strip.text = element_text(face='bold.italic',size=18),
        strip.background = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5))
ggsave("./Output/figs/plantgrowth_time_series_inoc.png",dpi=300,height = 6,width = 10)
plants$host
plants %>% 
  dplyr::filter(!is.na(host) & host == "Abies grandis") %>% 
  ggplot(aes(x=month,y=scaled_plantgrowth,color=drought_trt)) +
  geom_jitter(size=3,width = .2,alpha=.25) +
  geom_boxplot(alpha=.2,size=1.5,fill='black') +
  facet_wrap(~inoc_source) +
  scale_color_manual(values = drought_colors,labels=c("Low","High")) +
  labs(color="Moisture",x="\nTimepoint",y="Scaled plant growth") +
  theme_bw() +
  theme(strip.text = element_text(face='bold.italic',size=18),
        strip.background = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5))

plants %>% 
  dplyr::filter(!is.na(host) & host == "Ceanothus velutinus") %>% 
  ggplot(aes(x=month,y=scaled_plantgrowth,color=drought_trt,fill=drought_trt)) +
  geom_jitter(size=1,width = .2,alpha=.25) +
  geom_boxplot(alpha=.5,size=1.5) +
  facet_wrap(~inoc_source) +
  scale_color_manual(values = drought_colors,labels=c("Low","High")) +
  scale_fill_manual(values = drought_colors,guide='none') +
  labs(color="Moisture",x="\nTimepoint",y="Scaled plant growth") +
  theme_bw() +
  theme(strip.text = element_text(face='bold.italic',size=18),
        strip.background = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(face='bold',size=10,angle=90,hjust=1,vjust=.5))


ordplots <- 
(grandfir_fungi_ord_plot + snowbrush_fungi_ord_plot) / (grandfir_bact_ord_plot + snowbrush_bact_ord_plot) +
  plot_layout(guides = 'collect') # + plot_annotation(tag_levels = "A")


# rework...No idea where the code to build these plots went!?

# make backup copies
grandfir_fungi_ord_plot2 <- grandfir_fungi_ord_plot 
snowbrush_fungi_ord_plot2 <- snowbrush_fungi_ord_plot
grandfir_bact_ord_plot2 <- grandfir_bact_ord_plot 
snowbrush_bact_ord_plot2 <- snowbrush_bact_ord_plot 

# remove the ellipses
grandfir_fungi_ord_plot2$layers[[2]] <- NULL 
snowbrush_fungi_ord_plot2$layers[[2]] <- NULL
grandfir_bact_ord_plot2$layers[[2]] <- NULL
snowbrush_bact_ord_plot2$layers[[2]] <- NULL

# add new columns to data for separate ellipses
grandfir_fungi_ord_plot2$data$newgroup <- 
  ifelse(grandfir_fungi_ord_plot2$data$inoculum_site == "Sterile","Uninoculated","Inoculated")
snowbrush_fungi_ord_plot2$data$newgroup <- 
  ifelse(snowbrush_fungi_ord_plot2$data$inoculum_site == "Sterile","Uninoculated","Inoculated")
grandfir_bact_ord_plot2$data$newgroup <- 
  ifelse(grandfir_bact_ord_plot2$data$inoculum_site == "Sterile","Uninoculated","Inoculated")
snowbrush_bact_ord_plot2$data$newgroup <- 
  ifelse(snowbrush_bact_ord_plot2$data$inoculum_site == "Sterile","Uninoculated","Inoculated")

grandfir_fungi_ord_plot2$data$inoculum_site %>% unique()

# change color scheme to be colorblind friendly
grandfir_fungi_ord_plot2 <- grandfir_fungi_ord_plot2 +
  stat_ellipse(aes(linetype=newgroup,group=newgroup)) +
  scale_linetype_manual(values = c(1,212)) +
  labs(linetype = "Inoculation") +
  scale_color_manual(values = inoc_colors, labels = as.character(c(1:6,"Uninoculated"))) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.title.y = element_text(face='bold',size=18),legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14))

snowbrush_fungi_ord_plot2 <- snowbrush_fungi_ord_plot2 +
  stat_ellipse(aes(linetype=newgroup,group=newgroup)) +
  scale_linetype_manual(values = c(1,212)) +
  labs(linetype = "Inoculation") +
  scale_color_manual(values = inoc_colors, labels = as.character(c(1:6,"Uninoculated"))) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(face='bold',size=16), axis.title.y = element_blank(),
        legend.text = element_text(face='bold',size=14))

grandfir_bact_ord_plot2 <- grandfir_bact_ord_plot2 +
  stat_ellipse(aes(linetype=newgroup,group=newgroup)) +
  scale_linetype_manual(values = c(1,212)) +
  labs(linetype = "Inoculation", title = "") +
  scale_color_manual(values = inoc_colors, labels = as.character(c(1:6,"Uninoculated"))) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
        axis.title.y = element_text(face='bold',size=18),legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14))

snowbrush_bact_ord_plot2 <- snowbrush_bact_ord_plot2 +
  stat_ellipse(aes(linetype=newgroup,group=newgroup)) +
  scale_linetype_manual(values = c(1,212)) +
  labs(linetype = "Inoculation", title = "") +
  scale_color_manual(values = inoc_colors, labels = as.character(c(1:6,"Uninoculated"))) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank(),
       legend.title = element_text(face='bold',size=16), axis.title.y = element_blank(),
       legend.text = element_text(face='bold',size=14))



# see plot without ALL ellipses
ordplots <- (grandfir_fungi_ord_plot2 + snowbrush_fungi_ord_plot2) / (grandfir_bact_ord_plot2 + snowbrush_bact_ord_plot2) +
  plot_layout(guides = 'collect') # + plot_annotation(tag_levels = "A")


# make top half of figure (boxplots of plant growth)


# first, find significant differences
aug <- plants %>% 
  dplyr::filter(month == "Aug" & !is.na(unique_id))

gf <- 
  aug %>% 
  dplyr::filter(host == "Abies grandis")
sb <- 
  aug %>% 
  dplyr::filter(host == "Ceanothus velutinus")

# gf_lmer <- 
#   lmer(data = gf, formula = leaf_number ~ drought_trt + (1|inoc_source))
# summary(gf_lmer)

# sb_lmer <- 
#   lmer(data = sb, formula = leaf_number ~ drought_trt + (1|inoc_source))
# summary(sb_lmer)
# 
# tuk <- aov(data=sb,formula = leaf_number ~ drought_trt * inoc_source) %>% 
#   TukeyHSD()
# tuk$`drought_trt:inoc_source` %>% 
#   as.data.frame() %>% 
#   dplyr::filter(`p adj` < 0.05)




plants %>% names
plants$drought_trt

grandfir_plants_plot <- 
  plants %>% 
  mutate(inoc_source = case_when(inoc_source == "Sterile" ~ "Uninoculated",TRUE ~ inoc_source),
         burn_frequency = paste0(burn_frequency," Burn"),
         Moisture = case_when(drought_trt == "D" ~ "Low",drought_trt == "ND" ~ "High")) %>%
  mutate(combo = paste0(inoc_source,"_",Moisture)) %>% 
  dplyr::filter(seed_species == "GrandFir") %>% 
  dplyr::filter(month == "Aug") %>% 
  ggplot(aes(x=combo,y=scaled_plantgrowth,fill=inoc_source,color=inoc_source)) +
  geom_jitter(width = .1,alpha=.1) +
  geom_boxplot(color='black',linewidth=.2,alpha=c(rep(c(1,0.2),7))) +
  scale_color_manual(values = rep(inoc_colors,2),guide='none') +
  scale_fill_manual(values = rep(inoc_colors,2)) +
  geom_vline(linetype=221,xintercept = c(4.5,8.5,12.5)) +
  labs(fill="Inoculum source",
       y="Scaled\nplant growth",
       title = "") +
  annotate('text',x=2.5,y=6,label="0 Burn",size=6,fontface='bold') +
  annotate('text',x=6.5,y=6,label="1 Burn",size=6,fontface='bold') +
  annotate('text',x=10.5,y=6,label="3 Burn",size=6,fontface='bold') +
  coord_cartesian(ylim=c(min(plants$scaled_plantgrowth,na.rm=TRUE),6.5)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title.y = element_text(face='bold',size=18),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14,face='bold'),
        plot.title = element_text(face='bold.italic',size=18,hjust=.5),
        legend.position = 'none') 

snowbrush_plants_plot <- 
  plants %>% 
  mutate(inoc_source = case_when(inoc_source == "Sterile" ~ "Uninoculated",TRUE ~ inoc_source),
         burn_frequency = paste0(burn_frequency," Burn"),
         Moisture = case_when(drought_trt == "D" ~ "Low",drought_trt == "ND" ~ "High")) %>%
  mutate(combo = paste0(inoc_source,"_",Moisture)) %>% 
  dplyr::filter(seed_species == "Snowbrush") %>% 
  dplyr::filter(month == "Aug") %>% 
  ggplot(aes(x=combo,y=scaled_plantgrowth,fill=inoc_source,color=inoc_source)) +
  geom_jitter(width = .1,alpha=.1) +
  geom_boxplot(color='black',linewidth=.2,alpha=c(rep(c(1,0.2),7))) +
  scale_color_manual(values = rep(inoc_colors,2),guide='none') +
  scale_fill_manual(values = rep(inoc_colors,2)) +
  geom_vline(linetype=221,xintercept = c(4.5,8.5,12.5)) +
  labs(fill="Inoculum source",
       y="Scaled\nplant growth",
       title = "") +
  annotate('text',x=2.5,y=6,label="0 Burn",size=6,fontface='bold') +
  annotate('text',x=6.5,y=6,label="1 Burn",size=6,fontface='bold') +
  annotate('text',x=10.5,y=6,label="3 Burn",size=6,fontface='bold') +
  coord_cartesian(ylim=c(min(plants$scaled_plantgrowth,na.rm=TRUE),6.5)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        legend.position = 'none',
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(size=14,face='bold'),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face='bold.italic',size=18,hjust=.5)) 


barplots <- grandfir_plants_plot | snowbrush_plants_plot 


showtext::showtext_auto()
# showtext::showtext_opts(dpi = 96)
ordplots / barplots + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold',size=18,hjust = .5,vjust = 0))

ggsave("./Output/figs/Figure_3.png",dpi=300,height = 8, width = 12)
ggsave("./Output/figs/Figure_3.tiff",dpi=300,height = 8, width = 12)

m <- 
gf %>% 
  dplyr::filter(!is.na(burn_frequency)) %>% 
  aov(data=.,
       formula = leaf_number ~ drought_trt + (1|as.numeric(inoc_source)))
summary(m)

m <- 
  sb %>% 
  dplyr::filter(!is.na(burn_frequency)) %>% 
  aov(data=.,
      formula = leaf_number ~ drought_trt + (1|as.numeric(inoc_source)))
summary(m)

# load permanova tables for associated figure stats
perm_bact <- readRDS("./Output/16S_Permanova_Table.RDS") %>% 
  mutate(term=term %>% str_replace("drought","moisture"),amplicon="16S") %>% 
  dplyr::select(amplicon,everything())
perm_fung <- readRDS("./Output/ITS_Permanova_Table.RDS") %>% 
  mutate(term=term %>% str_replace("drought","moisture"),amplicon = "ITS2") %>% 
  dplyr::select(amplicon,everything())

# Table 3 - PermANOVA
bind_rows(perm_fung,perm_bact) %>% 
  dplyr::filter(p.value < 0.05) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic()

# Separated permanova for each amplicon and host

fung_gf_perm <- readRDS("./Output/ITS_Permanova_Table_gf-only.RDS") %>% dplyr::mutate(host="Grand fir",amplicon="ITS2")
fung_sb_perm <- readRDS("./Output/ITS_Permanova_Table_sb-only.RDS") %>% dplyr::mutate(host="Snowbrush",amplicon="ITS2")
bact_gf_perm <- readRDS("./Output/16S_Permanova_Table_gf-only.RDS") %>% dplyr::mutate(host="Grand fir",amplicon="16S")
bact_sb_perm <- readRDS("./Output/16S_Permanova_Table_sb-only.RDS") %>% dplyr::mutate(host="Snowbrush",amplicon="16S")

permanova_results_df <- 
  rbind(fung_gf_perm,fung_sb_perm,bact_gf_perm,bact_sb_perm) %>% 
  dplyr::select(amplicon,host,everything()) %>% 
  dplyr::filter(!is.na(p.value)) %>% 
  dplyr::mutate(term = term %>% str_replace("drought","moisture"))

names(permanova_results_df) <- c("Amplicon","Host","Term","DF","SumOfSqs","R2","Statistic","P value" )
sig_rows <- which(permanova_results_df$`P value` < 0.05)

permanova_results_df %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic() %>% 
  kableExtra::row_spec(sig_rows,bold = TRUE)

permanova_results_df %>% 
  write_csv("./Output/permanova_table_for_each_host_and_amplicon.csv")
