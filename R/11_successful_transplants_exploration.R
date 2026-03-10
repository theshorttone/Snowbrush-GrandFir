library(tidyverse)
library(ranger)
library(vip)

# inoculum data by site
meta <- read_csv("./Data/Sample_Metadata_Full.csv") %>% 
  dplyr::filter(!is.na(inoculum_site) & inoculum_site != "Sterile")
meta$site <- paste0("Site",meta$inoculum_site)


# read in successful taxa
fung <- read_csv("./Output/ITS_Successful_dataframe.csv") %>% 
  unique.data.frame() %>% 
  mutate(amplicon = "ITS")
bact <- read_csv("./Output/16S_Successful_dataframe.csv") %>% 
  unique.data.frame() %>% 
  mutate(amplicon = "16S")
amf <- read_csv("./Output/18S_Successful_dataframe.csv") %>% 
  unique.data.frame() %>% 
  mutate(amplicon = "18S")

# merge them and add "success" column FALSE
success <- fung %>% 
  rbind(bact) %>% 
  rbind(amf) %>% 
  mutate(success = TRUE)

# read in unsuccessful taxa
ufung <- read_csv("./Output/ITS_Unsuccessful_dataframe.csv") %>% 
  unique.data.frame() %>% 
  mutate(amplicon = "ITS")
ubact <- read_csv("./Output/16S_Unsuccessful_dataframe.csv") %>% 
  unique.data.frame() %>% 
  mutate(amplicon = "16S")
uamf <- read_csv("./Output/18S_Unsuccessful_dataframe.csv") %>% 
  unique.data.frame() %>% 
  mutate(amplicon = "18S")
# merge them and add "success" column FASLE
failure <- 
  ufung %>% 
  rbind(ubact) %>% 
  rbind(uamf) %>% 
  mutate(success = FALSE)

# merge both
dat <- full_join(success,failure)

# add inoculum metadata
dat <- meta %>% 
  dplyr::select(site,fire_freq) %>% 
  unique.data.frame() %>% 
  full_join(dat)


table(dat$success)

mod <- 
  dat %>% 
  dplyr::filter(amplicon == "18S") %>% 
  dplyr::select(success,fire_freq,Genus) %>% 
  glm(data=.,
    formula = success ~ .,
    family='binomial')
summary(mod)
dat$pred <- predict(mod,dat,type='response')

dat

dat %>% 
  dplyr::filter(amplicon == "18S") %>% 
  group_by(Kingdom,Phylum,Class,Order,Family,Genus,Species,amplicon) %>% 
  summarize(pred=mean(pred,na.rm=TRUE)) %>% 
  ggplot(aes(x=Family,y=pred)) +
  geom_boxplot() +
  facet_wrap(~Phylum,scales='free')


# rf models
rf <- ranger(data=dat[complete.cases(dat),],
             formula = success ~ .,
             importance = 'permutation')
vip(rf)
pred <- predict(rf,dat[complete.cases(dat),],type = 'response')
dat2 <- dat[complete.cases(dat),]
dat2$pred <- as.logical(pred$predictions)

dat2$agree <- ifelse(dat2$success == dat2$pred,TRUE,FALSE)
family_success_18S <- 
dat2 %>% 
  filter(amplicon == "18S") %>%
  group_by(Family,Genus) %>% 
  summarize(perc_success = sum(success)/n()) %>% 
  arrange(desc(perc_success))

family_success_18S$Genus <- family_success_18S$Genus %>% factor(levels = family_success_18S$Genus)  
family_success_18S %>% 
ggplot(aes(y=perc_success,x=Genus,fill=Family)) +
  geom_col(stat='count') +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  scale_fill_viridis_d(option = 'turbo') +
  labs(y="Transplat success rate")
ggsave("Output/figs/18S_Success_Rate_Genus.png",width = 6,height = 5)


phylum_success_16S <- 
  dat2 %>% 
  filter(amplicon == "16S") %>%
  group_by(Phylum,Class) %>% 
  summarize(perc_success = sum(success)/n()) %>% 
  arrange(desc(perc_success))

phylum_success_16S$Genus <- phylum_success_16S$Class %>% factor(levels = phylum_success_16S$Class %>% unique)  
phylum_success_16S %>% 
  ggplot(aes(y=perc_success,x=Class,fill=Phylum)) +
  geom_col(stat='count') +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  scale_fill_viridis_d(option = 'turbo') +
  labs(y="Transplat success rate")

phylum_success_ITS <- 
  dat2 %>% 
  filter(amplicon == "ITS") %>%
  group_by(Phylum,Class) %>% 
  summarize(perc_success = sum(success)/n()) %>% 
  arrange(desc(perc_success))

phylum_success_ITS$Genus <- phylum_success_ITS$Class %>% factor(levels = phylum_success_ITS$Class %>% unique)  

phylum_success_ITS %>% 
  ggplot(aes(y=perc_success,x=Class,fill=Phylum)) +
  geom_col(stat='count') +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  scale_fill_viridis_d(option = 'turbo') +
  labs(y="Transplat success rate")

