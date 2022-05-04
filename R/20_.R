
# Criteria #1: distribution

## Prep

library(tidyverse) 
library(here)
library(cowplot)


## Data --------------

# obtain taxon to filter the large dataset
focustax <- as.factor(invertdata$taxon_code) %>% levels()
# obtain Sites id from data
study_sites <- as.factor(invertdata$site_id) %>% levels()
# total number of sites
n_sites = length(study_sites)
# obtain taxon to filter the large dataset
focusffg <- as.factor(invertdata$ffg) %>% levels()


## Calculate distribution

# add all values of presence by the group taxon divide total number 
# of sites multiply by 100 to get percentage

# by taxa
Distr_dat <- invertdata %>%
  group_by(taxon_code, site_id) %>% 
  summarise(n = n()) %>%
  mutate(presence = 1) %>%
  group_by(taxon_code) %>%
  summarise(distr = (sum(presence)/n_sites)*100) %>%
  filter(taxon_code %in% focustax)

# by ffg
Distr_dat_FFG <- invertdata %>%
  filter(taxon_code %in% focustax) %>% 
  group_by(ffg, site_id) %>%
  summarise(n = n()) %>%
  mutate(presence = 1) %>%
  group_by(ffg) %>%
  summarise(distr = (sum(presence)/n_sites)*100)


## Plot

DistFigTaxon <- Distr_dat %>%
  filter(taxon_code %in% focustax) %>%
  ggplot(aes(x=fct_reorder(taxon_code, distr, median), y=distr)) +
  geom_point(size = 3,color = "black", shape = 21, fill = "#666666") + 
  #coord_flip()+
  geom_hline(yintercept = 75, color = "red", size =1, linetype = "dashed") +
  scale_y_continuous(limits = c(0,100)) +
  theme_classic() +
  labs(x = "Taxon", y = "Percent of Sites") +
  theme(
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    title = element_text(size = 16, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = rel(0.8)),
    axis.text.x.bottom = element_text(size = 10, face = "bold", color = "black")
    )

DistFigFFG <- Distr_dat_FFG %>% 
  ggplot(aes(x = fct_reorder(ffg, distr, median), y = distr))+
  geom_hline(yintercept = 75, color = "red", linetype = "dashed", size = 1)+
  geom_point(size = 3, shape = 21, color = "black",fill = "#666666")+
  #coord_flip()+
  scale_y_continuous(limits = c(0,100)) +
  labs(x = "FFG", y = "Percent of Sites")+
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, face = "bold", color = "black"),
    title = element_text(size = 16, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 90)
    )

# panel plot
distFigAll <- cowplot::plot_grid(
  DistFigTaxon, DistFigFFG, labels = c("A","B"), align = "hv",
  label_size = 18, ncol = 2, nrow = 1, label_x = 0.15, label_y = 0.99
  )

# save 
ggsave(here("out", "combined_distribution.png"), 
       distFigAll, device = ragg::agg_png,
       units = "in", width = 8.5, height = 6)
