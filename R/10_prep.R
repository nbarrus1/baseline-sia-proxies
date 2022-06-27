
# libraries
library(tidyverse)
library(here)


## Data --------------------

# bug field data
samplelog <- read_csv(here("data", "sample_log_bugs.csv"))

# taxon-ffg xref table
invert_group <- read_csv(here("data", "metadata_bugs.csv")) %>% 
  select(taxon, ffg) %>% 
  rename(taxon_code = taxon)

# C/N SI data
isodat <- read_csv(here("data", "sia_inverts_compiled.csv")) 

# longitudinal gradient proxy from PCA (Maitland 2020)
PCAdat <- read_csv(here("data", "data_PCA_results.csv")) %>% 
  # add 5 to gradient to aid interpretability
  mutate(PC1 = (PC1 + 5))

# land use, caluclated via GIS
LandUsedat <- read_csv(here("data", "land-use.csv"))

# add land use to site list / gradient
PCAdat <- PCAdat %>% 
  left_join(LandUsedat, by = "site_id")


## Explore reltionship between PC1 and natual land use ------------------
PCAdat %>% 
  ggplot(aes(x = PC1, y = lulc_nat*100))+
  geom_point(size = 3, shape = 21, color = "black", fill = "#666666")+
  geom_smooth(method = "lm", color = "black")+
  theme_classic()+
  labs(y = "Percent Natural Landcover",
       x = "PC1 (Environmental Gradient)")+
  geom_text(aes(label = site_id), vjust = 0, nudge_y = 1.5)+
  geom_text(aes(label = "r = - 0.700", x = 2.5,y = 25))+
  geom_text(aes(label = "p = 0.002", x = 2.5,y = 21))

cor.test(x = PCAdat$PC1, y = PCAdat$lulc_nat)
summary(lm(PCAdat$lulc_nat~PCAdat$PC1))


## Clean and combine data --------------------------

alldat <- isodat %>% 
  left_join(PCAdat, by = "site_id") %>% 
  # left_join(LandUsedat, by = "site_id") %>% 
  left_join(invert_group, by = "taxon_code") %>% 
  select(
    sample_year, stream_name, site_id,
    resource, taxon_code, d15N, PC1, ffg
    ) %>% 
  filter(site_id != "LR01") %>%  # remove test site
  filter(sample_year == "2016") %>%  # keep only 2016 data
  # keep only needed recourse groups
  filter(resource %in% c(
    "biofilm", "detritus" , "fish", "invert"))

# Subset data for analysis
invertdata <- alldat %>% 
  filter(resource == "invert") %>% 
  group_by(taxon_code) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5)

primprodat <- alldat %>% 
  filter(taxon_code == "Biofilm" | taxon_code == "Seston") %>% 
  group_by(taxon_code) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  group_by(taxon_code,site_id,PC1) %>% 
  summarise(mean_d15N = mean(d15N))






