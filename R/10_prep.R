# Prepare data sets for analyses


# libraries
library(tidyverse)
library(here)
library(patchwork)
library(broom)


# functions
# lower_ci <- function(mean, se, n, conf_level = 0.95){
#   lower_ci <- mean - qt(1-((1-conf_level)/2),n-1)*se
# }
# upper_ci <- function(mean, se, n, conf_level = 0.95){
#   upper_ci <- mean + qt(1-((1-conf_level)/2),n-1)*se
# }

## Data --------------------

# bug field data
samplelog <- read_csv(here("data", "sample_log_bugs.csv"))

# fuzzy coded taxon-ffg xref table
invert_group <- read_csv(here("data", "metadata_bugs.csv")) %>% 
  select(taxon, ffg) %>% 
  rename(taxon_code = taxon)

# C/N SI data from Maitland (2020)
isodat <- read_csv(here("data", "sia_inverts_compiled.csv")) %>% 
  # rename fish to common names
  mutate(taxon_code = if_else(taxon_code == "BNT",
                              true = "Brown Trout",
                              false = if_else(taxon_code == "CKC",
                                              true = "Creek Chub", 
                                              false = if_else(taxon_code == "LND",
                                                              true = "Longnose Dace",
                                                              false = if_else(taxon_code == "LNS",
                                                                              true = "Longnose Sucker",
                                                                              false = if_else(taxon_code == "WHS",
                                                                                              true = "White Sucker",
                                                                                              false = taxon_code)
                                                                              )))))

samplelog |> 
  filter(sample_year==2016) |> 
  left_join(isodat |> filter(compartment == "invert"), by = "unique_id", relationship = "many-to-many") |> 
  View()


# longitudinal gradient proxy from PCA (Maitland 2020)
PCAdat <- read_csv(here("data", "data_PCA_results.csv")) %>% 
  # add 5 to gradient to aid interpretability
  mutate(PC1 = (PC1 + 5))

# # land use, calculated via GIS
# LandUsedat <- read_csv(here("data", "land-use.csv"))
# 
# # add land use to site list / gradient
# PCAdat <- PCAdat %>% 
#   left_join(LandUsedat, by = "site_id")


## Explore relationship between PC1 and natual land use ------------------

# PCAdat %>% 
#   ggplot(aes(x = PC1, y = lulc_nat*100))+
#   geom_point(size = 3, shape = 21, color = "black", fill = "#666666")+
#   geom_smooth(method = "lm", color = "black")+
#   theme_classic()+
#   labs(y = "Percent Natural Landcover",
#        x = "PC1 (Environmental Gradient)")+
#   geom_text(aes(label = site_id), vjust = 0, nudge_y = 1.5)+
#   geom_text(aes(label = "r = - 0.700", x = 2.5,y = 25))+
#   geom_text(aes(label = "p = 0.002", x = 2.5,y = 21))
# 
# cor.test(x = PCAdat$PC1, y = PCAdat$lulc_nat)
# summary(lm(PCAdat$lulc_nat~PCAdat$PC1))


## Clean and combine data --------------------------

alldat <- isodat %>% 
  # add PC1 values
  left_join(PCAdat, by = "site_id") %>% 
  # left_join(LandUsedat, by = "site_id") %>% 
  # add ffgs to taxa
  left_join(invert_group, by = "taxon_code") %>% 
  select(
    unique_id, sample_year, hitch = sample_hitch, stream_name, site_id, compartment,
    resource, taxon_code, length_mm, d15N, d13C, CN, PC1, ffg
    ) %>% 
  filter(sample_year == "2016") %>%  # keep only 2016 data (only year bug data run)
  filter(resource %in% c("biofilm", "detritus" , "fish", "invert", "photo"))

# Subset data for analysis
invertdata <- alldat %>% 
  select(-length_mm) |> 
  filter(compartment == "invert") 
  # group_by(taxon_code) %>% 
  # mutate(n = n()) %>% 
  # ungroup() %>% 
  # filter(n >= 5)

# Why do we filter for n > 5 here after filtering our the inverts?
# This effectively removes 26 samples: 
# Corbiculidae, Odontoceridae, Orconectes, Ceratopogonidae, 
# Amphipoda, Hydrophilidae, Oligoneuriidae


primprodat <- alldat %>%   # bad object name
  select(-length_mm) |> 
  filter(taxon_code == "Biofilm" | taxon_code == "Seston") %>% 
  group_by(taxon_code, site_id, PC1) %>% 
  summarise(mean_d15N = mean(d15N), .groups = "drop")


# # fish species
# isodat |> filter(compartment == "fish") |> arrange(site_id, taxon_code) |> 
#   write_csv(here("data","fish.csv"))


# Samples per taxa by site and hitch
invertdata |> 
  group_by(site_id, hitch, taxon_code) |> 
  count() |> 
  pivot_wider(id_cols = c(site_id, hitch), names_from = taxon_code, values_from = n) |> 
  arrange(site_id, hitch) |> 
  View()

# I think the reason we are missing data from those hitches 

