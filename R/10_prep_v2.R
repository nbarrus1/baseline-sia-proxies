
# libraries
library(tidyverse)
library(here)
library(patchwork)
library(broom)
library(MetBrewer)
library(ggpubr)

## Data

# longitudinal gradient variable from PCA (Maitland 2020)
PCAdat <- read_csv(here("data", "data_PCA_results.csv")) %>% 
  mutate(PC1 = (PC1 + 5))  # add 5 to gradient to aid interpretability

# fuzzy coded taxon-ffg xref table
invert_group <- read_csv(here("data", "metadata_bugs.csv")) %>% 
  select(taxon, ffg) %>% 
  rename(taxon_code = taxon)

# bug field data (all records for observed species)
samplelog <- read_csv(here("data", "sample_log_bugs.csv")) |> 
  filter(sample_year==2016)  |> # Keep 2016 data only
  filter(site_id != "LR01") |>  # remove scouting site
  filter(gear_type == "dnet")  # keep dnet data only (no hess)

# Stable isotope data for inverts, fish, and baselines
isodat <- read_csv(here("data", "maitland_2020_SI_data.csv")) |> 
  filter(site_id != "LR01") |>  # remove scouting site
  filter(sample_year==2016) # Keep 2016 data only


# Subset tables for main groups
isodat_bug <- isodat |> filter(resource %in% c("invert"))
isodat_fish <-  isodat |> filter(resource %in% c("fish")) |> filter(taxon_code %in% c("BNT","CKC","LND","LNS","WHS"))
isodat_baseline <-  isodat |> filter(taxon_code %in% c("Biofilm","FBOM", "filimentous", "Seston"))

# Merge bug isotope data to sample log, PC1, and ffgs
isodat_bug <- samplelog |> 
  select(-gear_type, -coding_id, -length_mm, -sex, -date_sorted, -taxonomist, -bug_sample_status, -bug_samp_location, -notes) |> 
  left_join(
    isodat_bug |> select(unique_id, compartment, taxon_code, d13C, d15N, CN), 
    by = join_by(unique_id, taxon_code), relationship = "many-to-many") |> 
  left_join(invert_group, by = "taxon_code") |> 
  left_join(PCAdat, by = "site_id") |> 
  filter(taxon_code != "UNKN") |> 
  filter(! str_detect(ffg, "terrestrial")) |>
  mutate(compartment = if_else(is.na(compartment), "invert", compartment)) |> 
  arrange(site_id, sample_hitch, taxon_code)

# Merge PC1 values to fish and baseline tables
isodat_fish <- isodat_fish |> left_join(PCAdat, by = "site_id") 
isodat_baseline <- isodat_baseline |> left_join(PCAdat, by = "site_id") 

# Save prepared data for analysis
# save(isodat_bug, isodat_fish, isodat_baseline, file = here("data","data_analysis.Rdata"))

