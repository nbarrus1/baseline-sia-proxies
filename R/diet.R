

diet <- read_csv(here("data", "diet.csv"))


diet |> 
  left_join(PCAdat, by = "site_id") |> select(-PC2) |> 
  ggplot(aes(PC1, IRI_100, color = food_group)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(vars(taxon_code)) + 
  ylim(c(-10,100))


