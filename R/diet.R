

diet <- read_csv(here("data", "diet.csv"))|> 
  left_join(PCAdat, by = "site_id") |> select(-PC2) %>% 
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
                                                                                              false = taxon_code))))),
         food_group = if_else(food_group == "amphib",
                              true = "amphibian",
                              false = if_else(food_group == "benthicInv",
                                              true = "benthic invertebrate",
                                              false = if_else(food_group == "terrInv",
                                                              true = "terrestrial invertebrate",
                                                              false = food_group))))
  


diet |> 
  filter(taxon_code == "White Sucker") %>% 
  ggplot(aes(y = IRI_100, x = fct_reorder(food_group, IRI_100), fill = food_group)) + 
  geom_boxplot(show.legend = F) + 
 # scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_y_continuous(limits = c(-7, 100), breaks = seq(0,100,25))+
  coord_flip()+
  facet_wrap(~taxon_code)+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Stomach Content", y = "IRI (%)",
       x = "Food Group")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

p33 <- diet |> 
  ggplot(aes(PC1, IRI_100, fill = food_group)) + 
  geom_smooth(aes(color = food_group), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_y_continuous(limits = c(-7, 100), breaks = seq(0,100,25))+
  facet_wrap(~taxon_code)+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Food Group", y = "IRI (%)",
       x = "PC1")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )+
  guides(fill = guide_legend(override.aes = list(size=5)))

ggsave(here("out", "fig4_stomachcontent.pdf"), 
       p33, device = pdf,
       units = "in", width = 9
       , height = 5)


diet_reg <- diet |> 
  group_by(taxon_code,food_group) |> 
  nest() |> 
  mutate(fit_pc1 = map(data, ~lm(IRI_100~PC1, data = .x)),
         tidy = map(fit_pc1, tidy),
         glance = map(fit_pc1, glance),
         augment = map(fit_pc1, augment))


 diet_reg |> unnest(glance)
 diet_reg |> unnest(tidy) |> mutate(p.value = round(p.value, digits = 3))
