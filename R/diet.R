

diet <- read_csv(here("data", "diet.csv"))|> 
  left_join(PCAdat, by = "site_id") |> select(-PC2) 


diet |> 
  filter(taxon_code == "WHS") %>% 
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
  labs(fill = "Stomach Content", y = "TP",
       x = "PC1")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

ggsave(here("out", "fig5_stomachcontent.png"), 
       p33, device = ragg::agg_png,
       units = "in", width = 8, height = 5)


diet_ancova <- diet |>
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(IRI_100 ~ PC1*food_group, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2)))

geom_boxplot()

temp2 <- diet_ancova |> unnest(glance)
temp <- diet_ancova |> unnest(tidy)
