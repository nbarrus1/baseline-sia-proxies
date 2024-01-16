
diet <- read_csv(here("data", "diet.csv"))|> 
  left_join(PCAdat, by = "site_id") |> select(-PC2) %>% 
  mutate(taxon_code = case_when(taxon_code == "BNT" ~ "brown trout",
                                taxon_code == "CKC" ~ "creek chub",
                                taxon_code == "LND" ~ "longnose dace",
                                taxon_code == "LNS" ~ "longnose sucker",
                                taxon_code == "WHS" ~ "white sucker"),
         food_group = case_when(food_group == "amphib" ~ "amphibian",
                                food_group == "benthicInv" ~ "benthic invertebrate",
                                food_group == "terrInv" ~ "terrestrial invertebrate",
                                .default = as.character(food_group)))

diet_names <- levels(as.factor(diet$food_group))
palette_diet <- setNames(object = my_pallette[1:length(diet_names)], nm = diet_names) 


diet |> 
  filter(taxon_code == "white sucker") |>  
  ggplot(aes(y = IRI_100, x = fct_reorder(food_group, IRI_100), fill = food_group)) + 
  geom_boxplot(show.legend = F) + 
 # scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_y_continuous(limits = c(-7, 100), breaks = seq(0,100,25))+
  coord_flip()+
  facet_wrap(~taxon_code)+
  scale_color_manual(values=palette_diet) + 
  scale_fill_manual(values=palette_diet) + 
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
  scale_color_manual(values=palette_diet) + 
  scale_fill_manual(values=palette_diet) + 
  labs(fill = "Food Group", y = "IRI (%)",
       x = "PC1")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )+
  guides(fill = guide_legend(override.aes = list(size = 4)))
  

path <- here::here("out", "criteria-4-diet")
ggsave(glue::glue("{path}.pdf"), plot = p33, 
       width = 4, height = 2, scale = 2.5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)

#ggsave(here("out", "fig4_stomachcontent.png"), 
#       p33, device = ragg::agg_png,
#       units = "in", width = 9
#       , height = 5)


mods_diet <- diet |>
  group_by(taxon_code, food_group) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(IRI_100 ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2)))

mods_diet |> unnest(glance)

