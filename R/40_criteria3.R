
# Criteria #3: models of baseline, basal resources, and ffg
# do the different potential baselines responded similarly to PC1 changes?

# Plot model outputs for all the groups


## Visualize --------------------

# by fishes
p5 <- alldat |> 
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |>
  # filter(compartment=="fish") |> 
  ggplot(aes(PC1, d15N, fill = taxon_code)) + 
  geom_smooth(aes(color = taxon_code), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Fish Species", y = expression(' '~{delta}^15*N~' '),
       x = "Longitudinal Gradient (PC1)")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# resources
p6 <- alldat |> 
  filter(taxon_code %in% c("Biofilm","FBOM", "filimentous", "Seston")) |> 
  ggplot(aes(PC1, d15N, fill = taxon_code)) + 
  geom_smooth(aes(color = taxon_code), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Basal Resource", y = expression(' '~{delta}^15*N~' '),
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# inverts by ffg
p7 <- invertdata |> 
  filter(ffg != "Shredder") |>
  ggplot(aes(PC1, d15N, fill = ffg)) + 
  geom_smooth(aes(color = ffg), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Feeding Group", y = expression(' '~{delta}^15*N~' '),
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# inverts by taxa
p8 <- invertdata_sub |> 
  ggplot(aes(PC1, d15N, fill = taxon_code)) + 
  geom_smooth(aes(color = taxon_code), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Taxonomic Group", y = expression(' '~{delta}^15*N~' '),
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

## Fit models --------------------------

# ffgs
mods_ffg <- invertdata |> 
  group_by(ffg) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

# taxa
mods_taxa <- invertdata_sub |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

mods_fish <- alldat |> 
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |>
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 'PC1'))
  ) 

mods_baseline <- alldat |> 
  filter(taxon_code %in% c("Biofilm","FBOM", "filimentous", "Seston")) |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 


mods_ffg |> unnest(glance)


## Plot ---------


plot_sigcorr <- function(df, x_cat, alpha = 0.05) {
  df |>
    unnest(glance)|>
    mutate(sig = if_else(p.value < alpha, true = "y", false = "n")) |>
    ggplot(aes(x=fct_reorder({{ x_cat }}, adj.r.squared, median), y = adj.r.squared, fill = sig))+         #adj r sqr vs reorderd taxon   
    geom_segment(aes(xend = {{x_cat}}), yend = 0, colour = "grey50") +
    geom_point(size = 2, color = "black", shape = 21) + 
    coord_flip()+
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25)) +
    scale_fill_manual(values = c("midnightblue","darkred"))+
    labs(y = "", fill = "p < 0.05") + 
    theme_bw(base_size = 12) + 
    theme(
      # axis.line = element_line(linewidth = .5),
      axis.ticks.length = unit(.25, "cm"),
      axis.title.y = element_text(vjust = 2), 
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

p5 <- mods_ffg|>plot_sigcorr(x_cat = ffg)+labs(x = "Feeding Group")
p6 <- mods_taxa|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Taxonomic Group")
p7 <- mods_fish|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Fish Species")
p8 <- mods_baseline|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Basal Resource")


temp <- mods_fish %>% unnest(glance)
