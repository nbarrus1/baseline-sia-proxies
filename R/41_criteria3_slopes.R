# Criteria #3: models of baseline, basal resources, and ffg
# do the different potential baselines responded similarly to PC1 changes?

# Plot model outputs for all the groups


## Visualize --------------------

# by fishes
alldat |> 
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |>
  # filter(compartment=="fish") |> 
  ggplot(aes(PC1, d15N, color = taxon_code)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

# resources
alldat |> 
  filter(taxon_code %in% c("Biofilm","FBOM", "filimentous", "Seston")) |> 
  ggplot(aes(PC1, d15N, color = taxon_code)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

# inverts by ffg
invertdata_sub |> 
  ggplot(aes(PC1, d15N, color = ffg)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

# inverts by taxa
invertdata_sub |> 
  ggplot(aes(PC1, d15N, color = taxon_code)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


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
mods_taxa <- invertdata |> 
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

temp <- mods_ffg |> unnest(conf_int)

temp$conf_int
## Plot ---------

#plot f(x)n

plot_slopes <- function(df, x_cat){
  upp <- df |>
    unnest(tidy)|>
    unnest(conf_int)|>
    filter(term == "PC1")|>
    ggplot(aes(x = fct_reorder({{ x_cat }}, estimate, median), y = estimate)) +
    geom_linerange(ymin = conf_int, ymax = upp) +
    geom_point(size = 2, color = "black", shape = 21) + 
    coord_flip()+
    scale_y_continuous(limits = c(0.00,1.02), breaks = seq(0.00,1.00,.25)) +
    labs(y = "", fill = "P < 0.05") + 
    theme_bw(base_size = 12) + 
    theme(
      axis.line = element_line(linewidth = .5),
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

