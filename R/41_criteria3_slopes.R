# Criteria #3: models of baseline, basal resources, and ffg
# do the different potential baselines responded similarly to PC1 changes?

# Plot model outputs for all the groups


## Visualize --------------------

# by fishes
alldat |> 
  filter(taxon_code %in% c("Brown Trout", "Creek Chub", "White Sucker", "Longnose Dace", "Longnose Sucker")) |>
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
  filter(ffg != "Shredder")|>
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
  filter(taxon_code %in% c("Brown Trout", "Creek Chub", "White Sucker", "Longnose Dace", "Longnose Sucker")) |>
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

temp <- mods_ffg |>
  unnest(tidy)|>
  unnest(conf_int)|>
  filter(term == "PC1")

temp1 <- data.frame((unlist(temp$conf_int)))

get_cis <- function(data) {
  temp <- data |> 
    unnest(tidy)|>
    unnest(conf_int) |> 
    filter(term == "PC1")
  out <- data.frame((unlist(temp$conf_int)))
  print(out)
}

get_cis(mods_ffg)
get_cis(mods_baseline)

dfs <- list(mods_ffg, mods_taxa, mods_fish, mods_baseline)
cis <- map(dfs, get_cis)
str(cis)

mods_ffg |> 
  bind_cols(cis[[1]]) |> 
  unnest(tidy)|> 
  filter(term == "PC1") |> 
  ggplot(aes(x = fct_reorder(ffg, estimate, median), y = estimate)) +
  geom_point(size = 2, color = "black", shape = 21) + 
  geom_linerange(aes(ymin = X1, ymax = X2)) +
  coord_flip()+
  # scale_y_continuous(limits = c(0.00,1.02), breaks = seq(0.00,1.00,.25)) +
  labs(y = "", fill = "P < 0.05") + 
  theme_bw(base_size = 12)

mods_ffg |> 
  bind_cols(cis[[1]]) |> 
  unnest(tidy)|> 
  filter(term == "PC1") 

make_plot_df <- function(model, cis) {
  out <- model |> 
    bind_cols(cis) |> 
    unnest(tidy)|> 
    filter(term == "PC1") 
  out
}

make_plot_df(mods_ffg, cis[[1]])

plot_data <- map2(dfs, cis, make_plot_df)

plot_data[[1]] <- plot_data[[1]] |> rename(taxon_code = ffg)

plot_slopes <- function(df){
  df |>
    ggplot(aes(x = fct_reorder(taxon_code, estimate, median), y = estimate, fill = taxon_code)) +
    geom_linerange(aes(ymin = X1, ymax = X2)) +
    geom_point(size = 2, color = "black", shape = 21, show.legend = F) +
    geom_hline(yintercept = 0, color = "black", linewidth = .75, linetype = "dashed") +
    scale_fill_viridis_d()+
    coord_flip()+
    scale_y_continuous(limits = c(-.1,2), breaks = seq(0,2,.5))+
    labs(y = "") + 
    theme_bw(base_size = 12) + 
    theme(
      axis.ticks.length = unit(.25, "cm"),
      axis.title.y = element_text(vjust = 2), 
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

plot_slopes(plot_data[[2]])


p9 <- plot_data[[1]]|>
  filter(taxon_code != "Shredder")|>
  plot_slopes()+
  labs(x = "Feeding Group",
       y = NULL)
  
p10 <- plot_data[[2]]|>
  plot_slopes()+
  labs(x = "Taxonomic Group",
       y = NULL)

p11 <- plot_data[[3]]|>plot_slopes()+
  labs(x = "Fish Species",
       y = expression('slope ('~{beta}[1]~')'))

p12 <- plot_data[[4]]|>plot_slopes()+
  labs(x = "Basal Resource",
       y = NULL)



patch.crit3 <- p12+p6+p10+p8+p9+p7+p11+p5

patch.crit3.annote <- patch.crit3+
  plot_layout(ncol = 2,widths = c(1,2))& 
  plot_annotation(tag_levels = "A")

ggsave(here("out", "fig3_crit3.png"), 
       patch.crit3.annote, device = ragg::agg_png,
       units = "in", width = 8, height = 10)
