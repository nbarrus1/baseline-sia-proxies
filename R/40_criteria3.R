
# Criteria #3: models of baseline, basal resources, and ffg
# do the different potential baselines responded similarly to PC1 changes?

# Set factor names
fish_names <- c(
  "uncorrected",
  levels(as.factor(isodat_fish$taxon_code))
  )
fish_names_common <- c(
  "uncorrected",
  levels(as.factor(isodat_fish$common_name))
)
taxon_names <- c(
  "uncorrected",
  levels(as.factor(isodat_baseline$taxon_code)),
  levels(as.factor(isodat_bug_common$taxon_code))
  )

# set palettes
palette_fish <- setNames(object = my_pallette[1:length(fish_names)], nm = fish_names)
palette_fish_common <- setNames(object = my_pallette[1:length(fish_names_common)], nm = fish_names_common)
palette_taxon <-setNames(object = my_pallette[1:length(taxon_names)], nm = taxon_names)


## Visualize --------------------

# by fishes
p5 <- isodat_fish |> 
  ggplot(aes(PC1, d15N, fill = common_name)) + 
  geom_smooth(aes(color = common_name), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_fill_manual(values=palette_fish_common) + 
  scale_color_manual(values=palette_fish_common) + 
  labs(fill = "Fish Species", y = expression(' '~{delta}^15*N~' '),
       x = "Longitudinal Gradient (PC1)")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(linewidth = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# resources
p6 <- isodat_baseline |> 
  filter(taxon_code %in% c("Biofilm","FBOM", "filamentous", "Seston")) |>
  ggplot(aes(PC1, d15N, fill = taxon_code)) + 
  geom_smooth(
    data = isodat_baseline |> filter(taxon_code %in% c("Biofilm","FBOM", "Seston")),
    aes(color = taxon_code), method = "lm", se = FALSE, show.legend = FALSE
    )+
  geom_point(size = 2, color = "black", shape = 21, show.legend = FALSE) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_fill_manual(values=palette_taxon) + 
  scale_color_manual(values=palette_taxon) + 
  labs(fill = "Basal Resource", y = expression(' '~{delta}^15*N~' '),
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(linewidth = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# inverts by ffg
p7 <- isodat_bug_common |> 
  ggplot(aes(PC1, d15N, fill = ffg)) + 
  geom_smooth(aes(color = ffg), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_fill_manual(values=palette_ffgs) + 
  scale_color_manual(values=palette_ffgs) + 
  labs(fill = "Feeding Group", y = expression(' '~{delta}^15*N~' '),
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(linewidth = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# inverts by taxa
p8 <- isodat_bug_common |> 
  ggplot(aes(PC1, d15N, fill = taxon_code)) + 
  geom_smooth(aes(color = taxon_code), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_fill_manual(values=palette_taxon) + 
  scale_color_manual(values=palette_taxon) + 
  labs(fill = "Taxonomic Group", y = expression(' '~{delta}^15*N~' '),
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(linewidth = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )

## Fit models --------------------------

# ffgs
mods_ffg <- isodat_bug_common |> 
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
mods_taxa <- isodat_bug_common |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

mods_fish <- isodat_fish |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 'PC1'))
  ) 

mods_baseline <- isodat_baseline |> 
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

# plot_sigcorr <- function(df, x_cat, alpha = 0.05) {
#   df |>
#     unnest(glance)|>
#     mutate(sig = if_else(p.value < alpha, true = "y", false = "n")) |>
#     ggplot(aes(x=fct_reorder({{ x_cat }}, adj.r.squared, median), y = adj.r.squared, fill = sig))+         #adj r sqr vs reorderd taxon   
#     geom_segment(aes(xend = {{x_cat}}), yend = 0, colour = "grey50") +
#     geom_point(size = 2, color = "black", shape = 21) + 
#     coord_flip()+
#     scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25)) +
#     scale_fill_manual(values = c("midnightblue","darkred"))+
#     labs(y = "", fill = "p < 0.05") + 
#     theme_bw(base_size = 12) + 
#     theme(
#       # axis.line = element_line(linewidth = .5),
#       axis.ticks.length = unit(.25, "cm"),
#       axis.title.y = element_text(vjust = 2), 
#       panel.grid.minor.x = element_blank(),
#       panel.grid.minor.y = element_blank(),
#       panel.grid.major.y = element_blank(), 
#       plot.margin = unit(c(0,0,0,0), "cm")
#     )
# }

#p5 <- mods_ffg|>plot_sigcorr(x_cat = ffg)+labs(x = "Feeding Group")
#p6 <- mods_taxa|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Taxonomic Group")
#p7 <- mods_fish|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Fish Species")
#p8 <- mods_baseline|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Basal Resource")


# temp <- mods_fish %>% unnest(glance)
