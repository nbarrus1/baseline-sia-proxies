
# Criteria #1: distribution

## Data --------------

# obtain taxon to filter the large dataset
focustax <- as.factor(invertdata$taxon_code) %>% levels()
# obtain Sites id from data
study_sites <- as.factor(invertdata$site_id) %>% levels()
# total number of sites
n_sites = length(study_sites)
# obtain taxon to filter the large dataset
focusffg <- as.factor(invertdata$ffg) %>% levels()


## Calculate distribution

# add all values of presence by the group taxon divide total number 
# of sites multiply by 100 to get percentage

# by taxa
Distr_dat <- invertdata |> 
  distinct(site_id, taxon_code) |> 
  mutate(presence = 1) |> 
  group_by(taxon_code) |> 
  summarise(distr = (sum(presence)/n_sites)*100) %>%
  left_join(invert_group, by = "taxon_code") |> 
  rename(group = ffg)

# by ffg
Distr_dat_FFG <- invertdata %>%
  # filter(taxon_code %in% focustax) %>% 
  distinct(ffg, site_id) %>%
  mutate(presence = 1) %>%
  group_by(ffg) %>%
  summarise(distr = (sum(presence)/n_sites)*100)|> 
  mutate(group = ffg)


# plot function
plot_dist <- function(df, x_cat){
  df |>
    ggplot(aes(x = fct_reorder({{ x_cat }}, distr, median), y = distr, fill = group)) +
    geom_segment(aes(xend = {{x_cat}}), yend = 0, colour = "grey50") +
    geom_point(size = 2, color = "black", shape = 21) + 
    geom_hline(yintercept = 75, color = "black", size = .75, linetype = "dashed") +
    coord_flip()+
    scale_y_continuous(limits = c(0,102), breaks = seq(0,100,25)) +
    labs(y = "", fill = "Feeding Group") + 
    theme_bw(base_size = 12) + 
    theme(
      axis.line = element_line(size = .5),
      axis.ticks.length = unit(.25, "cm"),
      axis.title.y = element_text(vjust = 2), 
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

## Plot
p1 <- Distr_dat |> plot_dist(taxon_code) + labs(x = "Taxonomic Group")
p2 <- Distr_dat_FFG |> plot_dist(ffg) + labs(x = "Feeding Group", y = "Sites present (%)")
patch.dist <- p1 / p2
patch.dist.annote <- patch.dist + 
  plot_layout(heights = c(3, 1), guides = 'collect') & 
  plot_annotation(tag_levels = "A")
patch.dist.annote

# save plot
ggsave(here("out", "crit1_distributions.png"), 
       patch.dist.annote, device = ragg::agg_png,
       units = "in", width = 6, height = 6)
