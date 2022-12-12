
# Criteria #1: distribution

## Prep

# library(tidyverse) 
# library(here)
# library(patchwork)


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
  summarise(distr = (sum(presence)/n_sites)*100) 

# by ffg
Distr_dat_FFG <- invertdata %>%
  # filter(taxon_code %in% focustax) %>% 
  distinct(ffg, site_id) %>%
  mutate(presence = 1) %>%
  group_by(ffg) %>%
  summarise(distr = (sum(presence)/n_sites)*100)


# plot function
plot_dist <- function(df, x_cat){
  df |>
    ggplot(aes(x = fct_reorder({{ x_cat }}, distr, median), y = distr)) +
    geom_segment(aes(xend = {{x_cat}}), yend = 0, colour = "grey50") +
    geom_point(size = 3,color = "black", shape = 21, fill = "#666666") + 
    geom_hline(yintercept = 75, color = "red", size =1, linetype = "dashed") +
    scale_y_continuous(limits = c(0,100)) +
    ylab("Sites present (%)") + 
    theme_bw(base_size = 12) + 
    theme(
      axis.line = element_line(size = .5),
      axis.ticks.length = unit(.25, "cm"),
      axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
      panel.grid.major.x = element_blank(),   # No vert grid lines
      panel.grid.minor.y = element_blank(),   # No vert grid lines
    )
}



## Plot

p1 <- Distr_dat |> plot_dist(taxon_code) + xlab("Taxonomic Group")
p2 <- Distr_dat_FFG |> plot_dist(ffg) + xlab("Functional Feeding Group")

patch <- p1 | p2

patch_ann <- patch + 
  plot_layout(widths = c(2, 1)) & 
  plot_annotation(tag_levels = "A")

# save 
ggsave(here("out", "crit1_combined_distribution.png"), 
       patch_ann, device = ragg::agg_png,
       units = "in", width = 8.5, height = 4)
