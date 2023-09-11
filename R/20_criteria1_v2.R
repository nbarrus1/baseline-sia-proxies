
# Criteria #1: distribution

## Data --------------

# load(here("data","data_analysis.Rdata"))

# # obtain taxon to filter the large dataset
# focustax <- as.factor(invertdata$taxon_code) %>% levels()
# obtain taxon to filter the large dataset
# focusffg <- as.factor(isodat_bug$ffg) %>% levels()

# # obtain Sites id from data
study_sites <- as.factor(isodat_bug$site_id) %>% levels()
# total number of sites
n_sites = length(study_sites)


## Calculate occupancy (distributions) ------------------------------------


# by taxa
Distr_dat <- isodat_bug |> 
  distinct(site_id, taxon_code) |> 
  mutate(presence = 1) |> 
  group_by(taxon_code) |> 
  summarise(
    total_pres = sum(presence), 
    distr = (total_pres/n_sites)*100) %>%
  left_join(invert_group, by = "taxon_code") |> 
  rename(group = ffg)

# by ffg
Distr_dat_FFG <- isodat_bug %>%
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
    geom_hline(yintercept = 75, color = "black", linewidth = .75, linetype = "dashed") +
    coord_flip()+
    scale_y_continuous(limits = c(0,102), breaks = seq(0,102,25), expand = c(0,0)) +
    scale_fill_manual(values=met.brewer("Lakota", 6)) + 
    labs(y = "", fill = "Feeding Group") +
    theme_bw(8) +
    theme(
      axis.line = element_line(linewidth = .5),
      axis.title.y = element_text(vjust = 2),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

## Plot
p1 <- Distr_dat |> plot_dist(taxon_code) + labs(x = "Taxonomic Group", y = "Sites present (%)")
p2 <- Distr_dat_FFG |> plot_dist(ffg) + labs(x = "Feeding Group", y = "Sites present (%)")
patch.dist <- p1 / p2
patch.dist.annote <- patch.dist + 
  plot_layout(heights = c(5, 1)) + 
  plot_layout(guides = 'collect') & 
  plot_annotation(tag_levels = "A")
patch.dist.annote


# Save plot
path <- here::here("out", "criteria-1")
ggsave(glue::glue("{path}.pdf"), plot = last_plot(), 
       width = 3, height = 4, scale = 2, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)

