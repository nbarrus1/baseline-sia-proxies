
# Criteria 4: Trophic Estimates independent of Land Use

## Prep

###uncorrected

#calculate TP
TP_uncorrected <- alldat |> 
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |>
  mutate(discFactor = ((-0.281*d15N)+5.879),
         TP = d15N/discFactor,
         baseline = "uncorrected")

#plot
TP_uncorrected |> 
  ggplot(aes(x = PC1, y = TP, color = taxon_code))+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


###Corrected TP ####

###Taxa

#calculate mean d15N for each baseline at each site
taxa_correction <- invertdata_sub |>
  group_by(site_id,taxon_code)|>
  summarise(correction = mean(d15N, na.rm = T),
            correction_type = "taxa") %>% 
  rename(baseline = taxon_code)

#use mean d15N of the baseline (taxa) at each site to calculate corrected 
#Trophic positions

TP_taxa <- alldat |>
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |> 
  left_join(taxa_correction, by = "site_id")|>
  mutate(discFactor = ((-0.281*d15N)+5.879),
         TP = (((d15N-correction)/discFactor)+2)) 

#plot the TPs for each fish spp.

TP_taxa |> 
  ggplot(aes(x = PC1, y = TP, color = baseline))+
  facet_wrap(~taxon_code)+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

###ffg

#calculate mean d15N for each baseline at each site
ffg_correction <- invertdata |>
  filter(ffg != "Shredder")|>
  group_by(site_id,ffg)|>
  summarise(correction = mean(d15N, na.rm = T),
            correction_type = "ffg") %>% 
  rename(baseline = ffg)

#use mean d15N of the baseline (taxa) at each site to calculate corrected 
#Trophic positions
TP_ffg <- alldat |>
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |> 
  left_join(ffg_correction, by = "site_id")|>
  mutate(discFactor = ((-0.281*d15N)+5.879),
         TP = (((d15N-correction)/discFactor)+2))


TP_ffg |> 
  ggplot(aes(x = PC1, y = TP, color = baseline))+
  facet_wrap(~taxon_code)+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 



#####fit models-------------

# uncorrected
mods_uncorrected <- TP_uncorrected |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(TP ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2)),
    baseline = "uncorrected"
  ) 

##TPs corrected by baselines (ffgs)

mods_corrected_ffg <- TP_ffg |> 
  group_by(taxon_code,baseline) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(TP ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) %>% 
  bind_rows(mods_uncorrected)
##TPs corrected by baselines (taxa)

mods_corrected_taxa <- TP_taxa |> 
  group_by(taxon_code,baseline) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(TP ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) %>% 
  bind_rows(mods_uncorrected)

temp <- mods_corrected_ffg |> unnest(glance)
mods_corrected_taxa|> unnest(glance)

###slope plots####

dfs_crit4 <- list(mods_corrected_ffg,mods_corrected_taxa)
cis_crit4 <- map(dfs_crit4, get_cis)
plot_data_crit4 <- map2(dfs_crit4,cis_crit4,make_plot_df)

##function for plots

plot_slopes_TP <- function(df){
  df |>
    ggplot(aes(x = fct_reorder(baseline, estimate, median), y = estimate, fill = baseline)) +
    geom_linerange(aes(ymin = X1, ymax = X2)) +
    geom_point(size = 2, color = "black", shape = 21, show.legend = F) +
    geom_hline(yintercept = 0, color = "black", linewidth = .75, linetype = "dashed") +
    scale_fill_viridis_d()+
    coord_flip()+
    scale_y_continuous(limits = c(-.7,.7), breaks = seq(-.7,.7,.35))+
    labs(y = NULL) + 
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

###ffg plots
p13 <- plot_data_crit4[[1]]|>filter(taxon_code == "BNT")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Feeding Group")
p14 <- plot_data_crit4[[1]]|>filter(taxon_code == "CKC")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Feeding Group")
p15 <- plot_data_crit4[[1]]|>filter(taxon_code == "LND")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Feeding Group")
p16 <- plot_data_crit4[[1]]|>filter(taxon_code == "LNS")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Feeding Group")
p17 <- plot_data_crit4[[1]]|>filter(taxon_code == "WHS")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Feeding Group",
       y = "TP")

#taxa plots
p18 <- plot_data_crit4[[2]]|>filter(taxon_code == "BNT")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Taxonomic Group")
p19 <- plot_data_crit4[[2]]|>filter(taxon_code == "CKC")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Taxonomic Group")
p20 <- plot_data_crit4[[2]]|>filter(taxon_code == "LND")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Taxonomic Group")
p21 <- plot_data_crit4[[2]]|>filter(taxon_code == "LNS")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Taxonomic Group")
p22 <- plot_data_crit4[[2]]|>filter(taxon_code == "WHS")|>
  plot_slopes_TP()+
  facet_wrap(~taxon_code)+
  labs(x = "Taxonomic Group",
       y = "TP")

######scatter plots to show data


TP_ffg_plot <- TP_ffg %>% 
  bind_rows(TP_uncorrected)

TP_taxa_plot <- TP_taxa %>% 
  bind_rows(TP_uncorrected)

TP_ffg_plot |> 
  ggplot(aes(x = PC1, y = TP, color = baseline))+
  facet_wrap(~taxon_code)+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


TP_taxa_plot |> 
  ggplot(aes(x = PC1, y = TP, color = baseline))+
  facet_wrap(~taxon_code)+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

plot_TP_scatter <- function(df){
df |> 
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
  facet_wrap(~taxon_code)+
  scale_color_viridis_d()+ 
  scale_fill_viridis_d()+
  labs(fill = "Taxonomic Group", y = "TP",
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  )
}

#ffg scatters
p23 <- TP_ffg_plot |> filter(taxon_code == "BNT") |>plot_TP_scatter()
p24 <- TP_ffg_plot |> filter(taxon_code == "CKC") |>plot_TP_scatter()
p25 <- TP_ffg_plot |> filter(taxon_code == "LND") |>plot_TP_scatter()
p26 <- TP_ffg_plot |> filter(taxon_code == "LNS") |>plot_TP_scatter()
p27 <- TP_ffg_plot |> filter(taxon_code == "WHS") |>plot_TP_scatter()+
  labs(x = "PC1")

###taxa scatters

p28 <- TP_taxa_plot |> filter(taxon_code == "BNT") |>plot_TP_scatter()
p29 <- TP_taxa_plot |> filter(taxon_code == "CKC") |>plot_TP_scatter()
p30 <- TP_taxa_plot |> filter(taxon_code == "LND") |>plot_TP_scatter()
p31 <- TP_taxa_plot |> filter(taxon_code == "LNS") |>plot_TP_scatter()
p32 <- TP_taxa_plot |> filter(taxon_code == "WHS") |>plot_TP_scatter()+
  labs(x = "PC1")

patch.crit4 <- p18+p28+p13+p23+p19+p29+p14+p24+p20+p30+p15+p25+p21+p31+p16+p26+p22+p32+p17+p27

patch.crit4.annote <- patch.crit4+
  plot_layout(ncol = 4,widths = c(1,2,1,2))& 
  plot_annotation(tag_levels = "A")
patch.crit4.annote

ggsave(here("out", "fig4_crit4.png"), 
       patch.crit4.annote, device = ragg::agg_png,
       units = "in", width = 20, height = 20)
###ancovas###--------------------

#ancovas for Tp corrected by taxa

ancovataxa <- TP_taxa |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(TP_corrected ~ PC1*baseline, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

ancova_taxa |> unnest(glance)

#ancovas for Tp corrected by ffg

ancova_ffg <- TP_ffg |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(TP_corrected ~ PC1*baseline, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

temp <- ancova_ffg |> unnest(tidy)
