
# Criteria 4: Trophic Estimates independent of Land Use

## Prep

## Uncorrected TP ------------------

# calculate uncorrected TP
TP_uncorrected <- isodat_fish |> 
  mutate(
    discFactor = ((-0.281*d15N)+5.879),
    TP = d15N/discFactor,
    baseline = "uncorrected"
  )

# plot uncorrected TPs
TP_uncorrected |> 
  ggplot(aes(x = PC1, y = TP, color = taxon_code))+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


## Corrected TP ------------------

### Taxa --------------------

# calculate mean d15N (i.e., 'correction') for each baseline at each site

taxa_correction <- isodat_bug_common |>
  group_by(site_id, taxon_code) |>
  summarise(
    correction = mean(d15N, na.rm = TRUE),
    correction_type = "taxa", 
    .groups = "drop"
  ) |> 
  rename(baseline = taxon_code)

basal_taxa_correction <- isodat_baseline |>
  filter(taxon_code != "macrophtye") |>
  group_by(site_id, taxon_code) |>
  summarise(
    correction = mean(d15N,na.rm = TRUE),
    correction_type = "taxa", 
    .groups = "drop"
  ) |>  
  rename(baseline = taxon_code)

# use corrections (i.e., mean d15N) at each site to calculate corrected 
# trophic positions for each fish

TP_taxa <- isodat_fish |>
  left_join(taxa_correction, by = "site_id") |>
  mutate(
    discFactor = ((-0.281*d15N)+5.879),
    TP = (((d15N-correction)/discFactor)+2)
  ) 

TP_taxa_basal <- isodat_fish |>
  left_join(basal_taxa_correction, by = "site_id")|>
  mutate(
    discFactor = ((-0.281*d15N)+5.879),
    TP = (((d15N-correction)/discFactor)+1)
  ) 

# combine the two dfs
TP_taxa <- TP_taxa |> bind_rows(TP_taxa_basal)

# test plot
TP_taxa |> 
  ggplot(aes(x = PC1, y = TP, color = baseline))+
  facet_wrap(~taxon_code)+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 



### FFGs -----------------

# calculate mean d15N (i.e., 'correction') for each ffg baseline at each site

ffg_correction <- isodat_bug |>
  filter(ffg != "Shredder")|>
  group_by(site_id, ffg)|>
  summarise(correction = mean(d15N, na.rm = T),
            correction_type = "ffg") %>% 
  rename(baseline = ffg)


bulk_invert_correction <- isodat_bug |> 
  group_by(site_id)|>
  summarise(correction = mean(d15N,na.rm = T),
            correction_type = "ffg",
            baseline = "Average Invert")

ffg_correction <- ffg_correction |> bind_rows(bulk_invert_correction)

bulk_basal_correction <- isodat_baseline |>
  filter(compartment == "baseline")|>
  filter(taxon_code != "macrophtye")|>
  group_by(site_id)|>
  summarise(correction = mean(d15N,na.rm = T),
            correction_type = "ffg",
            baseline = "Average Basal")

#use mean d15N of the baseline (taxa) at each site to calculate corrected 
#Trophic positions

TP_ffg <- isodat_fish |> 
  left_join(ffg_correction, by = "site_id", relationship = "many-to-many")|>
  mutate(discFactor = ((-0.281*d15N)+5.879),
         TP = (((d15N-correction)/discFactor)+2))

TP_bulk_basal <- isodat_fish |>
  left_join(bulk_basal_correction, by = "site_id", relationship = "many-to-many")|> 
  mutate(discFactor = ((-0.281*d15N)+5.879),
         TP = (((d15N-correction)/discFactor)+1))

TP_ffg <- TP_ffg |> bind_rows(TP_bulk_basal)

# test plot
TP_ffg |> 
  ggplot(aes(x = PC1, y = TP, color = baseline))+
  facet_wrap(~taxon_code)+
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


## Fit models -----

### uncorrected ----
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

### TPs corrected by ffgs -------

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


### TPs corrected by taxa -------

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


## Slope plots -----

### Prep dataframes ------
get_cis <- function(data) {
  temp <- data |> 
    unnest(tidy)|>
    unnest(conf_int) |> 
    filter(term == "PC1")
  out <- data.frame((unlist(temp$conf_int)))
  print(out)
}

add_sig_col <- function(df_cis) {
  df_cis |> 
    mutate(sig = case_when(
      X1 < 0 & X2 > 0 ~ "no", 
      TRUE ~ "yes"
    ))
}

make_plot_df <- function(model, cis) {
  out <- model |> 
    bind_cols(cis) |> 
    unnest(tidy)|> 
    filter(term == "PC1") 
  out
}

# make list of mods
dfs_crit4 <- list(mods_corrected_ffg, mods_corrected_taxa)

# extract CIs
cis_crit4 <- map(dfs_crit4, get_cis)

# add significance if CIs don't overlap zero
# cis_crit4[[1]]
cis_crit4 <- map(cis_crit4, add_sig_col)

plot_data_crit4 <- map2(dfs_crit4, cis_crit4, make_plot_df)

##function for plots

plot_slopes_TP <- function(df){
  df |>
    # ggplot(aes(x = fct_reorder(baseline, estimate, median), y = estimate, fill = baseline)) +
    ggplot(aes(x = baseline, y = estimate, fill = baseline)) +
    geom_linerange(aes(ymin = X1, ymax = X2)) +
    geom_point(size = 2, color = "black", shape = 21, show.legend = F) +
    geom_hline(yintercept = 0, color = "black", linewidth = .75, linetype = "dashed") +
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
      plot.margin = unit(c(0,.1,0,0), "cm")
    )
}

### ffg plots ------

p13 <- plot_data_crit4[[1]] |> filter(taxon_code == "BNT") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "Average Basal",
        "Average Invert", 
        "Omnivore",
        "Collector", 
        "Filterer", 
        "Predator", 
        "Grazer"
        )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Feeding Group") +
  scale_fill_manual(values = palette_ffgs)

p14 <- plot_data_crit4[[1]] |> filter(taxon_code == "CKC") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected",
        "Average Basal", 
        "Average Invert", 
        "Grazer",
        "Omnivore",
        "Collector", 
        "Filterer", 
        "Predator"
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Feeding Group") +
  scale_fill_manual(values = palette_ffgs)

p15 <- plot_data_crit4[[1]] |> filter(taxon_code == "LND")|>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |>
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "Average Basal",
        "Average Invert", 
        "Filterer", 
        "Predator",
        "Omnivore", 
        "Collector", 
        "Grazer"
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Feeding Group") +
  scale_fill_manual(values = palette_ffgs)

p16 <- plot_data_crit4[[1]] |> filter(taxon_code == "LNS")|>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "Average Basal",
        "Average Invert", 
        "Predator", 
        "Omnivore", 
        "Collector", 
        "Filterer", 
        "Grazer"
      )))) |> 
  plot_slopes_TP() + 
  facet_wrap(~common_name) +
  labs(x = "Feeding Group") +
  scale_fill_manual(values = palette_ffgs)

p17 <- plot_data_crit4[[1]] |> filter(taxon_code == "WHS")|>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "Average Basal",
        "Average Invert", 
        "Grazer",
        "Omnivore", 
        "Filterer", 
        "Collector", 
        "Predator" 
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Feeding Group",
       y = expression('slope ('~{beta}[1]~')'))+
  scale_fill_manual(values = palette_ffgs)

### taxa plots -------

p18 <- plot_data_crit4[[2]] |> filter(taxon_code == "BNT") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "FBOM",
        "Seston", 
        "filamentous",
        "Biofilm", 
        "Ephemeridae", 
        "Elmidae-larvae", 
        "Dytiscidae", 
        "Hydropyschidae", 
        "Leptohyphidae",
        "Elmidae-adult", 
        "Simuliidae", 
        "Gomphidae", 
        "Baetidae", 
        "Perlidae", 
        "Heptaganeidae", 
        "Chironomidae"
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Taxonomic Group") +
  scale_fill_manual(values = palette_taxon)

p19 <- plot_data_crit4[[2]] |> filter(taxon_code == "CKC") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "FBOM",
        "filamentous",
        "Biofilm", 
        "Seston", 
        "Baetidae", 
        "Hydropyschidae", 
        "Chironomidae",
        "Heptaganeidae", 
        "Elmidae-adult", 
        "Simuliidae", 
        "Elmidae-larvae", 
        "Perlidae", 
        "Leptohyphidae",
        "Ephemeridae", 
        "Gomphidae", 
        "Dytiscidae"
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Taxonomic Group") +
  scale_fill_manual(values = palette_taxon)

p20 <- plot_data_crit4[[2]] |> filter(taxon_code == "LND") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "FBOM",
        "filamentous",
        "Biofilm", 
        "Seston", 
        "Dytiscidae",
        "Elmidae-larvae", 
        "Simuliidae", 
        "Hydropyschidae", 
        "Elmidae-adult", 
        "Perlidae", 
        "Gomphidae", 
        "Heptaganeidae", 
        "Baetidae", 
        "Ephemeridae", 
        "Chironomidae",
        "Leptohyphidae"
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Taxonomic Group") +
  scale_fill_manual(values = palette_taxon)

p21 <- plot_data_crit4[[2]] |> filter(taxon_code == "LNS") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "FBOM",
        "filamentous",
        "Seston", 
        "Biofilm", 
        "Gomphidae", 
        "Dytiscidae",
        "Elmidae-adult", 
        "Leptohyphidae",
        "Elmidae-larvae", 
        "Perlidae", 
        "Simuliidae", 
        "Hydropyschidae", 
        "Ephemeridae", 
        "Chironomidae",
        "Heptaganeidae", 
        "Baetidae" 
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(x = "Taxonomic Group") +
  scale_fill_manual(values = palette_taxon)

p22 <- plot_data_crit4[[2]] |> filter(taxon_code == "WHS") |>
  left_join(isodat_fish |> select(taxon_code, common_name), by = "taxon_code") |> 
  mutate(
    baseline = fct_rev(factor(
      baseline, 
      levels = c(
        "uncorrected", 
        "FBOM",
        "filamentous",
        "Seston", 
        "Biofilm", 
        "Chironomidae",
        "Baetidae",
        "Hydropyschidae", 
        "Elmidae-adult", 
        "Simuliidae", 
        "Heptaganeidae", 
        "Gomphidae", 
        "Perlidae", 
        "Elmidae-larvae", 
        "Leptohyphidae",
        "Dytiscidae",
        "Ephemeridae"
      )))) |> 
  plot_slopes_TP() +
  facet_wrap(~common_name) +
  labs(
    x = "Taxonomic Group",
    y = expression('slope ('~{beta}[1]~')')
    )+
  scale_fill_manual(values = palette_taxon)



## Scatter plots to show data -------


TP_ffg_plot <- TP_ffg |> bind_rows(TP_uncorrected)

TP_taxa_plot <- TP_taxa |> bind_rows(TP_uncorrected)

# TP_ffg_plot |> 
#   ggplot(aes(x = PC1, y = TP, color = baseline))+
#   facet_wrap(~taxon_code)+
#   geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) 
# 
# TP_taxa_plot |> 
#   ggplot(aes(x = PC1, y = TP, color = baseline))+
#   facet_wrap(~taxon_code)+
#   geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) 

# plot_TP_scatter <- function(df) {
#   df |> 
#     ggplot(aes(PC1, TP, fill = baseline)) + 
#     geom_smooth(aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
#     geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
#     scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
#     #scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
#     facet_wrap(~taxon_code)+  
#     labs(fill = "Taxonomic Group", y = "TP",
#          x = NULL)+
#     theme_bw(base_size = 12) + 
#     theme(
#       axis.line = element_line(size = .5),
#       axis.ticks.length = unit(.25, "cm"), 
#       panel.grid.minor.x = element_blank(), 
#       plot.margin = unit(c(0,0,0,0), "cm")
#     )
# }

### ffg scatters ------

plot_data_crit4[[1]] |> filter(taxon_code == "BNT")
p23 <- TP_ffg_plot |> filter(taxon_code == "BNT") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  #scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP",
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  scale_color_manual(values=palette_ffgs) + 
  scale_fill_manual(values=palette_ffgs) 

plot_data_crit4[[1]] |> filter(taxon_code == "CKC") |> select(1,2,15)
p24 <- TP_ffg_plot |> filter(taxon_code == "CKC") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_ffg_plot |> 
      filter(taxon_code == "CKC") |> 
      filter(! baseline %in% c("Collector","Filterer","Average Invert")),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  #scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP",
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  scale_color_manual(values=palette_ffgs) + 
  scale_fill_manual(values=palette_ffgs)  

plot_data_crit4[[1]] |> filter(taxon_code == "LND") |> select(1,2,15)
p25 <- TP_ffg_plot |> filter(taxon_code == "LND") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_ffg_plot |> 
      filter(taxon_code == "LND") |> 
      filter(! baseline %in% c(
        "Collector","Filterer","Average Invert","Grazer","Omnivore","Predator")),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  #scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP",
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  scale_color_manual(values=palette_ffgs) + 
  scale_fill_manual(values=palette_ffgs)  

plot_data_crit4[[1]] |> filter(taxon_code == "LNS") |> select(1,2,15)
p26 <- TP_ffg_plot |> filter(taxon_code == "LNS") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_ffg_plot |> 
      filter(taxon_code == "LNS") |> 
      filter(! baseline %in% c("Collector","Filterer","Average Invert","Grazer","Omnivore","Predator")),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  #scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP",
       x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  scale_color_manual(values=palette_ffgs) + 
  scale_fill_manual(values=palette_ffgs) 

plot_data_crit4[[1]] |> filter(taxon_code == "WHS") |> select(1,2,15)
p27 <- TP_ffg_plot |> filter(taxon_code == "WHS") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_ffg_plot |> 
      filter(taxon_code == "WHS") |> 
      filter(! baseline %in% c("Filterer","Grazer","Omnivore")),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F)+
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  #scale_y_continuous(limits = c(1, 8.3), breaks = seq(0,8,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP",
       x = "PC1")+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")
  ) + 
  scale_color_manual(values=palette_ffgs) + 
  scale_fill_manual(values=palette_ffgs) 

### taxa scatters ----

plot_data_crit4[[2]] |> filter(taxon_code == "BNT") |> select(1,2,15)
p28 <- TP_taxa_plot |> filter(taxon_code == "BNT") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_taxa_plot |> filter(taxon_code == "BNT") |>
      filter(! baseline %in% c(
        "Chironomidae",
        "Heptaganeidae"
        )),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F) +
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP", x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_color_manual(values=palette_taxon) + 
  scale_fill_manual(values=palette_taxon) 

plot_data_crit4[[2]] |> filter(taxon_code == "CKC") |> select(1,2,15)
p29 <- TP_taxa_plot |> filter(taxon_code == "CKC") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_taxa_plot |> filter(taxon_code == "CKC") |>
      filter(! baseline %in% c(
        "Elmidae-adult",
        "Elmidae-larvae",
        "Ephemeridae",
        "Heptaganeidae",
        "Leptohyphidae",
        "Perlidae",
        "Simuliidae",
        "Seston"
      )),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F) +
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP", x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_color_manual(values=palette_taxon) + 
  scale_fill_manual(values=palette_taxon)  

plot_data_crit4[[2]] |> filter(taxon_code == "LND") |> select(1,2,15)
p30 <- TP_taxa_plot |> filter(taxon_code == "LND") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_taxa_plot |> filter(taxon_code == "LND") |>
      filter(! baseline %in% c(
        "Baetidae",
        "Dytiscidae",
        "Elmidae-adult",
        "Elmidae-larvae",
        "Gomphidae",
        "Heptaganeidae",
        "Hydropyschidae",
        "Perlidae",
        "Simuliidae",
        "Biofilm",
        "Seston"
      )),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F) +
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP", x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_color_manual(values=palette_taxon) + 
  scale_fill_manual(values=palette_taxon)  

plot_data_crit4[[2]] |> filter(taxon_code == "LNS") |> select(1,2,15)
p31 <- TP_taxa_plot |> filter(taxon_code == "LNS") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_taxa_plot |> filter(taxon_code == "LNS") |>
      filter(! baseline %in% c(
        "Baetidae",
        "Chironomidae",
        "Elmidae-adult",
        "Elmidae-larvae",
        "Ephemeridae",
        "Heptaganeidae",
        "Hydropyschidae",
        "Perlidae",
        "Simuliidae",
        "Leptohyphidae"
      )),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F) +
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP", x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_color_manual(values=palette_taxon) + 
  scale_fill_manual(values=palette_taxon)  


plot_data_crit4[[2]] |> filter(taxon_code == "WHS") |> select(1,2,15)
p32 <- TP_taxa_plot |> filter(taxon_code == "WHS") |>
  ggplot(aes(PC1, TP, fill = baseline)) + 
  geom_smooth(
    data = TP_taxa_plot |> filter(taxon_code == "WHS") |>
      filter(! baseline %in% c(
        "Baetidae",
        "Chironomidae",
        "Elmidae-adult",
        "Hydropyschidae",
        "Simuliidae",
        "Biofilm",
        "Seston"
      )),
    aes(color = baseline), method = "lm", se = FALSE, show.legend = F) +
  geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
  scale_x_continuous(limits = c(1.4, 10), breaks = seq(0,10,2))+
  facet_wrap(~common_name)+  
  labs(fill = "Taxonomic Group", y = "TP", x = NULL)+
  theme_bw(base_size = 12) + 
  theme(
    axis.line = element_line(size = .5),
    axis.ticks.length = unit(.25, "cm"), 
    panel.grid.minor.x = element_blank(), 
    plot.margin = unit(c(0,0,0,0), "cm")) + 
  scale_color_manual(values=palette_taxon) + 
  scale_fill_manual(values=palette_taxon)  


## Combine plots --------
patch.crit4.taxa <- p18+p28+p19+p29+p20+p30+p21+p31+p22+p32
patch.crit4.ffg <- p13+p23+p14+p24+p15+p25+p16+p26+p17+p27

patch.crit4.taxa.annote <- patch.crit4.taxa+
  plot_layout(ncol = 2,widths = c(2,3))& 
  plot_annotation(tag_levels = "a")
patch.crit4.taxa.annote

patch.crit4.ffg.annote <- patch.crit4.ffg+
  plot_layout(ncol = 2,widths = c(2,3))& 
  plot_annotation(tag_levels = "a")
patch.crit4.ffg.annote


path <- here::here("out", "criteria-4-taxonomic_groups")
ggsave(glue::glue("{path}.pdf"), plot = patch.crit4.taxa.annote, 
       width = 4, height = 6, scale = 2.5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)

path <- here::here("out", "criteria-4-feeding_groups")
ggsave(glue::glue("{path}.pdf"), plot = patch.crit4.ffg.annote, 
       width = 4, height = 6, scale = 2.5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)
