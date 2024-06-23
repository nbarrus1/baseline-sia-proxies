
# Criteria #2: site-level d15N CV for most common taxa and ffgs


# Pull widely distributed  the data so only 75% distribution or higher is kept
focustaxon <- Distr_dat |> filter(distr >= 75) |> pull(taxon_code)
foucusffg <- Distr_dat_FFG |> filter(distr >= 75) |> pull(ffg)

# Only keep data from widely distributed taxa
isodat_bug_common <- isodat_bug |> filter(taxon_code %in% focustaxon)

# invert_group <- invert_group |> 
#   mutate(taxon_code = if_else(taxon_code == "Simulidae", true = "Simuliidae",
#                               false = taxon_code))

## Calculate d15N CV ----------

cv_taxa <- isodat_bug_common %>%
  group_by(site_id, taxon_code) %>%
  summarise(
    n = n(),
    mean_d15n = mean(d15N, na.rm = TRUE),
    sd_d15n = sd(d15N, na.rm = TRUE),
    CV_d15n = sd_d15n/mean_d15n,
    .groups = 'drop'
  ) |> 
  left_join(isodat_bug |> select(taxon_code, ffg), by = "taxon_code") |> 
  arrange(ffg, taxon_code) |> 
  drop_na(CV_d15n)|>
  mutate(log_CV = log(CV_d15n))


cv_ffg <- isodat_bug_common %>%
  group_by(site_id, ffg) %>%
  summarise(
    n = n(),
    mean_d15n = mean(d15N, na.rm = TRUE),
    sd_d15n = sd(d15N, na.rm = TRUE),
    CV_d15n = sd_d15n/mean_d15n,
    .groups = 'drop'
  ) |> 
  drop_na(CV_d15n) |> 
  mutate(log_CV = log(CV_d15n))


# plot data 

p3 <- cv_taxa |> 
    ggplot(aes(x = fct_reorder(taxon_code, CV_d15n, median), y = CV_d15n, fill = ffg)) +
    geom_boxplot() + 
    scale_fill_manual(values=palette_ffgs) + 
    labs(y=expression('CV (' ~{delta}^15*N~')'), x="", fill="Feeding group") + 
    theme_bw(base_size = 10) + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5))

p4 <- cv_ffg |>
    ggplot(aes(x = fct_reorder(ffg, CV_d15n, median), y = CV_d15n, fill = ffg)) +
    geom_boxplot() + 
  scale_fill_manual(values=palette_ffgs) + 
  labs(y=expression('CV (' ~{delta}^15*N~')'), x="", fill="Feeding group") + 
    theme_bw(base_size = 10)+ 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5))

patch.cv <- p3 | p4
patch.cv.annote <- patch.cv + 
  plot_layout(widths = c(1.75, 1))+ 
  plot_layout(guides = 'collect') & 
  plot_annotation(tag_levels = "a")
patch.cv.annote

# Save plot
path <- here::here("out", "criteria-2-raw")
ggsave(glue::glue("{path}.pdf"), plot = patch.cv.annote, 
       width = 3, height = 1.75, scale = 3, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)

## Compare CVs  ------------------------------

### Taxa -----------------------------------------

# normality
cv_taxa |> ggplot(aes(CV_d15n)) + geom_density()
cv_taxa |> ggplot(aes(log_CV)) + geom_density()
shapiro.test(cv_taxa$CV_d15n)  # nope
shapiro.test(cv_taxa$log_CV)  # yes
# logged data is normal

# test for heteroscedasticity
car::leveneTest(cv_taxa$log_CV, group = factor(cv_taxa$taxon_code)) 
# we have heteroscedasticity, but how bad?

# check variance of groups (max no more than 4x min var)
cv_taxa_var <- cv_taxa |> 
  group_by(taxon_code) |> 
  summarise(variance = var(log_CV)) |> 
  arrange(variance)

# ratios
max(cv_taxa_var$variance) / min(cv_taxa_var$variance)  
# ~15x so a lot
# kruskal best way to prevent type 1 error, compare with anova

# Fit models
kruskal.test(log_CV ~ taxon_code, data = cv_taxa)
compare_means(log_CV ~ taxon_code, data = cv_taxa, method = "kruskal")

compare_means(log_CV ~ taxon_code, data = cv_taxa, method = "anova")
fit_taxa <- lm(log(CV_d15n) ~ taxon_code, data = cv_taxa)
anova(fit_taxa)
agricolae::HSD.test(fit_taxa,'taxon_code',group = FALSE,console = TRUE, unbalanced = TRUE)


### FFGs -------------------------------------------

#  check normality
cv_ffg |> ggplot(aes(CV_d15n)) + geom_density()
cv_ffg |> ggplot(aes(log_CV)) + geom_density()
shapiro.test(cv_ffg$CV_d15n)
shapiro.test(cv_ffg$log_CV)  # good 

# check variances
car::leveneTest(cv_ffg$log_CV, group = factor(cv_ffg$ffg)) 
# good, but check variances anyway

# check variance of groups (max no more than 4x min var)
cv_ffg_var <- cv_ffg |> 
  group_by(ffg) |> 
  summarise(variance = var(log_CV)) |> 
  arrange(variance)
# ratio
max(cv_ffg_var$variance) / min(cv_ffg_var$variance)  
# ~6x so not bad
# ANOVA is fine but comapre to kruskal anyway

# Fit models
kruskal.test(log_CV ~ ffg, data = cv_ffg)
compare_means(log_CV ~ ffg, data = cv_ffg, method = "kruskal")

compare_means(log_CV ~ ffg, data = cv_ffg, method = "anova")
fit_ffg <- lm(log(CV_d15n) ~ ffg, data = cv_ffg)
anova(fit_ffg)
agricolae::HSD.test(fit_ffg,'ffg',group = FALSE,console = TRUE, unbalanced = TRUE)


### Visualize model results------------

# taxa
taxa_comparisons = list(
  c("Simuliidae", "Ephemeridae"),
  c("Simuliidae", "Perlidae"),
  c("Simuliidae", "Hydropyschidae"),
  c("Simuliidae", "Gomphidae"),
  c("Simuliidae", "Leptohyphidae"),
  c("Simuliidae", "Dytiscidae"),
  c("Simuliidae", "Chironomidae"),
  c("Simuliidae", "Elmidae-larvae"),
  c("Simuliidae", "Elmidae-adult"),
  c("Gomphidae", "Heptaganeidae"),
  c("Chironomidae", "Heptaganeidae"),
  c("Elmidae-larvae", "Heptaganeidae"),
  c("Elmidae-adult", "Heptaganeidae"),
  c("Heptaganeidae", "Leptohyphidae")
)

p3_stats <- cv_taxa |>
  mutate(ffg = factor(ffg, levels = c("Filterer","Grazer","Omnivore","Predator","Collector"))) |> 
  mutate(taxon_code = fct_reorder(taxon_code, log_CV, median)) |>
  ggboxplot(x = "taxon_code", y = "log_CV", fill = "ffg") +
  stat_compare_means(comparisons = taxa_comparisons, label.y = seq(0,10.5,.75),label = "p.signif")+
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 11)  +
  labs(y=expression('log [CV (' ~{delta}^15*N~')]'), x="", fill = "Feeding Group") +
  scale_fill_manual(values=palette_ffgs) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5))

# ffgs

ffg_comparisons = list(
  c("Filterer", "Omnivore"),
  c("Filterer", "Grazer"),
  c("Filterer", "Predator"),
  c("Filterer", "Collector"),
  c("Omnivore", "Collector"),
  c("Grazer", "Collector")
)

p4_stats <- cv_ffg |>
  mutate(ffg = fct_reorder(ffg, log_CV, median)) |>
  ggboxplot(x = "ffg", y = "log_CV", fill = "ffg") +
  stat_compare_means(comparisons = ffg_comparisons, label.y = seq(0,4.5,.75),label = "p.signif")+
  stat_compare_means(method = "anova", label.x = 1.5, label.y = 4.5)  +
  labs(
  y=expression('log [CV (' ~{delta}^15*N~')]'), x="", fill = "Feeding Group") +
  scale_fill_manual(values=palette_ffgs) + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5))

patch.cv.stats <- p3_stats | p4_stats
patch.cv.stats.annote <- patch.cv.stats + 
  plot_layout(widths = c(2, 1))+ 
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom') & 
  plot_annotation(tag_levels = "a")
patch.cv.stats.annote

# Save plot
path <- here::here("out", "criteria-2-stats")
ggsave(glue::glue("{path}.pdf"), plot = patch.cv.stats.annote, 
       width = 4, height = 3, scale = 2.5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)



## Old plts -----------
# prep data for plotting (summarize)
# 
# cv_taxa_means <- cv_taxa %>% 
#   # drop_na(CV_d15n) %>% # 112 to 90
#   group_by(taxon_code) %>% 
#   summarise(
#     n = n(),
#     mean_CV = mean(CV_d15n, na.rm = TRUE),
#     sd_CV = sd(CV_d15n, na.rm = TRUE)
#   ) %>% 
#   mutate(
#     # se = sd_CV/sqrt(n),
#     # low = lower_ci(mean_CV,se,n),
#     # high = upper_ci(mean_CV,se,n)
#     sd_low = mean_CV - sd_CV,
#     sd_high = mean_CV + sd_CV
#   ) %>%
#   left_join(invert_group, by = "taxon_code") |> 
#   rename(group = ffg)
# 
# cv_ffg_means <- cv_ffg %>%
#   drop_na(CV_d15n) %>%
#   group_by(ffg) %>% 
#   summarise(
#     n = n(),
#     mean_CV = mean(CV_d15n, na.rm = TRUE),
#     sd_CV = sd(CV_d15n, na.rm = TRUE)
#     ) %>% 
#   mutate(
#     # se = sd_CV/sqrt(n),
#     # low = lower_ci(mean_CV,se,n),
#     # high = upper_ci(mean_CV,se,n)
#     sd_low = mean_CV - sd_CV,
#     sd_high = mean_CV + sd_CV
#   ) |> 
#   mutate(group = ffg)
# 
# ## Plot
# 
# # plot function
# plot_cvs <- function(df, x_cat){
#   df |>
#     ggplot(aes(x = fct_reorder({{ x_cat }}, mean_CV, median), y = mean_CV, fill = group)) +
#     geom_linerange(ymin = df$sd_low, ymax = df$sd_high) +
#     geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
#     coord_flip()+
#     scale_y_continuous(limits = c(0, .63), breaks = seq(0,.6,.2))  +
#     scale_fill_manual(values = c("#440154FF","#404788FF","#287D8EFF",
#                                  "#29AF7FFF","#95D840FF"))+
#     labs(y = "", fill = "Feeding Group") + 
#     theme_bw(base_size = 12) + 
#     theme(
#       axis.line = element_line(size = .5),
#       axis.ticks.length = unit(.25, "cm"), 
#       panel.grid.minor.x = element_blank(), 
#       plot.margin = unit(c(0,0,0,0), "cm")
#     )
# }
# 
# p3 <- cv_taxa_means |> plot_cvs(taxon_code) + labs(x = NULL, y = expression('CV (' ~{delta}^15*N~')'))
# 
# p4 <- cv_ffg_means |> plot_cvs(ffg) + labs(x = NULL, y = expression('CV (' ~{delta}^15*N~')'))
# patch.cv <- p3 / p4
# 
# patch.cv.annote <- patch.cv + 
#   plot_layout(heights = c(3, 1))
# 
# patch.cv.annote
# 
# patch.fig2 <- patch.dist.annote | patch.cv.annote
# patch.fig2.annotate<- patch.fig2 + 
#   plot_layout(guides = 'collect') & 
#   plot_annotation(tag_levels = "A")
# 
# # save plot
# # ggsave(here("out", "fig2_crit1&crit2.png"), 
# #        patch.fig2.annotate, device = ragg::agg_png,
# #        units = "in", width = 10, height = 6)
# 
# 
# 
