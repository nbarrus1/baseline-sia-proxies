
# Criteria #2: site-level CV for most common taxa

## Data -----------

# filter the data so only 75% distribution or higher is kept
focustaxon <- Distr_dat %>%
  filter(distr >= 75)

foucusffg <- Distr_dat_FFG %>% 
  filter(distr >= 75)

invertdata_sub <- invertdata %>%
  filter(d15N >= 0) |>   # must remove negative data for CV calculations n=3
  filter(taxon_code %in% focustaxon$taxon_code)

## CV by site and group ----------

cv_taxa <- invertdata_sub %>%
  group_by(site_id, taxon_code) %>%
  summarise(
    n = n(),
    mean_d15n = mean(d15N, na.rm = TRUE),
    sd_d15n = sd(d15N, na.rm = TRUE),
    CV_d15n = sd_d15n/mean_d15n, 
    .groups = 'drop'
  ) |> 
  left_join(invert_group, by = "taxon_code") |> 
  arrange(ffg, taxon_code)

cv_ffg <- invertdata %>%
  filter(ffg != "Shredder") %>% 
  group_by(site_id, ffg) %>%
  summarise(
    n = n(),
    mean_d15n = mean(d15N, na.rm = TRUE),
    sd_d15n = sd(d15N, na.rm = TRUE),
    CV_d15n = sd_d15n/mean_d15n, 
    .groups = 'drop'
  ) 


### Summary tables for CV data ------------

cv_taxa |>
  group_by(taxon_code) %>% 
  summarise(CV = mean(CV_d15n),
            SD = sd(CV_d15n))

cv_ffg |>
  group_by(ffg) %>% 
  summarise(CV = mean(CV_d15n),
            SD = sd(CV_d15n))


## Summarize mean CV  ---------

cv_taxa_means <- cv_taxa %>% 
  drop_na(CV_d15n) %>% # 112 to 90
  # filter(CV_d15n > 0 & CV_d15n < 1) %>%  # 90 to 87 WHAT DOS THIS REMOVE????
  group_by(taxon_code) %>% 
  summarise(
    n = n(),
    mean_CV = mean(CV_d15n, na.rm = TRUE),
    sd_CV = sd(CV_d15n)
    ) %>% 
  mutate(
    se = sd_CV/sqrt(n),
    # low = lower_ci(mean_CV,se,n),
    # high = upper_ci(mean_CV,se,n)
    low = mean_CV - sd_CV,
    high = mean_CV + sd_CV
    ) %>%
  left_join(invert_group, by = "taxon_code") |> 
  rename(group = ffg)


cv_ffg_means <- cv_ffg %>%
  drop_na(CV_d15n) %>%
  group_by(ffg) %>% 
  summarise(
    n = n(),
    mean_CV = mean(CV_d15n, na.rm = TRUE),
    sd_CV = sd(CV_d15n)
    ) %>% 
  mutate(
    se = sd_CV/sqrt(n),
    # low = lower_ci(mean_CV,se,n),
    # high = upper_ci(mean_CV,se,n)
    low = mean_CV - sd_CV,
    high = mean_CV + sd_CV
    ) |> 
  mutate(group = ffg)

## Plot -----------

# plot function
plot_cvs <- function(df, x_cat){
  df |>
    ggplot(aes(x = fct_reorder({{ x_cat }}, mean_CV, median), y = mean_CV, fill = group)) +
    geom_linerange(ymin = df$low, ymax = df$high) +
    geom_point(size = 2, color = "black", shape = 21, show.legend = F) + 
    coord_flip()+
    scale_y_continuous(limits = c(0, .63), breaks = seq(0,.6,.2))  +
    scale_fill_manual(values = c("#440154FF","#404788FF","#287D8EFF",
                                            "#29AF7FFF","#95D840FF"))+
    labs(y = "", fill = "Feeding Group") + 
    theme_bw(base_size = 12) + 
    theme(
      axis.line = element_line(size = .5),
      axis.ticks.length = unit(.25, "cm"), 
      panel.grid.minor.x = element_blank(), 
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

p3 <- cv_taxa_means |> plot_cvs(taxon_code) + labs(x = NULL, y = expression('CV (' ~{delta}^15*N~')'))

p4 <- cv_ffg_means |> plot_cvs(ffg) + labs(x = NULL, y = expression('CV (' ~{delta}^15*N~')'))
patch.cv <- p3 / p4

patch.cv.annote <- patch.cv + 
  plot_layout(heights = c(3, 1))

patch.cv.annote

patch.fig2 <- patch.dist.annote | patch.cv.annote
patch.fig2.annotate<- patch.fig2 + 
  plot_layout(guides = 'collect') & 
  plot_annotation(tag_levels = "A")

# save plot
ggsave(here("out", "fig2_crit1&crit2.png"), 
       patch.fig2.annotate, device = ragg::agg_png,
       units = "in", width = 10, height = 6)




## ANOVA of CVs ------------------------------

# Taxa 


cv_taxa |> ggplot(aes(CV_d15n)) + geom_density()
cv_taxa |> ggplot(aes(log(CV_d15n))) + geom_density()
cv_taxa |> ggplot(aes(sqrt(CV_d15n))) + geom_density()


shapiro.test(cv_taxa$CV_d15n)
shapiro.test(log(cv_taxa$CV_d15n))
shapiro.test(sqrt(cv_taxa$CV_d15n))
car::leveneTest(log(cv_taxa$CV_d15n), group = cv_taxa$taxon_code)

fit_taxa <- lm(log(CV_d15n) ~ taxon_code, data = cv_taxa)
anova(fit_taxa)
agricolae::HSD.test(fit_taxa,'taxon_code',group = FALSE,console = TRUE, unbalanced = TRUE) 

aov_taxa <- aov(log(CV_d15n)~taxon_code, data=cv_taxa)
summary(aov_taxa)
TukeyHSD(aov_taxa)


# FFGs

cv_ffg |> ggplot(aes(CV_d15n)) + geom_density()
cv_ffg |> ggplot(aes(log(CV_d15n))) + geom_density()
cv_ffg |> ggplot(aes(sqrt(CV_d15n))) + geom_density()

shapiro.test(cv_ffg$CV_d15n)
shapiro.test(log(cv_ffg$CV_d15n))
shapiro.test(sqrt(cv_ffg$CV_d15n))
car::leveneTest(log(cv_ffg$CV_d15n), group = cv_ffg$ffg)

fit_ffg <- lm(log(CV_d15n) ~ ffg, data = cv_ffg)
anova(fit_ffg)
agricolae::HSD.test(fit_ffg,'ffg',group = FALSE,console = TRUE, unbalanced = TRUE) 

aov_ffg <- aov(log(CV_d15n)~ffg, data=cv_ffg)
summary(aov_ffg)
TukeyHSD(aov_ffg)






## OLD CODE BELOW ##################


# #####check normality for FFG first####
# 
# anovaCV_dat_ffg <- invertdata %>% 
#   #filter(ffg %in% focusffg) %>%
#   filter(ffg!= "Shredder") %>% 
#   group_by(ffg,site_id) %>% 
#   summarise(n = n(),
#             mean_d15n = mean(d15N, na.rm = TRUE),
#             sd_d15n = sd(d15N, na.rm = TRUE),
#             CV_d15n = sd_d15n/mean_d15n) %>%
#   drop_na()
# 
# 
# anovaCV_dat_ffg %>% 
#   # filter(CV_d15n < 6) %>% 
#   # filter(ffg != "Shredder") %>% 
#   mutate(CV_d15n = CV_d15n^(.5)) %>% 
#   ggplot(aes(CV_d15n)) + 
#   geom_density()+
#   facet_wrap(~ffg)
# 
# anovaCV_dat_ffg <- anovaCV_dat_ffg %>% 
#   filter(CV_d15n < 6) %>%
#   filter(ffg != "Shredder")
# 
# fit_ffg <- lm(CV_d15n ~ ffg, data = anovaCV_dat_ffg) 
# MASS::boxcox(fit_ffg, lambda = seq(-2,2,1/100))
# 
# 
# shapiro.test((anovaCV_dat_ffg$CV_d15n^(0.5)))
# car::leveneTest((anovaCV_dat_ffg$CV_d15n^(0.5)), group = anovaCV_dat_ffg$ffg)
# 
# fit_ffg <- lm(sqrt(CV_d15n)~ffg, data = anovaCV_dat_ffg)
# 
# anova(fit_ffg)
# library(agricolae)
# HSD.test(fit_ffg,'ffg',group = FALSE,console = TRUE, unbalanced = TRUE) 
# 
# 
# ####check for normality by taxa####
# 
# anovaCV_dat <- invertdata %>%#make a new data frame using metadat
#   filter(taxon_code %in% focustaxon$taxon_code) %>%#filter taxon using the high distribution
#   group_by(taxon_code,site_id) %>%#group by site and taxon
#   summarise(n = n(),#summarise get sample size, mean, standard deviation and CV
#             mean_d15n = mean(d15N, na.rm = TRUE),
#             sd_d15n = sd(d15N, na.rm = TRUE),
#             CV_d15n = sd_d15n/mean_d15n) %>% 
#   drop_na() %>% 
#   filter(CV_d15n > 0 & CV_d15n < 1)
# 
# anovaCV_dat%>% 
#   filter(CV_d15n > 0) %>% 
#   mutate(CV_d15n = CV_d15n^(0.2)) %>% 
#   ggplot(aes(CV_d15n)) + 
#   geom_histogram()+
#   facet_wrap(~taxon_code)
# 
# fit <- lm(CV_d15n ~ taxon_code, data = anovaCV_dat)
# boxcox(fit, lambda = seq(-2,2,1/100))
# 
# shapiro.test(anovaCV_dat$CV_d15n^(.2))
# leveneTest(anovaCV_dat$CV_d15n^(.2), group = anovaCV_dat$taxon_code)
# 
# fit <- lm((CV_d15n^.2)~taxon_code, data = anovaCV_dat)
# anova(fit)
# 
# HSD.test(fit,'taxon_code',group = FALSE,console = TRUE, unbalanced = TRUE)




