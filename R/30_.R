
# Criteria #2: low within-site variation

## Prep

# libraries
library(tidyverse) 
library(here)
library(cowplot)

# functions
lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1-((1-conf_level)/2),n-1)*se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1-((1-conf_level)/2),n-1)*se
}


## Data -----------

# filter the data so only 75% distribution or higher is kept
focustaxon <- Distr_dat %>%
  filter(distr >= 75)

foucusffg <- Distr_dat_FFG %>% 
  filter(distr >= 75)

## Summarize isotope data ---------

## THIS NEEDS TO BE FIXED
## TAKING SUMMARY OF SUMMARY - should we just analyze with LMEs?

lowCV_dat <- invertdata %>%
  filter(taxon_code %in% focustaxon$taxon_code) %>%
  # summarize by site and taxon
  group_by(taxon_code, site_id) %>%
  summarise(
    n = n(),
    mean_d15n = mean(d15N, na.rm = TRUE),
    sd_d15n = sd(d15N, na.rm = TRUE),
    CV_d15n = sd_d15n/mean_d15n, 
    .groups = 'drop'
    ) %>% 
  drop_na() %>% # 112 to 90
  filter(CV_d15n > 0 & CV_d15n < 1) %>%  # 90 to 87
  # now summarize by taxon
  group_by(taxon_code) %>% 
  summarise(
    n = n(),
    mean_CV = mean(CV_d15n, na.rm = TRUE),
    sd_CV = sd(CV_d15n)
    ) %>% 
  # add error
  mutate(
    se = sd_CV/sqrt(n),
    low = lower_ci(mean_CV,se,n),
    high = upper_ci(mean_CV,se,n)
    )


lowCV_dat_FFG <- invertdata %>% 
  # remove shredderd
  filter(ffg!= "Shredder")%>% 
  # summarize by site and ffg
  group_by(ffg,site_id) %>% 
  summarise(
    n = n(),
    mean_d15n = mean(d15N, na.rm = TRUE),
    sd_d15n = sd(d15N, na.rm = TRUE),
    CV_d15n = sd_d15n/mean_d15n
    ) %>%
  filter(CV_d15n < 6) %>% 
  drop_na() %>% 
  # now summarize by ffgs
  group_by(ffg) %>% 
  summarise(
    n = n(),
    mean_CV = mean(CV_d15n, na.rm = TRUE),
    sd_CV = sd(CV_d15n)
    ) %>% 
  drop_na() %>% 
  # add error
  mutate(
    se = sd_CV/sqrt(n),
    low = lower_ci(mean_CV,se,n),
    high = upper_ci(mean_CV,se,n)
    )

## Plot -----------

CV_taxon <- lowCV_dat %>%
  ggplot(aes(x=reorder(taxon_code, mean_CV, median), y = mean_CV))+
  geom_pointrange(ymin = lowCV_dat$low, ymax = lowCV_dat$high,
                  size = 1, color = "black", fill = "#666666", shape = 21)+
  theme_classic()+
  coord_flip()+
  labs(x = "Taxon", y = expression('CV (' ~{delta}^15*N~')'))+  
  theme(axis.title  = element_blank(), 
        axis.text = element_text(size = 12, color = "black"))+  
  scale_y_continuous(limits = c(0,0.60), breaks = c(0.00,0.20,0.40,0.60))


CV_FFG <- lowCV_dat_FFG %>%
  ggplot(aes(x=reorder(ffg, mean_CV, median), y = mean_CV))+    
  geom_pointrange(ymin = lowCV_dat_FFG$low, ymax = lowCV_dat_FFG$high,
                  size = 1, color = "black", fill = "#666666", shape = 21)+                                                  #use a boxplot
  theme_classic()+
  coord_flip()+
  labs(x = "FFG", y = expression('Mean CV (' ~{delta}^15*N~')'))+
  theme(axis.title  = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title.y  = element_blank())+ 
  scale_y_continuous(limits = c(0,0.60), breaks = c(0.00,0.20,0.40,0.60))


CVfig_All <- cowplot::plot_grid(
  CV_taxon,CV_FFG, align = "hv", labels = c("A","B"),
  ncol = 1, nrow = 2, label_x = 0.160, label_y = 0.99,label_size = 18
  )


ggsave(here("out", "combined_CVplot.png"), 
          plot = CVfig_All, ragg::agg_png,
          units = "in", width = 8.5, height = 8.5)

ggsave(here("out", "fig2_combined_CVplot.pdf"), 
          plot = CVfig_All, device = cairo_pdf,
          units = "in", width = 8.5, height = 8.5)


## ANOVA of CVs ------------------------------

#####check normality for FFG first####

anovaCV_dat_ffg <- invertdata %>% 
  #filter(ffg %in% focusffg) %>%
  filter(ffg!= "Shredder") %>% 
  group_by(ffg,site_id) %>% 
  summarise(n = n(),
            mean_d15n = mean(d15N, na.rm = TRUE),
            sd_d15n = sd(d15N, na.rm = TRUE),
            CV_d15n = sd_d15n/mean_d15n) %>%
  drop_na()



anovaCV_dat_ffg %>% 
  filter(CV_d15n < 6) %>% 
  filter(ffg != "Shredder") %>% 
  mutate(CV_d15n = CV_d15n^(.5)) %>% 
  ggplot(aes(CV_d15n)) + 
  geom_density()+
  facet_wrap(~ffg)

anovaCV_dat_ffg <- anovaCV_dat_ffg %>% 
  filter(CV_d15n < 6) %>%
  filter(ffg != "Shredder")

fit_ffg <- lm(CV_d15n ~ ffg, data = anovaCV_dat_ffg) 
MASS::boxcox(fit_ffg, lambda = seq(-2,2,1/100))


shapiro.test((anovaCV_dat_ffg$CV_d15n^(0.5)))
car::leveneTest((anovaCV_dat_ffg$CV_d15n^(0.5)), group = anovaCV_dat_ffg$ffg)

fit_ffg <- lm(sqrt(CV_d15n)~ffg, data = anovaCV_dat_ffg)

anova(fit_ffg)
library(agricolae)
HSD.test(fit_ffg,'ffg',group = FALSE,console = TRUE, unbalanced = TRUE) 


####check for normality by taxa####

anovaCV_dat <- invertdata %>%#make a new data frame using metadat
  filter(taxon_code %in% focustaxon$taxon_code) %>%#filter taxon using the high distribution
  group_by(taxon_code,site_id) %>%#group by site and taxon
  summarise(n = n(),#summarise get sample size, mean, standard deviation and CV
            mean_d15n = mean(d15N, na.rm = TRUE),
            sd_d15n = sd(d15N, na.rm = TRUE),
            CV_d15n = sd_d15n/mean_d15n) %>% 
  drop_na() %>% 
  filter(CV_d15n > 0 & CV_d15n < 1)

anovaCV_dat%>% 
  filter(CV_d15n > 0) %>% 
  mutate(CV_d15n = CV_d15n^(0.2)) %>% 
  ggplot(aes(CV_d15n)) + 
  geom_histogram()+
  facet_wrap(~taxon_code)

fit <- lm(CV_d15n ~ taxon_code, data = anovaCV_dat)
boxcox(fit, lambda = seq(-2,2,1/100))

shapiro.test(anovaCV_dat$CV_d15n^(.2))
leveneTest(anovaCV_dat$CV_d15n^(.2), group = anovaCV_dat$taxon_code)

fit <- lm((CV_d15n^.2)~taxon_code, data = anovaCV_dat)
anova(fit)

HSD.test(fit,'taxon_code',group = FALSE,console = TRUE, unbalanced = TRUE)




