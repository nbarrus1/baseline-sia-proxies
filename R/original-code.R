#Start

rm(list = ls())

#*******************************************************************************************************************************
#-------------------------------------------------------------------------------------------------------------------------------
#####libraries#####
#-------------------------------------------------------------------------------------------------------------------------------
#*******************************************************************************************************************************

library(tidyverse)                                                 #ggplot2, tibble, dplyr, tibble, forcats,purrr,stringr,tidyr
library(ggrepel)                                                   #ggrepel
library(gridExtra)                                                 #grid extra
library(broom)                                                     #broom
library(lemon)
library(MASS)
library(ggpubr)
library(car)
library(cowplot)

#************************************************************************************************
#------------------------------------------------------------------------------------------------
#####Load, clean, and combine data####
#------------------------------------------------------------------------------------------------
#************************************************************************************************

#----------------------------------------------------------
####load data####
#----------------------------------------------------------

invert_group <- read_csv("metadata_bugs.csv") %>% 
  dplyr::select(taxon, ffg) %>% 
  rename(taxon_code = taxon)

isodat <- read_csv("sia_inverts_compiled.csv") 

PCAdat <- read_csv("data_PCA_results.csv") %>% 
  mutate(PC1 = (PC1 + 5))#read in file data_PCA_results.csv

LandUsedat <- read_csv("land-use.csv")                            #read in file land-use.csv

samplelog <- read_csv("sample_log_bugs.csv")


PCAdat <- PCAdat %>% 
  left_join(LandUsedat, by = "site_id")

PCAdat %>% 
  ggplot(aes(x = PC1, y = lulc_nat*100))+
  geom_point(size = 3, shape = 21, color = "black", fill = "#666666")+
  geom_smooth(method = "lm", color = "black")+
  theme_classic()+
  labs(y = "Percent Natural Landcover",
       x = "PC1 (Environmental Gradient)")+
  geom_text(aes(label = site_id), vjust = 0, nudge_y = 1.5)+
  geom_text(aes(label = "r = - 0.700", x = 2.5,y = 25))+
  geom_text(aes(label = "p = 0.002", x = 2.5,y = 21))

cor.test(x = PCAdat$PC1, y = PCAdat$lulc_nat)
summary(lm(PCAdat$lulc_nat~PCAdat$PC1))

#----------------------------------------------------------
####Combine Data Clean Data####
#----------------------------------------------------------

alldat <- isodat %>% 
  left_join(PCAdat, by = "site_id") %>% 
  left_join(LandUsedat, by = "site_id") %>% 
  left_join(invert_group, by = "taxon_code") %>% 
  dplyr::select(sample_year, stream_name, site_id,
         resource, taxon_code, d15N, PC1, ffg) %>% 
  filter(site_id != "LR01") %>% 
  filter(sample_year == "2016") %>%
  filter(resource == "biofilm" | resource == "detritus" |
           resource == "fish" | resource == "invert") 
  
invertdata <- alldat %>% 
  filter(resource == "invert") %>% 
  group_by(taxon_code) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5)

primprodat <- alldat %>% 
  filter(taxon_code == "Biofilm" | taxon_code == "Seston") %>% 
  group_by(taxon_code) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 5) %>% 
  group_by(taxon_code,site_id,PC1) %>% 
  summarise(mean_d15N = mean(d15N))


#************************************************************************************************
#------------------------------------------------------------------------------------------------
######Analysis#####
#------------------------------------------------------------------------------------------------
#************************************************************************************************

#----------------------------------------------------------
#####Criteria 1: Distribution #####
#----------------------------------------------------------


focustax <- as.factor(invertdata$taxon_code) %>% levels()                   #obtain taxon to filter the large dataset
study_sites <- as.factor(invertdata$site_id) %>% levels()            #obtain Sites id from data
n_sites = length(study_sites)                                       #total number of sites


focusffg <- as.factor(invertdata$ffg) %>% levels()                   #obtain taxon to filter the large dataset


####calculate distribution####

Distr_dat <- invertdata %>%                                          #make new dataframe with dattaxall
  filter(sample_year == "2016") %>% 
  filter(site_id != "LR01") %>% 
  group_by(taxon_code, site_id) %>%                                 #group everything by taxon and site_id  
  summarise(n = n()) %>%                                            #using the groupings find the number of observations 
  mutate(presence = 1) %>%                                          #make new variable for presence 1
  group_by(taxon_code) %>%                                          #group again by taxon 
  summarise(distr = (sum(presence)/n_sites)*100) %>% 
  filter(taxon_code %in% focustax)                                  #Filter all the taxon and use focustax
  #add all values of presence by the group taxon divide total number of sites multiply by 100 to get percentage

Distr_dat_FFG <- invertdata %>%                                          #make new dataframe with dattaxall
  filter(taxon_code %in% focustax) %>% 
  group_by(ffg, site_id) %>%                                      #group everything by taxon and site_id  
  summarise(n = n()) %>%                                            #using the groupings find the number of observations 
  mutate(presence = 1) %>%                                          #make new variable for presence 1
  group_by(ffg) %>%                                               #group again by taxon 
  summarise(distr = (sum(presence)/n_sites)*100)                     #add all values of presence by the group taxon divide total number of sites multiply by 100 to get percentage
      
####plot distribution####

DistFigTaxon <- Distr_dat %>%                                                       #Take dist_r
  filter(taxon_code %in% focustax) %>%                                   #Filter all the taxon and use focustax
  ggplot( aes(x=fct_reorder(taxon_code, distr, median), y=distr))+       #plot using the reorder taxon based on the median distribution, by distribution
  geom_point(size = 3,color = "black", shape = 21, fill = "#666666")+          #Plot point, size of three, make the outline black, using shape circle and fill of grey color
  theme_classic()+                                                  #use the classic theme
  #coord_flip()+                                                     #flip the coordinates
  geom_hline(yintercept = 75, color = "red", size =1, linetype = "dashed")+    #create a size 1 dashed red horizontal line at 75% 
  theme(axis.text = element_text(size = 12, face = "bold", color = "black"),   #take axis text make it size 12 and bold it
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        title = element_text(size = 16, face = "bold", color = "black"),
        axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.text.x.bottom = element_text(size = 10, face = "bold", color = "black"))+ #take axis title make it size 16 and bold
  labs(x = "Taxon", y = "Percent of Sites")+                              #label  x axis Taxon and y % of sites
  theme(axis.text.y = element_text(size = rel(0.8)))+               #take the y axis text and rescale so it is a bit smaller
  scale_y_continuous(limits = c(0,100))                             #scale y so it starts at 0 and goes to 100

####plot FFG distribution data####

DistFigFFG <- Distr_dat_FFG %>% 
  ggplot(aes(x = fct_reorder(ffg, distr, median), y = distr))+
  geom_point(size = 3, shape = 21, color = "black",fill = "#666666")+
  theme_classic()+
  #coord_flip()+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10, face = "bold", color = "black"),
        title = element_text(size = 16, face = "bold", color = "black"),
        axis.text.x = element_text(angle = 90))+
  labs(x = "FFG", y = "Percent of Sites")+
  geom_hline(yintercept = 75, color = "red", linetype = "dashed", size = 1)+
  scale_y_continuous(limits = c(0,100))

#####plot together#####

distFigAll <-plot_grid(DistFigTaxon, DistFigFFG, labels = c("A","B"), align = "hv",
          label_size = 18, ncol = 2, nrow = 1, label_x = 0.15, label_y = 0.99)

####save plot####

save_plot(filename = "figures/combined_distribution.png", plot = distFigAll, device = png(),
       units = "in", base_width = 8.5, base_height = 6)

#------------------------------------------------------------
######Criteria 2: Low Within Site Variation ######
#------------------------------------------------------------

#####create focus taxon#####

focustaxon <- Distr_dat %>%                                        #create new data frame using Distr_dat
  filter(distr >= 75)                                              #filter the data so only 75% distribution or higher is kept

foucusffg <- Distr_dat_FFG %>% 
  filter(distr >= 75)

#####summarise isotope data####

lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1-((1-conf_level)/2),n-1)*se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1-((1-conf_level)/2),n-1)*se
}

lowCV_dat <- invertdata %>%                                           #make a new data frame using metadat
  filter(taxon_code %in% focustaxon$taxon_code) %>%                          #filter taxon using the high distribution
  group_by(taxon_code,site_id) %>%                                      #group by site and taxon
  summarise(n = n(),                                               #summarise get sample size, mean, standard deviation and CV
            mean_d15n = mean(d15N, na.rm = TRUE),
            sd_d15n = sd(d15N, na.rm = TRUE),
            CV_d15n = sd_d15n/mean_d15n) %>% 
  drop_na() %>% 
  filter(CV_d15n > 0 & CV_d15n < 1) %>% 
  group_by(taxon_code) %>% 
  summarise(n = n(),
            mean_CV = mean(CV_d15n, na.rm = TRUE),
            sd_CV = sd(CV_d15n)) %>% 
  drop_na() %>% 
  mutate(se = sd_CV/sqrt(n),
         low = lower_ci(mean_CV,se,n),
         high = upper_ci(mean_CV,se,n))

lowCV_dat_FFG <- invertdata %>% 
  filter(ffg!= "Shredder") %>% 
  group_by(ffg,site_id) %>% 
  summarise(n = n(),
            mean_d15n = mean(d15N, na.rm = TRUE),
            sd_d15n = sd(d15N, na.rm = TRUE),
            CV_d15n = sd_d15n/mean_d15n) %>%
  filter(CV_d15n < 6) %>% 
  drop_na() %>% 
  group_by(ffg) %>% 
  summarise(n = n(),
            mean_CV = mean(CV_d15n, na.rm = TRUE),
            sd_CV = sd(CV_d15n)) %>% 
  drop_na() %>% 
  mutate(se = sd_CV/sqrt(n),
         low = lower_ci(mean_CV,se,n),
         high = upper_ci(mean_CV,se,n))

#####Plot Summarised data#####

CV_taxon <- lowCV_dat %>%                                                      #use lowCV_dat
  ggplot(aes(x=reorder(taxon_code, mean_CV, median), y = mean_CV))+     #plot using reordered taxon with CV as the order vs CV                                               #use a boxplot
  geom_pointrange(ymin = lowCV_dat$low, ymax = lowCV_dat$high,
                  size = 1, color = "black", fill = "#666666", shape = 21)+
  theme_classic()+                                                 #themeclassic
  coord_flip()+                                                    #flip coordinates
  labs(x = "Taxon", y = expression('CV (' ~{delta}^15*N~')'))+     #label x taxon, and y CV with delta15N notation
  theme(axis.title  = element_blank(), #change axis titles to bold 16 black
        axis.text = element_text(size = 12, color = "black"))+  #change axis text to bold 12 black
  scale_y_continuous(limits = c(0,0.60), breaks = c(0.00,0.20,0.40,0.60))
  
#####Plot CV and FFG####

CV_FFG <- lowCV_dat_FFG %>%                                                      #use lowCV_dat
  ggplot(aes(x=reorder(ffg, mean_CV, median), y = mean_CV))+     #plot using reordered taxon with CV as the order vs CV
  geom_pointrange(ymin = lowCV_dat_FFG$low, ymax = lowCV_dat_FFG$high,
                  size = 1, color = "black", fill = "#666666", shape = 21)+                                                  #use a boxplot
  theme_classic()+                                                 #themeclassic
  coord_flip()+                                                    #flip coordinates
  labs(x = "FFG", y = expression('Mean CV (' ~{delta}^15*N~')'))+     #label x taxon, and y CV with delta15N notation
  theme(axis.title  = element_text(size = 16, color = "black"), #change axis titles to bold 16 black
        axis.text = element_text(size = 12, color = "black"),
        axis.title.y  = element_blank())+  #change axis text to bold 12 black
  scale_y_continuous(limits = c(0,0.60), breaks = c(0.00,0.20,0.40,0.60))

#####Plot together#####

CVfig_All <- plot_grid(CV_taxon,CV_FFG, align = "hv", labels = c("A","B"),
          ncol = 1, nrow = 2, label_x = 0.160, label_y = 0.99,
          label_size = 18)

#####save plot####

save_plot(filename = "figures/combined_CVplot.png", plot = CVfig_All, device = png(),
          units = "in", base_width = 8.5, base_height = 8.5)

save_plot(filename = "figures/pdfs/fig2_combined_CVplot.pdf", plot = CVfig_All, device = pdf(),
          units = "in", base_width = 8.5, base_height = 8.5)
#####ANOVA of CVs####

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
boxcox(fit_ffg, lambda = seq(-2,2,1/100))


shapiro.test((anovaCV_dat_ffg$CV_d15n^(0.5)))
leveneTest((anovaCV_dat_ffg$CV_d15n^(0.5)), group = anovaCV_dat_ffg$ffg)

fit_ffg <- lm(sqrt(CV_d15n)~ffg, data = anovaCV_dat_ffg)

anova(fit_ffg)
library(agricolae)
HSD.test(fit_ffg,'ffg',group = FALSE,console = TRUE, unbalanced = TRUE) 


####check for normality by taxa####

anovaCV_dat <- invertdata %>%                                           #make a new data frame using metadat
  filter(taxon_code %in% focustaxon$taxon_code) %>%                          #filter taxon using the high distribution
  group_by(taxon_code,site_id) %>%                                      #group by site and taxon
  summarise(n = n(),                                               #summarise get sample size, mean, standard deviation and CV
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

#------------------------------------------------------------------
#####Criteria 3#####
#------------------------------------------------------------------

#####Models of baseline, basal resources, and ffg using nest#####

pc1dat <- invertdata %>% 
  group_by(taxon_code, PC1) %>% 
  summarise(mean_d15n = mean(d15N))

pc1dat_ffg <- invertdata %>% 
  group_by(ffg, PC1) %>% 
  summarise(mean_d15n = mean(d15N))

pc1dat_pp <- primprodat %>% 
  group_by(taxon_code, PC1) %>% 
  summarise(mean_d15N = mean(mean_d15N))

###using pc1
pc1_corr <- pc1dat %>% 
  dplyr::select(taxon_code, PC1, mean_d15n) %>% 
  group_by(taxon_code) %>% 
  nest()

pc1_corr_ffg <- pc1dat_ffg %>% 
  dplyr::select(ffg, PC1, mean_d15n) %>% 
  group_by(ffg) %>% 
  nest()

pc1_corr_pp <- pc1dat_pp %>% 
  rename(mean_d15n = mean_d15N) %>% 
  dplyr::select(taxon_code, PC1, mean_d15n) %>% 
  group_by(taxon_code) %>% 
  nest()

#create fuction

lm_PC1 <- function(df) {lm(mean_d15n ~ PC1, data = df)}


###for PC1
corr_reg_pc1 <- pc1_corr %>% 
  mutate(
    fit_pc1 = map(data, lm_PC1),
    tidy_N = map(fit_pc1, tidy),
    glance_N = map(fit_pc1, glance),
    augment_N = map(fit_pc1, augment), 
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2)),
    conf_tidy = map(conf_int,tidy)
  )

corr_reg_pc1_ffg <- pc1_corr_ffg %>% 
  mutate(
    fit_pc1 = map(data, lm_PC1),
    tidy_N = map(fit_pc1, tidy),
    glance_N = map(fit_pc1, glance),
    augment_N = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2)),
    conf_tidy = map(conf_int,tidy)
  )

corr_reg_pc1_pp <- pc1_corr_pp %>% 
  mutate(
    fit_pc1 = map(data, lm_PC1),
    tidy_N = map(fit_pc1, tidy),
    glance_N = map(fit_pc1, glance),
    augment_N = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2)),
    conf_tidy = map(conf_int,tidy)
  )


#####parameter info####

#### for pc1

pc1_parainfo <- corr_reg_pc1 %>% 
  unnest(tidy_N, .drop = TRUE) %>% 
  rename(Baseline = taxon_code) %>% 
  mutate(group = "taxa") %>% 
  dplyr::select(-data, -fit_pc1, -glance_N, -augment_N)

pc1_parainfo_ffg <- corr_reg_pc1_ffg %>% 
  unnest(tidy_N, .drop = TRUE) %>% 
  rename(Baseline = ffg) %>% 
  mutate(group = "ffg") %>% 
  dplyr::select(-data, -fit_pc1, -glance_N, -augment_N)

pc1_parainfo_pp <- corr_reg_pc1_pp %>% 
  unnest(tidy_N, .drop = TRUE) %>% 
  rename(Baseline = taxon_code) %>% 
  mutate(group = "basal") %>% 
  dplyr::select(-data, -fit_pc1, -glance_N, -augment_N)


###combine pc1 parameter info

pc1_parainfo_all <- pc1_parainfo %>% 
  bind_rows(pc1_parainfo_ffg) %>% 
  bind_rows(pc1_parainfo_pp)

pc1_parainfo_all_int <- pc1_parainfo_all %>% 
  filter(term == "(Intercept)") %>% 
  rename(intercept = estimate)

pc1_parainfo_all_slp <- pc1_parainfo_all %>% 
  filter(term == "PC1") %>% 
  rename(slope = estimate)


###best model info using pc1

pc1_bestModel_info <- corr_reg_pc1 %>%                                      #place AIC, R2 and BIc in object
  unnest(glance_N, .drop = TRUE) %>% 
  mutate(group = "taxa") %>% 
  rename(Baseline = taxon_code) %>% 
  dplyr::select(-data, -fit_pc1, -tidy_N, -augment_N)

pc1_bestModel_info_ffg <- corr_reg_pc1_ffg %>%                                      #place AIC, R2 and BIc in object
  unnest(glance_N, .drop = TRUE) %>% 
  mutate(group = "ffg") %>% 
  rename(Baseline = ffg) %>% 
  dplyr::select(-data, -fit_pc1, -tidy_N, -augment_N)

pc1_bestModel_info_pp <- corr_reg_pc1_pp %>%                                      #place AIC, R2 and BIc in object
  unnest(glance_N, .drop = TRUE) %>% 
  mutate(group = "basal") %>% 
  rename(Baseline = taxon_code) %>% 
  dplyr::select(-data, -fit_pc1, -tidy_N, -augment_N)



###combe all of best model info using pc1

pc1_bestModel_info_all <- pc1_bestModel_info %>% 
  bind_rows(pc1_bestModel_info_ffg) %>% 
  bind_rows(pc1_bestModel_info_pp)

####slope confidence interval

conf_int_pc1 <- corr_reg_pc1 %>% 
  unnest(conf_tidy) 

conf_int_pc1 <- conf_int_pc1 %>% 
  dplyr::select(taxon_code, x) %>% 
  rename(Baseline = taxon_code) %>% 
  mutate(group = "taxa")

conf_int_pc1_ffg <- corr_reg_pc1_ffg %>% 
  unnest(conf_tidy)

conf_int_pc1_ffg <- conf_int_pc1_ffg %>% 
  dplyr::select(ffg, x) %>% 
  rename(Baseline = ffg) %>% 
  mutate(group = "ffg")

conf_int_pc1_pp <- corr_reg_pc1_pp %>% 
  unnest(conf_tidy) 

conf_int_pc1_pp <- conf_int_pc1_pp %>% 
  dplyr::select(taxon_code, x) %>% 
  rename(Baseline = taxon_code) %>% 
  mutate(group = "basal")

###combine all together

conf_int_pc1_all <- conf_int_pc1 %>% 
  bind_rows(conf_int_pc1_ffg) %>% 
  bind_rows(conf_int_pc1_pp)


######SUPPLEMENTAL TABLE CRITERIA 3 Regression Stats #####

Crit3_Table_pc1 <- pc1_bestModel_info_all %>% 
  dplyr::select(Baseline, p.value, adj.r.squared, group) %>% 
  left_join(conf_int_pc1_all, by = c("Baseline", "group")) %>% 
  left_join(pc1_parainfo_all_int, by = c("Baseline", "group")) %>% 
  left_join(pc1_parainfo_all_slp, by = c("Baseline", "group")) %>% 
  dplyr::select(Baseline, intercept, slope, x, p.value, adj.r.squared, group) %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p.value = if_else(condition = p.value == 0, true = paste("<0.001"), false = paste(p.value)))


#####Adjusted R squared plots of Baselines and Basal Resources####
sig_taxa_pc1 <- Crit3_Table_pc1 %>% 
  filter(group == "taxa") %>% 
  rename(taxon_code = Baseline) %>% 
  left_join(pc1dat, by = "taxon_code") %>% 
  filter(p.value <= 0.05) %>% 
  filter(taxon_code %in% focustaxon$taxon_code)

sig_FFG_pc1 <- Crit3_Table_pc1 %>% 
  filter(group== "ffg") %>% 
  rename(ffg = Baseline) %>% 
  left_join(pc1dat_ffg, by = "ffg") %>% 
  filter(p.value <= 0.05) %>% 
  filter(ffg != "Shredder")

sig_base_pc1 <- Crit3_Table_pc1 %>% 
  filter(group == "basal") %>% 
  rename(taxon_code = Baseline ) %>% 
  left_join(pc1dat_pp, by = "taxon_code") %>% 
  filter(p.value <= 0.05) %>% 
  rename(mean_d15N = mean_d15N)



#####plot of pc1 use and d15N with facet wrap of taxon####

#
format_pval <- function(pval){
  pval <- scales::pvalue(pval, accuracy= 0.001, add_p = TRUE)
  gsub(pattern = "(=|<)", replacement = " \\1 ", x = pval)
}
##pc1

pc1vsd15n <- pc1dat %>%                                                           #use correlation dat
  filter(taxon_code %in% focustaxon$taxon_code) %>%
  ggplot(aes(x=PC1, y = mean_d15n))+                                         #plot usein d15n as a function of lulc_nat
  geom_point()+                                                              #plot points
  geom_smooth(data = sig_taxa_pc1, method = lm, color = "black", fill = "#333333")+                                                  #plot best fit line using linear model
  facet_rep_wrap(~taxon_code)+                                                        #separate by taxon
  theme_classic()+                                                           #use theme classic
  labs(x = "Longitudinal Gradient (PC1)", y =  expression({delta}^15*N~'\u2030'),
       title = "C) Primary Consumers (Taxa)")+ #give labels expression tells it to use delta format
  theme(axis.title = element_text( size = 12, color = "black"),     #axis titles use bold 16 black 
        axis.text = element_text( size = 8, color = "black"),
        strip.text = element_text( size = 10, color = "black"),
        title = element_text( size = 12, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0))+
  scale_y_continuous(limits = c(-1,12), breaks = c(0,4,8,12))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 12, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)


#####Plot of pc1 and d15N with ffg as facet####
pc1dat_ffg.tmp <- pc1dat_ffg %>% 
  rename(taxon_code = ffg)

pc1dat_pp.tmp <- pc1dat_pp %>% 
  rename(mean_d15n = mean_d15N)

PC1dat_all <- pc1dat %>% 
  bind_rows(pc1dat_ffg.tmp) %>% 
  bind_rows(pc1dat_pp.tmp)

write_csv(PC1dat_all, "Crit3_dat.csv")

###using pc1
pc1dat_ffg <- pc1dat_ffg %>% 
  filter(ffg != "Shredder")

pc1vsd15n_ffg <- pc1dat_ffg %>%
  ggplot(aes(x=PC1, y = mean_d15n))+                                         #plot usein d15n as a function of lulc_nat
  geom_point()+                                                              #plot points
  geom_smooth(data = sig_FFG_pc1, method = lm, color = "black", fill = "#333333")+                                                  #plot best fit line using linear model
  facet_rep_wrap(~ffg)+                                                        #separate by taxon
  theme_classic()+                                                           #use theme classic
  labs(x = "Longitudinal Gradient (PC1)", y =  expression({delta}^15*N~'\u2030'),
       title = "B) Primary Consumers (FFG)")+ #give labels expression tells it to use delta format
  theme(axis.title = element_text( size = 12, color = "black"),     #axis titles use bold 16 black 
        axis.text = element_text(size = 8, color = "black"),
        strip.text = element_text( size = 10, color = "black"),
        title = element_text( size = 12, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0))+
  scale_y_continuous(limits = c(-1,12), breaks = c(0,4,8,12))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 12, aes(label = paste(..rr.label.., format_pval(..p..),
                                           sep = "*`,`~")), size = 2)


###plot using PC1

PC1vsd15Nprim <- pc1dat_pp %>% 
  ggplot( aes(x = PC1, y = mean_d15N))+
  geom_point()+
  geom_smooth(data = sig_base_pc1, method = lm, color = "black", fill = "#333333")+
  facet_rep_wrap(~taxon_code)+
  theme_classic()+
  labs(x = "Longitudinal Gradient (PC1)", y = expression({delta}^15*N~'\u2030'),
       title = "A) Basal Resources")+
  theme(axis.title = element_text(size = 12, color = "black"),
        axis.text = element_text( size = 8, color = "black"),
        strip.text = element_text( size = 10, color = "black"),
        title = element_text( size = 12, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0))+
  scale_y_continuous(limits = c(-1,8), breaks = c(0,2,4,6,8))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 8, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)



#####All Plots together #####


#pc1

pc1vsd15n_all <- grid.arrange(PC1vsd15Nprim,pc1vsd15n_ffg,pc1vsd15n, ncol = 9, nrow = 5,
                              layout_matrix = rbind(c(1,1,1,1,1,2,2,2,2),c(1,1,1,1,1,2,2,2,2),
                                                    c(3,3,3,3,3,3,3,3,3),c(3,3,3,3,3,3,3,3,3),
                                                    c(3,3,3,3,3,3,3,3,3)))

save_plot(filename = "figures/pc1vsd15N_all.png", plot = pc1vsd15n_all, device = png(),
          units = "in", base_width = 8.5, base_height = 9.5)

save_plot(filename = "figures/pdfs/fig3_pc1vsd15N_all.pdf", plot = pc1vsd15n_all, device = pdf(),
          units = "in", base_width = 8.5, base_height = 9.5)

#####plot function R2 values#####


#for pc1

pc1_bestModel_info_all <- pc1_bestModel_info_all %>% #use this data
  drop_na() %>% 
  mutate(sig = if_else(condition = p.value <= 0.05, true = "significant", false = "nonsignificant", missing = NULL)) %>% 
  filter(adj.r.squared <= 0.9) %>% 
  rename(model = Baseline)

pc1_bestModel_info_all %>% 
  ggplot(aes(x=reorder(model,adj.r.squared, median), y = adj.r.squared))+         #adj r sqr vs reorderd taxon   
  geom_point(aes(color = sig),size = 3)+              #plot points with shape 21, color balck, grey fill, size 3
  theme_classic()+                                                                #use theme classic
  facet_grid(~group, scales = "free", space = "free_x")+
  labs(x = "Model", y = "Adjusted R squared")+                     #rename labels
  theme(axis.title = element_text(size = 14, face = "bold", color = "black"),     #adjust axis titles
        axis.text = element_text(size = 10, face ="bold", color ="black"),        #adjust axis text
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, face = "bold", color = "black", angle = 90))+                                  #adjust x axist text so it is vertical not horzontal
  scale_shape_discrete(breaks = c("A","B","C"),labels = c("Taxa","FFG", "Basal Resources"))+    #add red line at 0.6 that is dashed and size 1
  scale_y_continuous(limits = c(-0.2, 1.0))+
  scale_color_manual(values = c("midnightblue", "darkred")) 

ggsave(filename = "figures/adjRsqrvsmodel_pc1.png", plot = last_plot(), device = png(),
       units = "in", width = 8, height = 6)

#----------------------------------------------------------------
#######Criteria 4: Trophic Estimates independent of Land Use#####
#----------------------------------------------------------------

######Retrieve Baseline data#####

meanNdatall <- invertdata %>% 
  group_by(taxon_code,site_id) %>%                                      #group by site and taxon
  summarise(                                              #summarise get sample size, mean, standard deviation and CV
            mean_d15n = mean(d15N, na.rm = TRUE)) %>% 
  spread(taxon_code,mean_d15n)

ffg_meanN <- invertdata %>% 
  group_by(ffg,site_id) %>% 
  summarise(mean_d15n = mean(d15N, na.rm = TRUE)) %>% 
  spread(ffg, mean_d15n)

#####Trophic Position Estimate #####
fishdat <- isodat %>% 
  filter(resource == "fish") %>% 
  filter(taxon_code == "CKC"|
           taxon_code == "WHS"|
           taxon_code == "BNT") 

uncorrected_data<- fishdat %>% 
  left_join(meanNdatall,by = "site_id") %>% 
  left_join(ffg_meanN, by = "site_id") %>% 
  left_join(PCAdat, by = "site_id") %>% 
  mutate(discFactor = ((-0.281*d15N)+5.879)) %>% 
  dplyr::select(-sample_year, -d13C, -stream_name, -unique_id, -sia_sample_id, -sample_hitch,
         -compartment, -resource, -length_mm, -d13C_scl, -d15N_scl, -CN, -lulc_del, 
         -lulc_ag, -PC2) %>% 
  gather(key = "model", value = "correction", 4:30) %>% 
  drop_na(correction) %>% 
  mutate(TP_nobase = d15N/discFactor) 

bnt_uncorrected_data <- uncorrected_data %>% 
  filter(model == "Predator" & taxon_code == "BNT") %>% 
  dplyr::select(PC1,TP_nobase,site_id,taxon_code)
ckc_uncorrected_data <- uncorrected_data %>% 
  filter(model == "Predator" & taxon_code == "CKC")%>% 
  dplyr::select(PC1,TP_nobase,site_id,taxon_code)
whs_uncorrected_data <- uncorrected_data %>% 
  filter(model == "Predator" & taxon_code == "WHS")%>% 
  dplyr::select(PC1,TP_nobase,site_id,taxon_code)

comb_uncorrected_data <- bnt_uncorrected_data %>% 
  bind_rows(ckc_uncorrected_data) %>% 
  bind_rows(whs_uncorrected_data)

unique(uncorrected_data$model)
unique(uncorrected_data$taxon_code) 
unique(uncorrected_data$PC1)

ucor_data <- expand.grid(taxon_code = unique(uncorrected_data$taxon_code),
                         model = unique(uncorrected_data$model),
                         PC1 = unique(uncorrected_data$PC1))
ucor_data <- ucor_data %>% 
  left_join(comb_uncorrected_data, by = c("taxon_code", "PC1"))


corrected_data<- fishdat %>% 
  left_join(meanNdatall,by = "site_id") %>% 
  left_join(ffg_meanN, by = "site_id") %>% 
  left_join(PCAdat, by = "site_id") %>% 
  mutate(discFactor = ((-0.281*d15N)+5.879)) %>% 
  dplyr::select(-sample_year, -d13C, -stream_name, -unique_id, -sia_sample_id, -sample_hitch,
                -compartment, -resource, -length_mm, -d13C_scl, -d15N_scl, -CN, -lulc_del, 
                -lulc_ag, -PC2) %>% 
  gather(key = "model", value = "correction", 4:30) %>% 
  drop_na(correction) %>% 
  mutate(TP = (((d15N-correction)/discFactor)+2)) 

fishdataall <- corrected_data %>% 
  bind_rows(ucor_data)

#####uncorrected TP estimate Regressions####
###regression for BNT noncorrected###

fit_uncorrect <- lm(TP_nobase ~ PC1, data = bnt_uncorrected_data) 
summary(fit_uncorrect)
#regression equation BNT is TP = 0.5057*PC1 + 0.8598, Adj R square = 0.2268, pval < 0.001

###regression for CKC noncorrected###

fit_uncorrect <- lm(TP_nobase ~ PC1, data = ckc_uncorrected_data) 
summary(fit_uncorrect)
#regression equation CKC is TP = 0.2629*PC1 + 2.0782, Adj R square = 0.02, pval < 0.002

###regression for WHS noncorrected###

fit_uncorrect <- lm(TP_nobase ~ PC1, data = whs_uncorrected_data) 
summary(fit_uncorrect)
#regression equation WHS is TP = 0.2372*PC1 + 1.3266, Adj R square = 0.1034, pval < 0.001

write_csv(fishdataall, "fishdataall.csv")

uncor_dat <- fishdataall %>% 
  dplyr::select(PC1,TP_nobase,taxon_code)
  
uncor_b1_BNT <- uncor_b1 %>% 
  filter(taxon_code == "BNT")
uncor_b1_CKC <- uncor_b1 %>% 
  filter(taxon_code == "CKC")
uncor_b1_WHS <- uncor_b1 %>% 
  filter(taxon_code == "WHS")

fit_BNT <- lm(TP_nobase~PC1, data = uncor_b1_BNT)
summary(fit_BNT)
#b1 = 0.497
fit_CKC <- lm(TP_nobase~PC1, data = uncor_b1_CKC)
summary(fit_CKC)
#b1 = 0.186
fit_WHS <- lm(TP_nobase~PC1, data = uncor_b1_WHS)
summary(fit_WHS)
#b1 = 0.230

uncor_b1 <- tibble(taxon_code = c("BNT", "CKC", "WHS"),
                   uncor_b1 = c("0.497", "0.186", "0.230"))

uncor_b1$uncor_b1 <- as.numeric(uncor_b1$uncor_b1)

######Plot Trophic Position vs. gradient ######

###plots them all together hard to see any relationship


#pc1

TPplot_pc1 <- fishdataall %>%                                             #Use Trophic position dat
  ggplot(aes(x = PC1, y = TP_nobase))+                       #plot TP vs LUlc
  geom_point(data = uncorrected_data, aes(x = PC1, y = TP_nobase), color = "grey")+                                                    #plot points
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(data = uncorrected_data, aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+                                        #Plot linear model
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  facet_grid(model~taxon_code)+                                         #separately based on taxon_code
  theme_classic()+                                                 #use theme clasic
  labs( y = "Trophic Position", x = "Longitudinal Gradient (PC1)", title = "TP Estimates By Baseline")+                                       #Labels 
  theme(axis.title = element_text(size = 10, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        title = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0))   

ggsave(filename = "figures/tpvsnatall_pc1.pdf", plot = last_plot(), 
       device = "pdf", width = 52, height = 52, units = "in", limitsize = FALSE)


####BNT focus taxon

#pc1

BNT_TPplot_pc1 <- fishdataall %>% 
  filter(model %in% focustaxon$taxon_code) %>% 
  filter(taxon_code == "BNT") %>% 
  ggplot(aes(x = PC1, y = TP_nobase))+                       #plot TP vs LUlc
  geom_point(aes(x = PC1, y = TP_nobase), color = "grey")+                                                    #plot points
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+                                        #Plot linear model
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  facet_rep_wrap(~model)+                                         #separately based on taxon_code
  theme_classic()+                                                 #use theme clasic
  labs( y = "Trophic Position", x = "Longitudinal Gradient (PC1)")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 9, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)


ggsave(filename = "figures/pdfs/BNT_tpvsnatfocustax_pc1.pdf", plot = last_plot(), 
       device = "pdf", width = 8, height = 8, units = "in", limitsize = FALSE)

####CKC focus taxon

#pc1

CKC_TPplot_pc1 <- fishdataall %>% 
  filter(model %in% focustaxon$taxon_code) %>% 
  filter(taxon_code == "CKC") %>% 
  ggplot(aes(x = PC1, y = TP_nobase))+                       #plot TP vs LUlc
  geom_point(aes(x = PC1, y = TP_nobase), color = "grey")+                                                    #plot points
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+                                        #Plot linear model
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  facet_rep_wrap(~model)+                                         #separately based on taxon_code
  theme_classic()+                                                 #use theme clasic
  labs( y = "Trophic Position", x = "Longitudinal Gradient (PC1)")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 9, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)


ggsave(filename = "figures/pdfs/CKC_tpvsnatfocustax_pc1.pdf", plot = last_plot(), 
       device = "pdf", width = 8, height = 8, units = "in", limitsize = FALSE)

####WHS focus taxon

#pc1

WHS_TPplot_pc1 <- fishdataall %>% 
  filter(model %in% focustaxon$taxon_code) %>% 
  filter(taxon_code == "WHS") %>% 
  ggplot(aes(x = PC1, y = TP_nobase))+                       #plot TP vs LUlc
  geom_point(aes(x = PC1, y = TP_nobase), color = "grey")+                                                    #plot points
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+                                        #Plot linear model
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  facet_rep_wrap(~model)+                                         #separately based on taxon_code
  theme_classic()+                                                 #use theme clasic
  labs( y = "Trophic Position", x = "Longitudinal Gradient (PC1)")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 9, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)


ggsave(filename = "figures/pdfs/WHS_tpvsnatfocustax_pc1.pdf", plot = last_plot(), 
       device = "pdf", width = 8, height = 8, units = "in", limitsize = FALSE)
###just ffg

focusffg <- corr_reg_pc1_ffg$ffg

#pc1

TPplot_pc1 <- fishdataall %>% 
  filter(model %in% focusffg) %>% 
  ggplot(aes(x = PC1, y = TP_nobase))+                       #plot TP vs LUlc
  geom_point(aes(x = PC1, y = TP_nobase), color = "grey")+                                                    #plot points
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+                                        #Plot linear model
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  facet_rep_grid(model~taxon_code)+                                         #separately based on taxon_code
  theme_classic()+                                                 #use theme clasic
  labs( y = "Trophic Position", x = "Longitudinal Gradient (PC1)")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 10, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8))+
  scale_x_continuous(limits = c(0,10), breaks = c(0,3,6,9))+
  stat_cor(label.y = 9, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)


ggsave(filename = "figures/tpvsnatfocusffg_pc1.png", plot = last_plot(), 
       device = "png", width = 8, height = 8, units = "in", limitsize = FALSE)

ggsave(filename = "figures/pdfs/tpvsnatfocusffg_pc1.pdf", plot = last_plot(), 
       device = "pdf", width = 8, height = 8, units = "in", limitsize = FALSE)

#####MODEL INFORMATON USING NESTED DATA#####

TP_fnx_pc1 <- function(df) {lm(TP~PC1, data = df)}


#using pc1

TP_nest_pc1 <- fishdataall %>% 
  dplyr::select(TP,PC1,model,taxon_code) %>% 
  group_by(taxon_code,model) %>% 
  nest() %>% 
  mutate(TP_lm = map(data,TP_fnx_pc1),
         TP_tidy = map(TP_lm,tidy),
         TP_glance = map(TP_lm, glance),
         conf_int = map(TP_lm, ~ confint(.x, parm = 2)),
         conf_tidy = map(conf_int,tidy))

TP_tidy_b0_pc1 <- TP_nest_pc1 %>% 
  unnest(TP_tidy, .drop = TRUE) %>% 
  filter(term == "(Intercept)") %>% 
  mutate(Intercept = estimate) %>% 
  dplyr::select(taxon_code,model,Intercept)

TP_tidy_b1_pc1 <- TP_nest_pc1 %>% 
  unnest(TP_tidy, .drop = TRUE) %>% 
  filter(term == "PC1") %>% 
  mutate(slope = estimate) %>% 
  dplyr::select(taxon_code, model, slope)

TP_glance_info_pc1 <- TP_nest_pc1 %>% 
  unnest(TP_glance, .drop = TRUE) %>% 
dplyr::select(taxon_code,model,p.value, adj.r.squared) 
######Supplemental Table 2######

TP_confint<- TP_nest_pc1 %>% 
  unnest(conf_tidy) 

TP_confint <- TP_confint %>% 
  dplyr::select(model,taxon_code, X2.5.., X97.5..) %>% 
  rename(b1_low = X2.5..,
         b1_high = X97.5..)

####using pc1

TP_table_dat_pc1 <- TP_glance_info_pc1 %>% 
  left_join(TP_tidy_b0_pc1, by = c("model", "taxon_code"))%>% 
  left_join(TP_tidy_b1_pc1, by = c("model","taxon_code")) %>% 
  left_join(TP_confint, by = c("model", "taxon_code")) %>% 
  dplyr::select(taxon_code, model, adj.r.squared, Intercept, slope,b1_low,b1_high, p.value)%>%
  mutate_if(is.numeric, round, 3) %>% 
  unite(conf_int, b1_low, b1_high, sep = ",") %>% 
  unite(b1_confint, slope, conf_int, sep = "(") %>% 
  unite(modelinfo,adj.r.squared,Intercept,b1_confint,p.value, sep = "/") %>% 
  spread(key = taxon_code, value = modelinfo) %>% 
  separate(BNT, c("BNT_adj.r.squared", "BNT_Intercept", "BNT_slope", "BNT_p.value"), sep = "/") %>% 
  separate(CKC, c("CKC_adj.r.squared", "CKC_Intercept", "CKC_slope", "CKC_p.value"), sep = "/") %>% 
  separate(WHS, c("WHS_adj.r.squared", "WHS_Intercept", "WHS_slope", "WHS_p.value"), sep = "/") %>% 
  mutate_at(as.numeric, .vars = c(5,9,13)) %>% 
  mutate(BNT_p.value = if_else(condition = BNT_p.value==0, true = paste("<0.001"), false = paste(BNT_p.value))) %>% 
  mutate(CKC_p.value = if_else(condition = CKC_p.value==0, true = paste("<0.001"), false = paste(CKC_p.value))) %>% 
  mutate(WHS_p.value = if_else(condition = WHS_p.value==0, true = paste("<0.001"), false = paste(WHS_p.value)))

write_csv(TP_table_dat_pc1, "TP_lm_table_pc1.csv")

slopechange <- TP_tidy_b1_pc1 %>% 
  left_join(uncor_b1, by = "taxon_code") %>% 
  mutate(b1chan = uncor_b1-slope)

range(slopechange$b1chan)
mean(slopechange$b1chan)

hist(slopechange$b1chan)

sc_ft <- slopechange %>% 
  filter(model %in% focustaxon$taxon_code) %>% 
  group_by(taxon_code) %>% 
  summarise(mean_change = mean(b1chan, na.rm = T),
            low_change = min(b1chan),
            high_change = max(b1chan))

sc_ffg <- slopechange %>% 
  filter(model %in% focusffg) %>% 
  filter(model != "Shredder") %>% 
  group_by(taxon_code) %>% 
  summarise(mean_change = mean(b1chan, na.rm = T),
            low_change = min(b1chan),
            high_change = max(b1chan))

sc_ephsimfilt <- slopechange %>% 
  filter(model == "Simulidae" | model == "Ephemeridae" | model == "Filterer")




#----------------------------------------------------------------------------#
#############ANCOVAS for all plots##########
#----------------------------------------------------------------------------#

uncorrecteddata <- fishdataall %>% 
  filter(model == "Predator") %>% 
  mutate(model = "Uncorrected") %>% 
  dplyr::select(- TP) %>% 
  rename(TP = TP_nobase)

ancovadata <- fishdataall %>% 
  bind_rows(uncorrecteddata) %>% 
  dplyr::select(-TP_nobase)

library(car)
library(tidyverse)
#####BNT by baseline model####

#Predator vs uncorrected 
BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "BNT") %>% 
  filter(model == "Uncorrected"| model == "Predator" | model == "Collector"
         |model == "Filterer"| model == "Grazer")

BNT_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

save_plot(filename = "figures/pdfs/current.pdf", plot = last_plot(), device = pdf(),
          units = "in", base_width = 8.5, base_height = 8.5)

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "BNT") %>% 
  filter(model == "Predator" | model == "Collector"
         |model == "Filterer"| model == "Grazer")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

#BNT~uncorrected vs predator
BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "BNT") %>% 
  filter(model == "Uncorrected"| model == "Predator" | model == "Collector"
         |model == "Filterer"| model == "Grazer")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "BNT") %>% 
  filter(model == "Uncorrected"| model == "Predator" | model == "Collector"
         |model == "Filterer"| model == "Grazer")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")
#####CKC by baseline model####

#all vs uncorrected 
BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "CKC") %>% 
  filter(model == "Uncorrected"| model == "Predator" | model == "Collector"
         |model == "Filterer"| model == "Grazer")

BNT_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

save_plot(filename = "figures/pdfs/current.pdf", plot = last_plot(), device = pdf(),
          units = "in", base_width = 8.5, base_height = 8.5)

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "CKC") %>% 
  filter(model == "Uncorrected"| model == "Collector"
         |model == "Filterer"| model == "Grazer")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

#####WHS by baseline model####

#Predator vs uncorrected 
BNT_ancovadata <- ancovadata %>% 
  filter(taxon_code == "WHS") %>% 
  filter(model == "Uncorrected"| model == "Predator" | model == "Collector"
         |model == "Filterer"| model == "Grazer")

BNT_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

save_plot(filename = "figures/pdfs/current.pdf", plot = last_plot(), device = pdf(),
          units = "in", base_width = 8.5, base_height = 8.5)

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")


#########BNT and Taxa###########

taxa_ancovadata <- ancovadata %>% 
  filter(taxon_code == "BNT") %>% 
  filter(model == "Uncorrected" | model == "Baetidae" | model == "Chironomidae" |
           model == "Dytiscidae" | model == "Elmidae-adult" | model == "Elmidae-larvae" |
           model == "Heptaganeidae" | model == "Hydropyschidae" | model == "Simulidae")

taxa_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

taxa_ancovadata <- ancovadata %>% 
  filter(taxon_code == "BNT") %>% 
  filter( model == "Baetidae" | model == "Chironomidae" |
           model == "Dytiscidae" | model == "Elmidae-adult" | model == "Elmidae-larvae" |
           model == "Heptaganeidae" | model == "Hydropyschidae" | model == "Simulidae")

taxa_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

####wHS and Taxa####

taxa_ancovadata <- ancovadata %>% 
  filter(taxon_code == "WHS") %>% 
  filter(model == "Uncorrected" | model == "Baetidae" | model == "Chironomidae" |
           model == "Dytiscidae" | model == "Elmidae-adult" | model == "Elmidae-larvae" |
           model == "Heptaganeidae" | model == "Hydropyschidae" | model == "Simulidae")

taxa_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

taxa_ancovadata <- ancovadata %>% 
  filter(taxon_code == "WHS") %>% 
  filter( model == "Baetidae" | model == "Chironomidae" |
            model == "Dytiscidae" | model == "Elmidae-adult" | model == "Elmidae-larvae" |
            model == "Heptaganeidae" | model == "Hydropyschidae" | model == "Simulidae")

taxa_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

####CKC and Taxa####

taxa_ancovadata <- ancovadata %>% 
  filter(taxon_code == "CKC") %>% 
  filter(model == "Uncorrected" | model == "Baetidae" | model == "Chironomidae" |
           model == "Dytiscidae" | model == "Elmidae-adult" | model == "Elmidae-larvae" |
           model == "Heptaganeidae" | model == "Hydropyschidae" | model == "Simulidae")

taxa_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")

taxa_ancovadata <- ancovadata %>% 
  filter(taxon_code == "CKC") %>% 
  filter( model == "Baetidae" | model == "Chironomidae" |
            model == "Dytiscidae" | model == "Elmidae-adult" | model == "Elmidae-larvae" |
            model == "Heptaganeidae" | model == "Hydropyschidae" | model == "Simulidae")

taxa_ancovadata %>% 
  ggplot(aes(x = PC1, y = TP, shape = model))+
  geom_point(aes(color = model))+
  geom_smooth(aes(color = model), method = "lm")

hoRS <- aov(TP ~ PC1 * model, data = BNT_ancovadata)
Anova(hoRS, type="III")


read.csv()