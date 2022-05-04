
# Criteria 4: Trophic Estimates independent of Land Use

## Prep

# libraries
library(tidyverse) 
library(here)


# Data ----------------

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
