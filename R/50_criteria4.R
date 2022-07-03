
# Criteria 4: Trophic Estimates independent of Land Use

## Prep

# libraries
library(tidyverse) 
library(gridExtra)


# Data ----------------
#find the mean d15N for each baseline proxy at our sites across our three sampling periods


meanNdatall <- invertdata %>% 
  group_by(taxon_code,site_id) %>%                        #group by site and taxon
  summarise(                                              #summarise get sample size, mean, standard deviation and CV
    mean_d15n = mean(d15N, na.rm = TRUE)) %>% 
  spread(taxon_code,mean_d15n)

ffg_meanN <- invertdata %>% 
  group_by(ffg,site_id) %>% 
  summarise(mean_d15n = mean(d15N, na.rm = TRUE)) %>% 
  spread(ffg, mean_d15n)

#####Trophic Position Estimate #####

#get the isotope data for our three target fish species

fishdat <- isodat %>% 
  filter(resource == "fish") %>% 
  filter(taxon_code == "CKC"|
           taxon_code == "WHS"|
           taxon_code == "BNT") 

#estimate the TP of each of our target fish sample without correction from baseline

uncorrected_data<- fishdat %>% 
  left_join(PCAdat, by = "site_id") %>% 
  mutate(discFactor = ((-0.281*d15N)+5.879)) %>% 
  dplyr::select(-sample_year, -d13C, -stream_name, -unique_id, -sia_sample_id, -sample_hitch,
                -compartment, -resource, -length_mm, -d13C_scl, -d15N_scl, -CN, -lulc_del, 
                -lulc_ag, -PC2) %>% 
  mutate(TP_nobase = d15N/discFactor) 


#quick plot to view the uncorrected TP data

uncorrected_data %>% 
  ggplot(aes (x = PC1, y = TP_nobase, color = taxon_code))+
  geom_point()+
  geom_smooth(method = "lm")

#get regression statistics for uncorrected TP vs PC1 for each fish

fit.bnt <- lm(TP_nobase~PC1, data = uncorrected_data[uncorrected_data$taxon_code == "BNT",])
summary(fit.bnt)

     #int = 0.860, slope = 0.506, R2 = 0.227, F(1,371) = 110.1, p < 0.001

fit.ckc <- lm(TP_nobase~PC1, data = uncorrected_data[uncorrected_data$taxon_code == "CKC",])
summary(fit.ckc)

     #int = 2.679, slope = 0.141, R2 = 0.024, F(1,358) = 9.798, p = 0.002

fit.whs <- lm(TP_nobase~PC1, data = uncorrected_data[uncorrected_data$taxon_code == "WHS",])
summary(fit.whs)

     #int = 1.4739, slope = 0.201, R2 = 0.103, F(1,289) = 34.43, p < 0.001

#estimate the TP while correcting values with our baseline proxies

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

######Plot Trophic Position vs. gradient ######

#####  taxa proxies.............

#brown trout plots

bnt_uncorrect <- uncorrected_data %>% 
  filter(taxon_code == "BNT")

#bnt corrected by baetidae

fit.bnt <- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "BNT"& corrected_data$model == "Baetidae",])
summary(fit.bnt)

bnt.baetidae <- corrected_data %>% 
  filter(model == "Baetidae") %>%
  filter(taxon_code == "BNT") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = bnt_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = " ", x = NULL, title = "A) BNT by Baetidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(1,9), breaks = c(1,3,5,7,9))+
  scale_x_continuous(limits = c(1,9), breaks = c(0,2,4,6,8))+
  annotate(geom = "text", x = 3.2, y = 9, label = expression(UNC:~R^2~"="~0.227:~italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 3.2, y = 8.5, label = expression(COR:~R^2~"="~0.024:~italic(P)~"="~0.004), size = 2.5)
  
#bnt corrected by ephemeridae

fit.bnt <- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "BNT"& corrected_data$model == "Ephemeridae",])
summary(fit.bnt)

bnt.ephemeridae<- corrected_data %>% 
  filter(model == "Ephemeridae") %>%
  filter(taxon_code == "BNT") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = bnt_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  theme_classic()+ 
  labs( y = " ", x = NULL, title = "D) BNT by Ephemeridae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(1,9), breaks = c(1,3,5,7,9))+
  scale_x_continuous(limits = c(1,9), breaks = c(0,2,4,6,8))+
  annotate(geom = "text", x = 3.2, y = 9, label = expression(UNC:~R^2~"="~0.227:~italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 1.4, y = 8.4, label = "COR: ns",size = 2.5)

#bnt corrected by Heptaganeidae

fit.bnt <- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "BNT"& corrected_data$model == "Heptaganeidae",])
summary(fit.bnt)

bnt.heptaganeidae <- corrected_data %>% 
  filter(model == "Heptaganeidae") %>%
  filter(taxon_code == "BNT") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = bnt_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  theme_classic()+ 
  labs( y = "Trophic Position", x = NULL, title = "G) BNT by Heptaganeidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(1,9), breaks = c(1,3,5,7,9))+
  scale_x_continuous(limits = c(1,9), breaks = c(0,2,4,6,8))+
  annotate(geom = "text", x = 3.2, y = 9, label = expression(UNC:~R^2~"="~0.227:~italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 1.4, y = 8.4, label = "COR: ns",size = 2.5)

expression(UNC:~R^2~"="~0.227:~italics(P)~"<"~0.001)

#bnt corrected by Hydropyschidae

fit.bnt <- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "BNT"& corrected_data$model == "Hydropyschidae",])
summary(fit.bnt)

bnt.hydropyschidae<- corrected_data %>% 
  filter(model == "Hydropyschidae") %>%
  filter(taxon_code == "BNT") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = bnt_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = " ", x = NULL, title = "J) BNT by Hydropyschidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(1,9), breaks = c(1,3,5,7,9))+
  scale_x_continuous(limits = c(1,9), breaks = c(0,2,4,6,8))+
  annotate(geom = "text", x = 3.2, y = 9, label = expression(UNC:~R^2~"="~0.227:~italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 3.2, y = 8.5, label = expression(COR:~R^2~"="~0.126:~italic(P)~"<"~ 0.001),size = 2.5)

#bnt corrected by Simulidae

fit.bnt <- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "BNT"& corrected_data$model == "Simulidae",])
summary(fit.bnt)

bnt.simulidae<- corrected_data %>% 
  filter(model == "Simulidae") %>%
  filter(taxon_code == "BNT") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = bnt_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = " ", x = " ", title = "M) BNT by Simulidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(1,9), breaks = c(1,3,5,7,9))+
  scale_x_continuous(limits = c(1,9), breaks = c(0,2,4,6,8))+
  annotate(geom = "text", x = 3.2, y = 9, label = expression(UNC:~R^2~"="~0.227:~italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 3.2, y = 8.4, label = expression(COR:~R^2~"="~0.072:~italic(P)~"<"~0.001),size = 2.5)

#creek chub plots

ckc_uncorrect <- uncorrected_data %>% 
  filter(taxon_code == "CKC")

#ckc corrected by baetidae

fit.ckc<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "CKC"& corrected_data$model == "Baetidae",])
summary(fit.ckc)

ckc.baetidae <- corrected_data %>% 
  filter(model == "Baetidae") %>%
  filter(taxon_code == "CKC") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "B) CKC by Baetidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.024:~italic(P)~"="~0.002),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.5, label = expression(COR:~R^2~"="~0.084:~italic(P)~"<"~0.001),size = 2.5)

#ckc corrected by ephemeridae

fit.ckc<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "CKC"& corrected_data$model == "Ephemeridae",])
summary(fit.ckc)

ckc.ephemeridae <- corrected_data %>% 
  filter(model == "Ephemeridae") %>%
  filter(taxon_code == "CKC") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "E) CKC by Ephemeridae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.024:~italic(P)~"="~0.002),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.5, label = expression(COR:~R^2~"="~0.033:~ italic(P)~"="~0.003),size = 2.5)

#ckc corrected by heptaganeidae

fit.ckc<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "CKC"& corrected_data$model == "Heptaganeidae",])
summary(fit.ckc)

ckc.heptaganeidae<- corrected_data %>% 
  filter(model == "Heptaganeidae") %>%
  filter(taxon_code == "CKC") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "H) CKC by Heptaganeidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.024:~italic(P)~"="~0.002),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.5, label = expression(COR:~R^2~"="~0.023:~ italic(P)~"="~0.004),size = 2.5)

#ckc corrected by hydropyschidae

fit.ckc<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "CKC"& corrected_data$model == "Hydropyschidae",])
summary(fit.ckc)

ckc.hydropyschidae<- corrected_data %>% 
  filter(model == "Hydropyschidae") %>%
  filter(taxon_code == "CKC") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "K) CKC by Hydropyschidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.024:~italic(P)~"="~0.002),size = 2.5)+
  annotate(geom = "text", x = 4.28, y = 7.4, label = "COR: ns",size = 2.5)

#ckc corrected by simulidae

fit.ckc<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "CKC"& corrected_data$model == "Simulidae",])
summary(fit.ckc)

ckc.simulidae<- corrected_data %>% 
  filter(model == "Simulidae") %>%
  filter(taxon_code == "CKC") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  theme_classic()+ 
  labs( y = NULL, x = "Longitudinal Gradient (PC1)", title = "N) CKC by Simulidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.024:~italic(P)~"="~0.002),size = 2.5)+
  annotate(geom = "text", x = 4.28, y = 7.3, label = "COR: ns",size = 2.5)

#white sucker plots

ckc_uncorrect <- uncorrected_data %>% 
  filter(taxon_code == "WHS")

#whs corrected by baetidae

fit.whs<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "WHS"& corrected_data$model == "Baetidae",])
summary(fit.whs)

whs.baetidae <- corrected_data %>% 
  filter(model == "Baetidae") %>%
  filter(taxon_code == "WHS") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "C) WHS by Baetidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.103:~ italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 4.28, y = 7.4, label = "COR: ns",size = 2.5)

#whs corrected by ephemeridae

fit.whs<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "WHS"& corrected_data$model == "Ephemeridae",])
summary(fit.whs)

whs.ephemeridae <- corrected_data %>% 
  filter(model == "Ephemeridae") %>%
  filter(taxon_code == "WHS") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "F) WHS by Ephemeridae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.103:~ italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.5, label = expression(COR:~R^2~"="~0.213:~italic(P)~"<"~0.001),size = 2.5)

#whs corrected by heptaganeidae

fit.whs<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "WHS"& corrected_data$model == "Heptaganeidae",])
summary(fit.whs)

whs.heptaganeidae <- corrected_data %>% 
  filter(model == "Heptaganeidae") %>%
  filter(taxon_code == "WHS") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs(y = NULL, x = NULL, title = "I) WHS by Heptaganeidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.103:~ italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.5, label = expression(COR:~R^2~"="~0.082:~italic(P)~"<"~0.001),size = 2.5)

#whs corrected by Hydropyschidae

fit.whs<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "WHS"& corrected_data$model == "Hydropyschidae",])
summary(fit.whs)

whs.hydropyschidae<- corrected_data %>% 
  filter(model == "Hydropyschidae") %>%
  filter(taxon_code == "WHS") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = NULL, x = NULL, title = "L) WHS by Hydropyschidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.103:~ italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.5, label = expression(COR:~R^2~"="~0.042:~italic(P)~"<"~0.001),size = 2.5)

#whs corrected by Simulidae

fit.whs<- lm(TP ~ PC1, data = corrected_data[corrected_data$taxon_code == "WHS"& corrected_data$model == "Simulidae",])
summary(fit.whs)

whs.simulidae<- corrected_data %>% 
  filter(model == "Simulidae") %>%
  filter(taxon_code == "WHS") %>% 
  ggplot(aes(x = PC1, y = TP))+                       
  geom_point(data = ckc_uncorrect, aes(x = PC1, y = TP_nobase), color = "grey",
             shape = 17)+
  geom_smooth(data = uncorrected_data,aes(x = PC1, y = TP_nobase),method = lm, color = "grey", fill = "#CCCCCC")+ 
  geom_point(aes( x = PC1, y = TP))+
  geom_smooth(aes(x = PC1, y = TP),method = lm, color = "black", fill = "#333333")+ 
  theme_classic()+ 
  labs( y = NULL, x = " ", title = "O) WHS by Simulidae")+                                       #Labels 
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 8, color = "black"),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, color = "black"))+
  scale_y_continuous(limits = c(0,8), breaks = c(1,3,5,7))+
  scale_x_continuous(limits = c(4,9.5), breaks = c(5,7,9))+
  annotate(geom = "text", x = 5.4, y = 8, label = expression(UNC:~R^2~"="~0.103:~ italic(P)~"<"~0.001),size = 2.5)+
  annotate(geom = "text", x = 5.4, y = 7.4, label = expression(COR:~R^2~"="~0.003:~italic(P)~"="~0.003),size = 2.5)

#combine all plots together

cor.taxa <- grid.arrange(bnt.baetidae,ckc.baetidae,whs.baetidae,
             bnt.ephemeridae, ckc.ephemeridae, whs.ephemeridae,
             bnt.heptaganeidae, ckc.heptaganeidae, whs.heptaganeidae,
             bnt.hydropyschidae, ckc.hydropyschidae, whs.hydropyschidae,
             bnt.simulidae,ckc.simulidae,whs.simulidae,
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6),
                                   c(7,8,9),
                                   c(10,11,12),
                                   c(13,14,15)))

ggsave("out/correctedtaxa.png", plot = cor.taxa, device = png, width = 8.5, height = 11, units = "in")

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


