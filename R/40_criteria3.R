
# Criteria #3: models of baseline, basal resources, and ffg
# do the different potential baselines responded similarly to PC1 changes?

# Plot model outputs for all the groups


## Visualize --------------------

# by fishes
alldat |> 
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |>
  # filter(compartment=="fish") |> 
  ggplot(aes(PC1, d15N, color = taxon_code)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

# resources
alldat |> 
  filter(taxon_code %in% c("Biofilm","FBOM", "filimentous", "Seston")) |> 
  ggplot(aes(PC1, d15N, color = taxon_code)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

# inverts by ffg
invertdata_sub |> 
  ggplot(aes(PC1, d15N, color = ffg)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 

# inverts by taxa
invertdata_sub |> 
  ggplot(aes(PC1, d15N, color = taxon_code)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) 


## Fit models --------------------------

# ffgs
mods_ffg <- invertdata |> 
  group_by(ffg) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

# taxa
mods_taxa <- invertdata |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

mods_fish <- alldat |> 
  filter(taxon_code %in% c("BNT", "CKC", "WHS", "LND", "LNS")) |>
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

mods_baseline <- alldat |> 
  filter(taxon_code %in% c("Biofilm","FBOM", "filimentous", "Seston")) |> 
  group_by(taxon_code) |> 
  nest() |> 
  mutate(
    fit_pc1 = map(data, ~ lm(d15N ~ PC1, data = .x)),
    tidy = map(fit_pc1, tidy),
    glance = map(fit_pc1, glance),
    augment = map(fit_pc1, augment),
    conf_int = map(fit_pc1, ~ confint(.x, parm = 2))
  ) 

mods_ffg |> unnest(glance)
mods_taxa |> unnest(glance) |> print(n=100)
mods_fish |> unnest(glance)


temp <- mods_baseline |> 
  unnest(tidy)|>
  unnest(conf_int)|>
  filter(term == "PC1")|>
  select(11)



temp2 <- as.numeric(temp$conf_int[,1])

temp[[11]]

temp$conf_int[1:4]
## Plot ---------

#plot f(x)n

plot_sigcorr <- function(df, x_cat){
  upp <- df |>
    unnest(tidy)|>
    unnest(conf_int)|>
    filter(term == "PC1")|
    select(Conf)
           
           sig = if_else(low > 0, true = "yes", false = "no"))|>
    ggplot(aes(x = fct_reorder({{ x_cat }}, estimate, median), y = estimate, fill = sig)) +
    geom_linerange(ymin = low, ymax = upp) +
    geom_point(size = 2, color = "black", shape = 21) + 
    coord_flip()+
    scale_y_continuous(limits = c(0.00,1.02), breaks = seq(0.00,1.00,.25)) +
    labs(y = "", fill = "P < 0.05") + 
    theme_bw(base_size = 12) + 
    theme(
      axis.line = element_line(linewidth = .5),
      axis.ticks.length = unit(.25, "cm"),
      axis.title.y = element_text(vjust = 2), 
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

p5 <- mods_ffg|>plot_sigcorr(x_cat = ffg)+labs(x = "Feeding Group")
p6 <- mods_taxa|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Taxonomic Group")
p7 <- mods_fish|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Fish Species")
p8 <- mods_baseline|>plot_sigcorr(x_cat = taxon_code)+labs(x = "Basal Resource")


plot_sigcorr <- function(df, x siglevel = 0.05) {
  df |>
    unnest(glance)|>
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
}



#### OLD CODE #######



## Data ----------------------------------

pc1dat <- invertdata %>% 
  group_by(taxon_code, PC1) %>% 
  summarise(mean_d15n = mean(d15N))

pc1dat_ffg <- invertdata %>% 
  group_by(ffg, PC1) %>% 
  summarise(mean_d15n = mean(d15N))

pc1dat_pp <- primprodat %>% 
  group_by(taxon_code, PC1) %>% 
  summarise(mean_d15N = mean(mean_d15N))


## pc1 --------------

pc1_corr <- pc1dat %>% 
  select(taxon_code, PC1, mean_d15n) %>% 
  group_by(taxon_code) %>% 
  nest()

pc1_corr_ffg <- pc1dat_ffg %>% 
  select(ffg, PC1, mean_d15n) %>% 
  group_by(ffg) %>% 
  nest()

pc1_corr_pp <- pc1dat_pp %>% 
  rename(mean_d15n = mean_d15N) %>% 
  select(taxon_code, PC1, mean_d15n) %>% 
  group_by(taxon_code) %>% 
  nest()


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
  lemon::facet_rep_wrap(~taxon_code)+                                                        #separate by taxon
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
  ggpubr::stat_cor(label.y = 12, aes(label = paste(..rr.label.., format_pval(..p..), sep = "*`,`~")), size = 3)


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