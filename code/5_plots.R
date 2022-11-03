# Scortichini model ----

risultati_rip_eta
risultati_rip_eta_11_19

dati_baseline_Scort <- bind_rows(risultati_rip_eta,
                                 risultati_rip_eta_11_19) %>%
  tibble::add_column(baseline = rep(c("15-21", "11-21"), each = 20)) %>%
  pivot_longer(!c(Rip, baseline),
               names_to = "anno",
               values_to = "eccesso") 
dati_baseline_Scort$eccesso_est = NA
dati_baseline_Scort$LB = NA
dati_baseline_Scort$UB = NA
for(i in 1:nrow(dati_baseline_Scort)){
  
  dati_baseline_Scort$eccesso_est[i] = strsplit(as.character(str_replace_all(dati_baseline_Scort$eccesso,
                                                                          "[(\\,)]",
                                                                          "")),
                                             split = " ")[[i]][1]
  dati_baseline_Scort$LB[i] = strsplit(as.character(str_replace_all(dati_baseline_Scort$eccesso,
                                                                 "[(\\,)]",
                                                                 "")),
                                    split = " ")[[i]][2]
  dati_baseline_Scort$UB[i] = strsplit(as.character(str_replace_all(dati_baseline_Scort$eccesso,
                                                                 "[(\\,)]",
                                                                 "")),
                                    split = " ")[[i]][3]
}
  
dati_baseline_Scort = dati_baseline_Scort %>%
  dplyr::select(-eccesso) %>%
  mutate(eccesso_est = as.numeric(eccesso_est),
         LB = as.numeric(LB),
         UB = as.numeric(UB))

  # pivot_longer(!c(Rip, baseline, anno),
  #              names_to = "tipo_stima",
  #              values_to = "eccesso")


comparison_baseline_20_Scortichini <- ggplot(dati_baseline_Scort %>%
                                            filter(anno == 2020),
                                          aes(x = baseline,
                                              color = baseline))+
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 1.1)+
  facet_wrap(vars(Rip), scales = "free_y")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#bf0603"))+
  labs(caption = "Year 2020, model Scortichini",
       y = "Mortality excess")+
  theme_minimal()

comparison_baseline_21_Scortichini <- ggplot(dati_baseline_Scort %>%
                                               filter(anno == 2021),
                                             aes(x = baseline,
                                                 color = baseline))+
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 1.1)+
  facet_wrap(vars(Rip), scales = "free_y")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#bf0603"))+
  labs(caption = "Year 2021, model Scortichini",
       y = "Mortality excess")+
  theme_minimal()


ggsave(plot = comparison_baseline_20_Scortichini,
       filename = "C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/baseline_Scortichini_20.png",
       width = 12.5, height = 7.5,
       units = "in",
       bg = "white")
ggsave(plot = comparison_baseline_21_Scortichini,
       filename = "C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/baseline_Scortichini_21.png",
       width = 12.5, height = 7.5,
       units = "in",
       bg = "white")



# Maruotti model ----

dati_baseline_Maruotti <- bind_rows(as.data.frame(results_area_age_Maruotti),
                                    as.data.frame(results_area_age_Maruotti_1519)) %>%
  `rownames<-`(NULL) %>%
  dplyr::select(-c(`% 2020`, `% 2021`)) %>%
  tibble::add_column(Rip = rep(risultati_rip_eta$Rip, 2),
                     baseline = rep(c("11-19", "15-19"), each = 20)) %>%
  pivot_longer(!c(Rip, baseline),
               names_to = "anno",
               values_to = "eccesso") 
dati_baseline_Maruotti$eccesso_est = NA
dati_baseline_Maruotti$LB = NA
dati_baseline_Maruotti$UB = NA
for(i in 1:nrow(dati_baseline_Maruotti)){
  
  dati_baseline_Maruotti$eccesso_est[i] = strsplit(as.character(str_replace_all(dati_baseline_Maruotti$eccesso,
                                                                             "[(\\;)]",
                                                                             "")),
                                                split = " ")[[i]][1]
  dati_baseline_Maruotti$LB[i] = strsplit(as.character(str_replace_all(dati_baseline_Maruotti$eccesso,
                                                                    "[(\\;)]",
                                                                    "")),
                                       split = " ")[[i]][2]
  dati_baseline_Maruotti$UB[i] = strsplit(as.character(str_replace_all(dati_baseline_Maruotti$eccesso,
                                                                    "[(\\;)]",
                                                                    "")),
                                       split = " ")[[i]][3]
}
dati_baseline_Maruotti = dati_baseline_Maruotti %>%
  dplyr::select(-eccesso) %>%
  mutate(eccesso_est = as.numeric(eccesso_est),
         LB = as.numeric(LB),
         UB = as.numeric(UB))

comparison_baseline_20_Maruotti <- ggplot(dati_baseline_Maruotti %>%
                                               filter(anno == 2020),
                                             aes(x = baseline,
                                                 color = baseline))+
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 1.1)+
  facet_wrap(vars(Rip), scales = "free_y")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#bf0603"))+
  labs(caption = "Year 2020, model Maruotti",
       y = "Mortality excess")+
  theme_minimal()

comparison_baseline_21_Maruotti <- ggplot(dati_baseline_Maruotti %>%
                                               filter(anno == 2021),
                                             aes(x = baseline,
                                                 color = baseline))+
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 1.1)+
  facet_wrap(vars(Rip), scales = "free_y")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#bf0603"))+
  labs(caption = "Year 2021, model Maruotti",
       y = "Mortality excess")+
  theme_minimal()


ggsave(plot = comparison_baseline_20_Maruotti,
       filename = "C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/baseline_Maruotti_20.png",
       width = 12.5, height = 7.5,
       units = "in",
       bg = "white")
ggsave(plot = comparison_baseline_21_Maruotti,
       filename = "C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/baseline_Maruotti_21.png",
       width = 12.5, height = 7.5,
       units = "in",
       bg = "white")




# Dorrucci model ----

dati_baseline_Dorrucci <- bind_rows(as.data.frame(results_area_age_Dorrucci),
                                    as.data.frame(results_area_age_Dorrucci_1119)) %>%
  `rownames<-`(NULL) %>%
  dplyr::select(-c(`% 2020`, `% 2021`)) %>%
  tibble::add_column(Rip = rep(risultati_rip_eta$Rip, 2),
                     baseline = rep(c("11-19", "15-19"), each = 20)) %>%
  pivot_longer(!c(Rip, baseline),
               names_to = "anno",
               values_to = "eccesso") 
dati_baseline_Dorrucci$eccesso_est = NA
dati_baseline_Dorrucci$LB = NA
dati_baseline_Dorrucci$UB = NA
for(i in 1:nrow(dati_baseline_Dorrucci)){
  
  dati_baseline_Dorrucci$eccesso_est[i] = strsplit(as.character(str_replace_all(dati_baseline_Dorrucci$eccesso,
                                                                                "[(\\,)]",
                                                                                "")),
                                                   split = " ")[[i]][1]
  dati_baseline_Dorrucci$LB[i] = strsplit(as.character(str_replace_all(dati_baseline_Dorrucci$eccesso,
                                                                       "[(\\,)]",
                                                                       "")),
                                          split = " ")[[i]][2]
  dati_baseline_Dorrucci$UB[i] = strsplit(as.character(str_replace_all(dati_baseline_Dorrucci$eccesso,
                                                                       "[(\\,)]",
                                                                       "")),
                                          split = " ")[[i]][3]
}
dati_baseline_Dorrucci = dati_baseline_Dorrucci %>%
  dplyr::select(-eccesso) %>%
  mutate(eccesso_est = as.numeric(eccesso_est),
         LB = as.numeric(LB),
         UB = as.numeric(UB))

comparison_baseline_20_Dorrucci <- ggplot(dati_baseline_Dorrucci %>%
                                            filter(anno == 2020),
                                          aes(x = baseline,
                                              color = baseline))+
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 1.1)+
  facet_wrap(vars(Rip), scales = "free_y")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#bf0603"))+
  labs(caption = "Year 2020, model Dorrucci",
       y = "Mortality excess")+
  theme_minimal()

comparison_baseline_21_Dorrucci <- ggplot(dati_baseline_Dorrucci %>%
                                            filter(anno == 2021),
                                          aes(x = baseline,
                                              color = baseline))+
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 1.1)+
  facet_wrap(vars(Rip), scales = "free_y")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#bf0603"))+
  labs(caption = "Year 2021, model Dorrucci",
       y = "Mortality excess")+
  theme_minimal()


ggsave(plot = comparison_baseline_20_Dorrucci,
       filename = "C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/baseline_Dorrucci_20.png",
       width = 12.5, height = 7.5,
       units = "in",
       bg = "white")
ggsave(plot = comparison_baseline_21_Dorrucci,
       filename = "C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/baseline_Dorrucci_21.png",
       width = 12.5, height = 7.5,
       units = "in",
       bg = "white")



# Comparison Italy all 3 models ----

plot_confronto_baseline <- bind_rows(dati_baseline_Maruotti,
                                     dati_baseline_Dorrucci, 
                                     dati_baseline_Scort)%>%
  filter(Rip == "Italy All ages") %>%
  tibble::add_column(model = rep(c("Model 1: \ntime-series model",
                                   "Model 2: \nepidemiological model",
                                   "Model 3: \ntime-series and temperatures model"),
                                 each = 4)) %>%
  mutate(model = ordered(model,
                         levels = c("Model 1: \ntime-series model",
                                    "Model 2: \nepidemiological model",
                                    "Model 3: \ntime-series and temperatures model"))) %>%
  ggplot(data = .,
         aes(x = baseline,
             color = baseline)) +
  geom_point(aes(y = eccesso_est),
             size = 3)+
  geom_errorbar(aes(ymin = LB,
                    ymax = UB),
                size = 0.8,
                width = 0.7)+
  facet_wrap(vars(anno, model), scales = "free_x")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("#219ebc", "#0077b6", "#bf0603", "#9d0208"))+
  labs(y = "Mortality excess")+
  theme_minimal() +
  ggtitle("Comparison excess mortality in 2020 and 2021 across three models and two baselines")
plot_confronto_baseline
ggsave(plot = plot_confronto_baseline,
       filename = paste0("C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/plot_confronto_baseline.png"),
       width = 8, height = 5,
       units = "in",
       bg = "white")



grid.arrange(grobs = list(list_models_area_age[[20]]$`Plot excess 20-21`+
                             labs(caption = "Model 1: time-series model"),
                           list_models_area_age_Dorrucci[[20]]$epi_plot +
                             labs(caption = "Model 2: epidemiological model"),
                           plot_baseline_Italy_Scortichini+
                             labs(caption = "Model 3: time-series and temperatures model")),
              ncol = 1,
              top = "Italy all ages")
ggsave(plot = arrangeGrob(grobs = list(list_models_area_age[[20]]$`Plot excess 20-21`+
                                         labs(caption = "Model 1: time-series model"),
                                       list_models_area_age_Dorrucci[[20]]$epi_plot +
                                         labs(caption = "Model 2: epidemiological model"),
                                       plot_baseline_Italy_Scortichini+
                                         labs(caption = "Model 3: time-series and temperatures model")),
                          ncol = 1,
                          top = "Italy all ages"),
       filename = paste0("C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/Articolo eccesso di mortalità/Grafici/grafico_baseline.png"),
       width = 8, height = 10,
       units = "in",
       bg = "white")



