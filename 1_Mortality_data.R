library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(merTools)
library(stringr)
library(lme4)
library(lmerTest)
require(merTools)
library(doParallel)
library(foreach)
library(splines)
library(MASS)
library(gridExtra)
library("ciTools")
library(dlnm)
library(pbs)
library(mixmeta)
library(tidyverse)

library(curl)
library(tsModel) 
library(scales)
library(sp)
library(sf)
library(grid)
library(MASS)
library(tidyr)
library(forcats)

setwd("C:/Users/ceccarelli_emiliano/OneDrive - Istituto Superiore di Sanità/ISS/Stima eccesso")
memory.limit(size = 99999999999)


mortality_data <- read.table("data/c_index.csv", sep = ",", header = T, fileEncoding = "Latin1") %>%
  dplyr::select("PRO_COM_T", "COD_REG20", "COD_REG21", "NOME_NUTS2")%>%
  mutate(COD_REG21 = COD_REG20/10) %>%
  left_join(., read.csv("data/comuni_giornaliero_31dicembre21.csv", sep = ","),
            by = c("PRO_COM_T" = "COD_PROVCOM")) %>%
  dplyr::select(NOME_NUTS2, REG, GE, CL_ETA, T_11, T_12, T_13, T_14, T_15, T_16, T_17, T_18, T_19, T_20, T_21) %>%
  mutate(across(everything(), ~replace(., . ==  "n.d." , 0)),
         across(T_11:T_21, ~as.numeric(.))) %>%
  mutate(CL_ETA = ifelse(CL_ETA %in% 0:10,
                         "0-49",
                         ifelse(CL_ETA %in% 11:13,
                                "50-64",
                                ifelse(CL_ETA %in% 14:16,
                                       "65-79",
                                       "80+"))),
         macroarea = ifelse(REG %in% 1:8,
                            "North Italy",
                            ifelse(REG %in% 9:12,
                                   "Center Italy",
                                   "South Italy and islands")))%>%
  group_by(NOME_NUTS2, macroarea, CL_ETA, GE) %>%
  summarise(across(T_11:T_21, list(sum)))%>%
  `colnames<-`(c("Regione", "Macroarea", "CL_ETA", "Data", paste0("a", 2011:2021))) %>%
  mutate(Data = ifelse(Data==229, 228, Data))%>%
  mutate(Data = as.character(Data))%>%
  mutate(Data = ifelse(str_count(Data)==3,
                       paste0("0", Data),
                       Data)) %>%
  group_by(Regione, Macroarea, CL_ETA, Data)%>%
  summarise(across(a2011:a2021, list(sum)))

mortality_data <- bind_cols(Data = as.Date(paste0(rep(c(2011:2021),
                                                           nrow(mortality_data)),
                                                       rep(mortality_data$Data,
                                                           each = 11)),
                                                "%Y%m%d"),
                                 region = rep(mortality_data$Regione,
                                               each = 11),
                                 area = rep(mortality_data$Macroarea,
                                                 each = 11),
                                 age_class = rep(mortality_data$CL_ETA,
                                              each = 11),
                                 deaths = c(t(as.matrix(mortality_data[,5:15]))))%>%
  arrange(Data) %>%
  mutate(week = week(Data),
         month = month(Data),
         year = year(Data)) %>%
  group_by(week, month, year, region, area, age_class) %>%
  summarise(deaths = sum(deaths)) %>%
  mutate(region = ifelse(region == "Provincia Autonoma di Bolzano/Bozen",
                          "P.A. Bolzano",
                          ifelse(region == "Provincia Autonoma di Trento",
                                 "P.A. Trento",
                                 ifelse(region == "Valle d'Aosta/Vallée d'Aoste",
                                        "Valle d'Aosta",region))))





theme_epicurve <- function(...) {
  theme(axis.text.x = element_text(size=11, angle = 90, color = "black"),
        axis.text.y = element_text(size=11, color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom")
}




 
  

library(lubridate)
ese = data.frame(data = seq(dmy("01/01/2021"),dmy("31/12/2021"),by="day"))
ese$week = week(ese$data)
ese$isoweek = isoweek(ese$data)
head(ese)
plot(ese$week, ese$isoweek)

