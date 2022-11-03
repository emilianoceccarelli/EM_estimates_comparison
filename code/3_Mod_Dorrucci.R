
# Function and knots position research ------------------------------------

memory.limit(size = 99999999999)
all_comb_knots = combn(53, 4)

registerDoParallel(cores = detectCores() -1)

data_under_review = mortality_data %>%
  filter(year>2014) %>%
  group_by(week, year) %>%
  summarise(deaths = sum(deaths)) %>%
  mutate(year = ifelse(year %in% 2015:2019,
                       "baseline",
                       year)) %>%
  group_by(year, week) %>%
  summarise(deaths = mean(deaths))%>%
  mutate(year = factor(year),
         deaths = round(deaths, 0)) %>%
  mutate(year = relevel(year,  ref = "baseline"))

aic_knots = vector()

for(i in 1:dim(all_comb_knots)[2]){
  aic_knots[i] = glm.nb(data = data_under_review,
                        formula = deaths ~ year + year*bs(data_under_review$week,
                                                          degree = 2,
                                                          knots = c(all_comb_knots[,i])),
                        link = log)$aic
}


best_knots = all_comb_knots[,which.min(aic_knots)]



excess_Dorrucci_fn <- function(dataset, baseline, cv = F){
  
  if(cv){
    RMSE = vector()
    for(i in baseline){
      
      data_under_review = dataset %>%
        filter(year %in% baseline) %>%
        mutate(year = ifelse(year == i,
                             i,
                             "baseline")) %>%
        group_by(year, week) %>%
        summarise(deaths = mean(deaths))%>%
        mutate(year = factor(year),
               deaths = round(deaths, 0)) %>%
        mutate(year = relevel(year,  ref = "baseline"))
      
        model <- glm.nb(data = data_under_review,
                      formula = deaths ~ year + year*bs(week,
                                                        degree = 2,
                                                        knots = best_knots),
                      link = log)
      prediction = add_ci(data_under_review,
                          model,
                          names = c("LB", "UB")) %>%
        mutate(pred = round(pred, 0),
               LB = round(LB, 0),
               UB = round(UB, 0))
      
      
      RMSE[i-min(baseline)+1] = sqrt(mean((prediction$pred[prediction$year == i]- prediction$pred[prediction$year == "baseline"])^2))
      
      
    }
    return(RMSE)
  }else{
    
    data_under_review = dataset %>%
      filter(year >= min(baseline)) %>%
      mutate(year = ifelse(year %in% baseline,
                           "baseline",
                           year)) %>%
      group_by(year, week) %>%
      summarise(deaths = mean(deaths))%>%
      mutate(year = factor(year),
             deaths = round(deaths, 0)) %>%
      mutate(year = relevel(year,  ref = "baseline"))
    
    model <- glm.nb(data = data_under_review,
                    formula = deaths ~ year + year*bs(week,
                                                      degree = 2,
                                                      knots = best_knots),
                    link = log)
    prediction = add_ci(data_under_review, alpha = 0.05,
                        model,
                        names = c("LB", "UB"),
                        type = "parametric") %>%
      mutate(pred = round(pred, 0),
             LB = round(LB, 0),
             UB = round(UB, 0))
    
    excess_2020 = sum(prediction$pred[prediction$year == 2020]- prediction$pred[prediction$year == "baseline"])
    excess_2020_UB = sum(prediction$UB[prediction$year == 2020]- prediction$UB[prediction$year == "baseline"])
    excess_2020_LB = sum(prediction$LB[prediction$year == 2020]- prediction$LB[prediction$year == "baseline"])
    
    excess_2021 = sum(prediction$pred[prediction$year == 2021]- prediction$pred[prediction$year == "baseline"])
    excess_2021_UB = sum(prediction$UB[prediction$year == 2021]- prediction$UB[prediction$year == "baseline"])
    excess_2021_LB = sum(prediction$LB[prediction$year == 2021]- prediction$LB[prediction$year == "baseline"])
    
    excess_2020_ci = paste0(excess_2020,
                            " (",
                            min(excess_2020_UB, excess_2020_LB),
                            ", ",
                            max(excess_2020_UB, excess_2020_LB),
                            ")")
    excess_2021_ci = paste0(excess_2021,
                            " (",
                            min(excess_2021_UB, excess_2021_LB),
                            ", ",
                            max(excess_2021_UB, excess_2021_LB),
                            ")")
    mean_15_19 = sum(prediction$pred[prediction$year == "baseline"])
    
    data_epi_plot = prediction %>%
      filter(week != 53)
    data_epi_plot = data.frame(week = c(seq(ymd("2020-01-01"),
                                            ymd("2020-12-31"),
                                            by = '1 week')[1:52],
                                        seq(ymd("2021-01-01"),
                                            ymd("2021-12-31"),
                                            by = '1 week')[2:53]),
                               baseline = rep(prediction%>%
                                                filter(year == "baseline", week != 53) %>%
                                                .$pred,
                                              2),
                               LB = rep(prediction%>%
                                          filter(year == "baseline", week != 53) %>%
                                          .$LB,
                                        2),
                               UB = rep(prediction%>%
                                          filter(year == "baseline", week != 53) %>%
                                          .$UB,
                                        2),
                               real = c(prediction%>%
                                          filter(year == "2020", week != 53) %>%
                                          .$death,
                                        prediction%>%
                                          filter(year == "2021", week != 53) %>%
                                          .$death))
    epi_plot = ggplot(data_epi_plot,
                      aes(x = week))+
      geom_path(aes(y = real, color = "Number of deaths"), lwd=.8)+
      geom_line(aes(y = baseline, color = "Expected number of deaths"), lwd=.8)+
      geom_ribbon(aes(ymin=LB,ymax=UB),alpha=0.3)+
      theme_epicurve()+
      scale_color_manual(name =" ",
                         values = c(`Number of deaths` = "Red",
                                    `Expected number of deaths` = "Blue"))+
      ylab("Weekly number of deaths")+
      scale_x_date(date_breaks = "1 months", date_labels = "%b %Y")
    
    
    return(list(fit_2020 = excess_2020,
                fit_2021 = excess_2021,
                UB_2020 = excess_2020_UB,
                LB_2020 = excess_2020_LB,
                UB_2021 = excess_2021_UB,
                LB_2021 = excess_2021_LB,
                `Excess 2020` = excess_2020_ci,
                `Excess 2021` = excess_2021_ci,
                mean_15_19 = mean_15_19,
                epi_plot = epi_plot))
  }
  
}



# Estimate Italy and cross-validation analysis ----------------------------


mod_Dorrucci_Italy <- excess_Dorrucci_fn(mortality_data %>%
                                           group_by(week, year) %>%
                                           summarise(deaths = sum(deaths)),
                                         baseline = 2015:2019,
                                         cv = F)
mod_Dorrucci_Italy$epi_plot

mean(excess_Dorrucci_fn(mortality_data %>%
                          group_by(week, year) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2011:2019,
                        cv = T))


# Excess area and age -----------------------------------------------------

results_area_age_Dorrucci <- matrix(NA, 20, 4)
rownames(results_area_age_Dorrucci) <- paste0(rep(area, each = 5),
                                              " ",
                                              rep(age, 4))
colnames(results_area_age_Dorrucci) <- c("2020", "2021", "% 2020", "% 2021")
list_models_area_age_Dorrucci = list()
z = 1
for(i in area){
  for(j in age){
    
    list_models_area_age_Dorrucci[[z]] = if(i == "Italy" & j == "All ages"){
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(year>2014) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2015:2019,
                         cv = F)
      
    }else if(i == "Italy" & j != "All ages"){
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(age_class == j, year>2014) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2015:2019,
                         cv = F)
    }else if(i != "Italy" & j == "All ages"){
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(area == i, year>2014) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2015:2019,
                         cv = F)
    }else{
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(area == i, age_class == j, year>2014) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2015:2019,
                         cv = F)
      
    }
    z = z + 1
  }
}


for (i in 1:20) {
  results_area_age_Dorrucci[i,1] = list_models_area_age_Dorrucci[[i]]$`Excess 2020`
  results_area_age_Dorrucci[i,2] = list_models_area_age_Dorrucci[[i]]$`Excess 2021`
  results_area_age_Dorrucci[i,3] = round(list_models_area_age_Dorrucci[[i]]$fit_2020/list_models_area_age_Dorrucci[[i]]$mean_15_19, 4)*100
  results_area_age_Dorrucci[i,4] = round(list_models_area_age_Dorrucci[[i]]$fit_2021/list_models_area_age_Dorrucci[[i]]$mean_15_19, 4)*100
}
results_area_age_Dorrucci

# Excess region -----------------------------------------------------------

results_region_Dorrucci <- matrix(NA, 21, 4)
rownames(results_region_Dorrucci) <- unique(mortality_data$region)
colnames(results_region_Dorrucci) <- c("2020", "2021", "% 2020", "% 2021")

list_models_region_Dorrucci = list()

z = 1
for(i in unique(mortality_data$region)){
  
  list_models_region_Dorrucci[[z]] =  excess_Dorrucci_fn(mortality_data %>%
                                                           filter(region == i, year>2014) %>%
                                                           group_by(week, year) %>%
                                                           summarise(deaths = sum(deaths)),
                                                         baseline = 2015:2019,
                                                         cv = F)
  z = z + 1
  
}


for (i in 1:21) {
  results_region_Dorrucci[i,1] = list_models_region_Dorrucci[[i]]$`Excess 2020`
  results_region_Dorrucci[i,2] = list_models_region_Dorrucci[[i]]$`Excess 2021`
  results_region_Dorrucci[i,3] = round(list_models_region_Dorrucci[[i]]$fit_2020/list_models_region_Dorrucci[[i]]$mean_15_19, 4)*100
  results_region_Dorrucci[i,4] = round(list_models_region_Dorrucci[[i]]$fit_2021/list_models_region_Dorrucci[[i]]$mean_15_19, 4)*100
}
results_region_Dorrucci

# Excess area and age change of baseline ----------------------------------

results_area_age_Dorrucci_1119 <- matrix(NA, 20, 4)
rownames(results_area_age_Dorrucci_1119) <- paste0(rep(area, each = 5),
                                              " ",
                                              rep(age, 4))
colnames(results_area_age_Dorrucci_1119) <- c("2020", "2021", "% 2020", "% 2021")
list_models_area_age_Dorrucci_1119 = list()
z = 1
for(i in area){
  for(j in age){
    
    list_models_area_age_Dorrucci_1119[[z]] = if(i == "Italy" & j == "All ages"){
      
      excess_Dorrucci_fn(mortality_data  %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2011:2019,
                         cv = F)
      
    }else if(i == "Italy" & j != "All ages"){
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(age_class == j) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2011:2019,
                         cv = F)
    }else if(i != "Italy" & j == "All ages"){
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(area == i) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2011:2019,
                         cv = F)
    }else{
      
      excess_Dorrucci_fn(mortality_data %>%
                           filter(area == i, age_class == j) %>%
                           group_by(week, year) %>%
                           summarise(deaths = sum(deaths)),
                         baseline = 2011:2019,
                         cv = F)
      
    }
    z = z + 1
  }
}


for (i in 1:20) {
  results_area_age_Dorrucci_1119[i,1] = list_models_area_age_Dorrucci_1119[[i]]$`Excess 2020`
  results_area_age_Dorrucci_1119[i,2] = list_models_area_age_Dorrucci_1119[[i]]$`Excess 2021`
  results_area_age_Dorrucci_1119[i,3] = round(list_models_area_age_Dorrucci_1119[[i]]$fit_2020/list_models_area_age_Dorrucci_1119[[i]]$mean_15_19, 4)*100
  results_area_age_Dorrucci_1119[i,4] = round(list_models_area_age_Dorrucci_1119[[i]]$fit_2021/list_models_area_age_Dorrucci_1119[[i]]$mean_15_19, 4)*100
}
results_area_age_Dorrucci_1119
