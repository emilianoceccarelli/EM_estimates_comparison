
# Function ----------------------------------------------------------------


stima_eccesso_fne <- function(dataset, baseline, label_area, label_age, cv = F){
  require(dplyr)
  require(ggplot2)
  require(lme4)
  require(lubridate)
  require(merTools)
  data_under_review = dataset%>%
    filter(., week!=53)
  
  
  
  if(cv){
    RMSE = vector()
    for (j in baseline) {
      new_baseline = baseline[baseline != j]
      data_under_review_baseline = data_under_review %>%
        filter(., year %in% new_baseline)
      
      X1 <- X2 <- matrix(NA,nrow = length(data_under_review_baseline$deaths),ncol=8)
      for (i in 1:8) {
        X1[,i] <- sin(2*i*pi*(1:length(data_under_review_baseline$deaths))/52)
        X2[,i] <- cos(2*i*pi*(1:length(data_under_review_baseline$deaths))/52)
      }
      
      data <- data.frame(data_under_review_baseline$deaths,X1,X2,rep(new_baseline, each=52))
      names(data) <- c("Deaths","sin.1","sin.2","sin.3","sin.4","sin.5","sin.6","sin.7","sin.8",
                       "cos.1","cos.2","cos.3","cos.4","cos.5","cos.6","cos.7","cos.8","year")
      data$time <- 1:length(new_baseline)
      
      set.seed(1603)
      model <- glmer.nb(Deaths~sin.1+sin.2+sin.3+sin.4+sin.5+sin.6+sin.7+sin.8+
                          cos.1+cos.2+cos.3+cos.4+cos.5+cos.6+cos.7+cos.8+(1|year), data=data)
      
      weekly_series_i = data_under_review %>%
        filter(., year == j)
      #
      Y1 <- Y2 <- matrix(NA,nrow = 52,ncol=8)
      for (i in 1:8) {
        Y1[,i] <- sin(2*i*pi*(469:520)/52)
        Y2[,i] <- cos(2*i*pi*(469:520)/52)
      }
      data.Y <- data.frame(Y1,Y2,rep(new_baseline[length(new_baseline)] ,52))
      names(data.Y) <-  c("sin.1","sin.2","sin.3","sin.4","sin.5","sin.6","sin.7","sin.8",
                          "cos.1","cos.2","cos.3","cos.4","cos.5","cos.6","cos.7","cos.8","year")
      require(merTools)
      boot <- predictInterval(merMod = model, newdata = data.Y,
                              level = 0.95, n.sims = 2000,
                              stat = "median", type="probability",
                              include.resid.var = FALSE)
      RMSE[j-baseline[1]+1] = sqrt(mean((weekly_series_i$deaths - boot$fit)^2))
      
      
    }
    
    return(RMSE)
  }else{
    
    data_under_review_baseline = data_under_review %>%
      filter(., year %in% baseline)
    
    X1 <- X2 <- matrix(NA,nrow = length(data_under_review_baseline$deaths),ncol=8)
    for (i in 1:8) {
      X1[,i] <- sin(2*i*pi*(1:length(data_under_review_baseline$deaths))/52)
      X2[,i] <- cos(2*i*pi*(1:length(data_under_review_baseline$deaths))/52)
    }
    
    data <- data.frame(data_under_review_baseline$deaths, X1,X2,rep(baseline, each=52))
    names(data) <- c("Deaths","sin.1","sin.2","sin.3","sin.4","sin.5","sin.6","sin.7","sin.8",
                     "cos.1","cos.2","cos.3","cos.4","cos.5","cos.6","cos.7","cos.8","year")
    data$time <- 1:nrow(data)
    # fit model ----
    
    set.seed(1603)
    
    model <- glmer.nb(Deaths~sin.1+sin.2+sin.3+sin.4+sin.5+sin.6+sin.7+sin.8+
                        cos.1+cos.2+cos.3+cos.4+cos.5+cos.6+cos.7+cos.8+(1|year),
                      data=data)
    
    plot_baseline <- ggplot(data = reshape::melt(data.frame(week = seq(ymd(paste0(baseline[1],
                                                                                  "-01-01")),
                                                                       ymd(paste0(baseline[length(baseline)],
                                                                                  "-12-31")),
                                                                       by = '1 week')[1:length(data$Deaths)],
                                                            deaths = data$Deaths,
                                                            fit = fitted(model)),
                                                 id.var = "week"),
                            aes(x = week,
                                y = value,
                                colour = variable))+
      geom_path(lwd=.8)+
      theme_epicurve()+
      theme(legend.position = "bottom",
            legend.title = element_blank())+
      scale_colour_manual(values = c("Red", "Blue"),
                          labels = c("Real values", "Predicted Values"))+
      labs(title = paste0("Baseline 2011-2019 ", label_area),
           subtitle = label_age)+
      ylab("")+
      scale_x_date(date_breaks = "1 year", date_labels = "%b %Y")
    
    # Prediction Interval 2020----
    
    weekly_series_2020 = data_under_review %>%
      filter(., year == 2020)
    #
    Y1 <- Y2 <- matrix(NA,nrow = 52,ncol=8)
    for (i in 1:8) {
      Y1[,i] <- sin(2*i*pi*(469:520)/52)
      Y2[,i] <- cos(2*i*pi*(469:520)/52)
    }
    data.Y <- data.frame(Y1,Y2,rep(2019,52)) 
    names(data.Y) <-  c("sin.1","sin.2","sin.3","sin.4","sin.5","sin.6","sin.7","sin.8",
                        "cos.1","cos.2","cos.3","cos.4","cos.5","cos.6","cos.7","cos.8","year")
    
    boot20 <- predictInterval(merMod = model, newdata = data.Y,
                              level = 0.95, n.sims = 2000,
                              stat = "median", type="probability",
                              include.resid.var = FALSE)
    fit_2020 <- round(sum(weekly_series_2020$deaths) - sum(boot20$fit), 0)
    mean_15_19 <- sum(boot20$fit)
    UB_2020 <- round(sum(weekly_series_2020$deaths) - sum(boot20$upr), 0)
    LB_2020 <- round(sum(weekly_series_2020$deaths) - sum(boot20$lwr), 0)
    
    ci_2020 <- paste0(fit_2020, " (", UB_2020, "; ", LB_2020, ")")
    
    
    
    # Prediction Interval 2021----
    
    weekly_series_2021 = data_under_review %>%
      filter(., year == 2021)
    Y1 <- Y2 <- matrix(NA,nrow = 52,ncol=8)
    for (i in 1:8) {
      Y1[,i] <- sin(2*i*pi*(469:520)/52)
      Y2[,i] <- cos(2*i*pi*(469:520)/52)
    }
    data.Y <- data.frame(Y1,Y2,rep(2019,52))
    names(data.Y) <-  c("sin.1","sin.2","sin.3","sin.4","sin.5","sin.6","sin.7","sin.8",
                        "cos.1","cos.2","cos.3","cos.4","cos.5","cos.6","cos.7","cos.8","year")
    
    boot21 <- predictInterval(merMod = model, newdata = data.Y,
                              level = 0.95, n.sims = 2000,
                              stat = "median", type="probability",
                              include.resid.var = FALSE)
    fit_2021 <- round(sum(weekly_series_2021$deaths) - sum(boot21$fit), 0)
    mean_15_19 <- sum(boot21$fit)
    UB_2021 <- round(sum(weekly_series_2021$deaths) - sum(boot21$upr), 0)
    LB_2021 <- round(sum(weekly_series_2021$deaths) - sum(boot21$lwr), 0)
    
    ci_2021 <- paste0(fit_2021, " (", UB_2021, "; ", LB_2021, ")")
    
    data_epi_plot = data.frame(week = c(seq(ymd('2020-01-01'),
                                            ymd('2020-12-31'),
                                            by = '1 week')[1:52],
                                        seq(ymd('2021-01-01'),
                                            ymd('2021-12-31'),
                                            by = '1 week')[2:53]),
                               real = c(weekly_series_2020$deaths,
                                        weekly_series_2021$deaths),
                               baseline = c(boot20$fit,
                                            boot21$fit),
                               LB = c(boot20$lwr,
                                      boot21$lwr),
                               UB = c(boot20$upr,
                                      boot21$upr))
    
    epi_plot = ggplot(data_epi_plot,
                      aes(x = week))+
      geom_path(aes(y = real, color = "Number of deaths"), lwd=.8)+
      geom_line(aes(y = baseline, color = "Expected number of deaths"), lwd=.8)+
      geom_ribbon(aes(ymin=LB,ymax=UB),alpha=0.3)+
      theme_epicurve()+
      scale_color_manual(name =" ",
                         values = c(`Number of deaths` = "Red",
                                    `Expected number of deaths` = "Blue"))+
      scale_x_date(date_breaks = "1 months", date_labels = "%b %Y") +
      labs(y = "Weekly number of deaths",
           caption = paste0(label_area, " ", label_age))
    
    
    return(list(model,
                `Plot baseline` = plot_baseline,
                mean_15_19 = mean_15_19,
                `Excess 2020` = ci_2020,
                fit_2020 = fit_2020,
                UB_2020 = UB_2020,
                LB_2020 = LB_2020,
                `Excess 2021` = ci_2021,
                fit_2021 = fit_2021,
                UB_2021 = UB_2021,
                LB_2021 = LB_2021,
                `Plot excess 20-21` = epi_plot,
                label = paste0(label_area, " Italia, ", label_age),
                area = label_area,
                age = label_age))
  }
  
}


# Estimate Italy and cross-validation analysis ----------------------------
dataset = mortality_data %>%
  group_by(year, week) %>%
  summarise(deaths = sum(deaths))


mod_Maruotti_Italy <- stima_eccesso_fne(mortality_data %>%
                                          group_by(year, week) %>%
                                          summarise(deaths = sum(deaths)),
                                        baseline = 2011:2019,
                                        label_area = "All Italy",
                                        label_age = "All ages",
                                        cv = F)
mod_Maruotti_Italy$`Plot baseline`
mod_Maruotti_Italy$`Excess 2020`
mod_Maruotti_Italy$`Excess 2021`
mod_Maruotti_Italy$`Plot excess 20-21`

RMSE_Italy_Maruotti <- vector()
RMSE_Italy_Maruotti <- stima_eccesso_fne(mortality_data %>%
                                           group_by(year, week) %>%
                                           summarise(deaths = sum(deaths)),
                                         baseline = 2011:2019,
                                         label_area = "All Italy",
                                         label_age = "All ages",
                                         cv = T)
mean(RMSE_Italy_Maruotti)



# Excess area and age -----------------------------------------------------


area = c("North Italy", "Center Italy", "South Italy and islands", "Italy")
age = c("0-49", "50-64", "65-79", "80+", "All ages")

results_area_age_Maruotti <- matrix(NA, 20, 4)
rownames(results_area_age_Maruotti) <- paste0(rep(area, each = 5),
                                              " ",
                                              rep(age, 4))
colnames(results_area_age_Maruotti) <- c("2020", "2021", "% 2020", "% 2021")

list_models_area_age = list()
z = 1
for(i in area){
  for(j in age){
    
    list_models_area_age[[z]] = if(i == "Italy" & j == "All ages"){
      stima_eccesso_fne(mortality_data %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2011:2019,
                        label_area = "All Italy",
                        label_age = "All ages",
                        cv = F)
    }else if(i == "Italy" & j != "All ages"){
      stima_eccesso_fne(mortality_data %>%
                          filter(age_class == j) %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2011:2019,
                        label_area = "All Italy",
                        label_age = j,
                        cv = F)
    }else if(i != "Italy" & j == "All ages"){
      stima_eccesso_fne(mortality_data %>%
                          filter(area == i) %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2011:2019,
                        label_area = i,
                        label_age = "All ages",
                        cv = F)
    }else {
      stima_eccesso_fne(mortality_data %>%
                          filter(area == i & age_class == j) %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2011:2019,
                        label_area = i,
                        label_age = j,
                        cv = F)
    }
    z = z + 1
    
  }
}


for (i in 1:20) {
  results_area_age_Maruotti[i,1] = list_models_area_age[[i]]$`Excess 2020`
  results_area_age_Maruotti[i,2] = list_models_area_age[[i]]$`Excess 2021`
  results_area_age_Maruotti[i,3] = round(list_models_area_age[[i]]$fit_2020/list_models_area_age[[i]]$mean_15_19, 4)*100
  results_area_age_Maruotti[i,4] = round(list_models_area_age[[i]]$fit_2021/list_models_area_age[[i]]$mean_15_19, 4)*100
}
results_area_age_Maruotti




# Excess region -----------------------------------------------------------


results_region_Maruotti <- matrix(NA, 21, 4)
rownames(results_region_Maruotti) <- unique(mortality_data$region)
colnames(results_region_Maruotti) <- c("2020", "2021", "% 2020", "% 2021")

list_models_region = list()

z = 1
for(i in unique(mortality_data$region)){
  
  list_models_region[[z]] =  stima_eccesso_fne(mortality_data %>%
                                                 filter(region == i) %>%
                                                 group_by(year, week) %>%
                                                 summarise(deaths = sum(deaths)),
                                               baseline = 2011:2019,
                                               label_area = i,
                                               label_age = "All ages",
                                               cv = F)
  z = z + 1
  
}


for (i in 1:21) {
  results_region_Maruotti[i,1] = list_models_region[[i]]$`Excess 2020`
  results_region_Maruotti[i,2] = list_models_region[[i]]$`Excess 2021`
  results_region_Maruotti[i,3] = round(list_models_region[[i]]$fit_2020/list_models_region[[i]]$mean_15_19, 4)*100
  results_region_Maruotti[i,4] = round(list_models_region[[i]]$fit_2021/list_models_region[[i]]$mean_15_19, 4)*100
}
results_region_Maruotti



# Excess area and age change of baseline ----------------------------------


results_area_age_Maruotti_1519 <- matrix(NA, 20, 4)
rownames(results_area_age_Maruotti_1519) <- paste0(rep(area, each = 5),
                                              " ",
                                              rep(age, 4))
colnames(results_area_age_Maruotti_1519) <- c("2020", "2021", "% 2020", "% 2021")

list_models_area_age_1519 = list()
z = 1
for(i in area){
  for(j in age){
    
    list_models_area_age_1519[[z]] = if(i == "Italy" & j == "All ages"){
      stima_eccesso_fne(mortality_data %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2015:2019,
                        label_area = "All Italy",
                        label_age = "All ages",
                        cv = F)
    }else if(i == "Italy" & j != "All ages"){
      stima_eccesso_fne(mortality_data %>%
                          filter(age_class == j) %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2015:2019,
                        label_area = "All Italy",
                        label_age = j,
                        cv = F)
    }else if(i != "Italy" & j == "All ages"){
      stima_eccesso_fne(mortality_data %>%
                          filter(area == i) %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2015:2019,
                        label_area = i,
                        label_age = "All ages",
                        cv = F)
    }else {
      stima_eccesso_fne(mortality_data %>%
                          filter(area == i & age_class == j) %>%
                          group_by(year, week) %>%
                          summarise(deaths = sum(deaths)),
                        baseline = 2015:2019,
                        label_area = i,
                        label_age = j,
                        cv = F)
    }
    z = z + 1
    
  }
}


for (i in 1:20) {
  results_area_age_Maruotti_1519[i,1] = list_models_area_age_1519[[i]]$`Excess 2020`
  results_area_age_Maruotti_1519[i,2] = list_models_area_age_1519[[i]]$`Excess 2021`
  results_area_age_Maruotti_1519[i,3] = round(list_models_area_age_1519[[i]]$fit_2020/list_models_area_age_1519[[i]]$mean_15_19, 4)*100
  results_area_age_Maruotti_1519[i,4] = round(list_models_area_age_1519[[i]]$fit_2021/list_models_area_age_1519[[i]]$mean_15_19, 4)*100
}
results_area_age_Maruotti_1519


