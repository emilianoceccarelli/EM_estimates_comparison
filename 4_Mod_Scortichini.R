
# Data preparation and parameters setting ---------------------------------


load("data/regioni_temp_medie.RData")


temp_media_reg <- temp_media_reg %>%
  mutate(COD_REG21 = as.numeric(COD_REG21))
mortality_data_regionale <- read.table("data/c_index.csv", sep = ",", header = T, fileEncoding = "Latin1") %>%
  dplyr::select("PRO_COM_T", "COD_REG21")%>%
  left_join(., read.csv("data/comuni_giornaliero_31dicembre21.csv", sep = ","),
            by = c("PRO_COM_T" = "COD_PROVCOM"))%>%
  dplyr::select(COD_REG21, GE, CL_ETA, T_11, T_12, T_13, T_14, T_15, T_16, T_17, T_18, T_19, T_20, T_21) %>%
  mutate(across(everything(), ~replace(., . ==  "n.d." , 0)),
         across(T_11:T_21, ~as.numeric(.)),
         CL_ETA = ifelse(CL_ETA %in% 0:10,
                         "0-49",
                         ifelse(CL_ETA %in% 11:13,
                                "50-64",
                                ifelse(CL_ETA %in% 14:16,
                                       "65-79",
                                       "80+")))) %>%
  group_by(COD_REG21, CL_ETA, GE) %>%
  summarise(across(T_11:T_21, list(sum)))%>%
  `colnames<-`(c("COD_REG21", "CL_ETA", "Data", paste0("T", 2011:2021))) %>%
  mutate(Data = ifelse(Data==229, 228, Data))%>% #sostituisco il 29/02 con il 28/02 per non avere problemi
  mutate(Data = as.character(Data))%>%
  mutate(Data = ifelse(stringr::str_count(Data)==3,
                       paste0("0", Data),
                       Data)) %>% #Aggiungo uno 0 davanti alla data per poterla trasformare in formato Date
  group_by(COD_REG21, CL_ETA, Data)%>%
  summarise(across(T2011:T2021, list(sum)))

mortality_data_regionale <- bind_cols(data = as.Date(paste0(rep(c(2011:2021),
                                                                nrow(mortality_data_regionale)),
                                                            rep(mortality_data_regionale$Data,
                                                                each = 11)),
                                                     "%Y%m%d"),
                                      COD_REG21 = rep(mortality_data_regionale$COD_REG21,
                                                      each = 11),
                                      classe_eta = rep(mortality_data_regionale$CL_ETA,
                                                       each = 11),
                                      decessi = c(t(as.matrix(mortality_data_regionale[,4:14]))))%>% #Trasformo la data in Date e incolonno 
  arrange(data) %>%
  group_by(COD_REG21, data, classe_eta) %>%
  summarise(decessi = sum(decessi)) %>%
  mutate(area = ifelse(COD_REG21 %in% c(10:80),
                       "North Italy",
                       ifelse(COD_REG21 %in% c(90:120),
                              "Center Italy",
                              "South Italy and islands")))

mortality_data_regionale = mortality_data_regionale%>%
  left_join(., temp_media_reg, by = c("COD_REG21", "data")) 


datamodel2 = mortality_data_regionale %>%
  filter(data>="2015-01-01")%>%
  group_by(COD_REG21, data) %>%
  summarise(decessi = sum(decessi),
            temp = mean(temp)) %>%
  mutate(y = decessi)


# DEFINE START AND END DAY FOR POST-PERIOD AND COVID PERIOD
startdate <- dmy(01012020)
coviddate <- dmy(31122019)
enddate <- dmy(31122021)

# DEFINE WEEK PERIODS FOR FEB-APR
seqpost <- seq(startdate, enddate, 1)
#cutdate <- unique(c(startdate-1,seqpost[seqpost %in% tapply(seqpost, week(seqpost), last)]))
cutdate <- as.Date(c("2019-12-31", "2020-12-31", "2021-12-31"),
                   format = "%Y-%m-%d")
labperiod1 <- sapply(seq(length(cutdate)-1), function(i) 
  paste(paste0(day(cutdate[i]+1), month(cutdate[i]+1,lab=T)),
        paste0(day(cutdate[i+1]), month(cutdate[i+1],lab=T)), sep="-"))
labperiod2 <- paste("year", unique(year(seqpost)))
datamodel2$tspost <- pmax(as.numeric(datamodel2$data-startdate),0)

# DEFINE PERIODS
seqperiod <- cut(unique(datamodel2$tspost), cutdate-startdate, labels=labperiod2)


mformula <- y ~ bpost + data + bseas + factor(wday(data)) + cbtmean

# DEFINE THE KNOTS FOR THE SPLINES FOR MODELLING EXCESS IN POST-PERIOD
nkpost <- 4

# DEFINE THE DF FOR THE CYCLIC SPLINE FOR SEASONALITY
dfseas <- 5

# DEFINE PARAMETERS OF CROSS-BASIS FOR TEMPERATURE
lagtmean <- 21
kpertmean <- c(10,75,90)
nklagtmean <- 3


seqreg <- unique(datamodel2$COD_REG21)
lab_reg <- read.csv("data/COD_REG.csv",
                    sep = ",", header = T, fileEncoding = "Latin1") %>%
  arrange(COD_REG21) %>%
  .$NOME_NUTS2
#labprov <- sapply(strsplit(unique(datamodel2$regione), "/"), "[[", 1)


# Excess regions  --------------------------------------------------------


# LIST TO STORE COEF/VCOV AND CONVERGENCE INDICATOR
stage1list <- vector("list", length(seqreg))
names(stage1list) <- lab_reg

#i = 1
for(i in seq(seqreg)) {
  
  # PRINT
  cat(lab_reg[i],"")
  
  # EXTRACT THE DATA AND COMPLETE IT
  dd <- subset(datamodel2, COD_REG21 ==seqreg[i])
  
  # DEFINE BASIS FUNCTIONS FOR POST-PERIOD, SEASONALITY, TEMPERATURE, FLU
  # NB: USE onebasis TO SIMPLIFY PREDICTIONS AND PLOTTING
  kpost <- equalknots(dd$tspost, nkpost)
  bpost <- onebasis(dd$tspost, fun="bs", degree=2, knots=kpost)
  kseas <- equalknots(yday(dd$data), dfseas)
  bseas <- onebasis(yday(dd$data), fun="pbs", knots=kseas)
  cbtmean <- crossbasis(dd$temp, lag=lagtmean,
                        argvar=list(fun="bs", degree=2, knots=quantile(dd$temp, kpertmean/100)),
                          arglag=list(knots=logknots(lagtmean, nklagtmean)))
  #flu013 <- runMean(dd$flu, 0:13)
  
  # RUN THE MODEL
  # NB: PRESERVE MISSING TO COMPUTE RESIDUALS LATER
  mod <- glm(mformula, data=dd, family=quasipoisson, na.action="na.exclude")
  
  # SAVE THE RESULTS: COEF/VCOV, RESIDUALS, OVERDISPERSION
  loglik <- sum(dpois(mod$y,mod$fitted.values,log=TRUE))
  disp <- sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.res
  stage1list[[i]] <- list(coef=coef(mod), vcov=vcov(mod), dispersion=disp,
                          residuals=residuals(mod, type="deviance"))
}


summary(mod)
# SECOND-STAGE META-ANALYSES

# MULTIVARIATE META-ANALYSIS OF COEFFICIENTS OF POST-PERIOD EXCESS
indpost <- grep("bpost", names(stage1list[[1]]$coef))
coefpost <- t(sapply(stage1list, function(x) x$coef[indpost]))
Scov <- lapply(stage1list, function(x) x$vcov[indpost,indpost])
metapost <- mixmeta(coefpost, Scov)
bluppost <- blup(metapost, vcov=T)


# COMPUTE EXCESS MORTALITY

# REDEFINE BASIS 
kpost <- equalknots(datamodel2$tspost, nkpost)
bpost <- onebasis(unique(datamodel2$tspost), fun="bs", degree=2, knots=kpost)

# DEFINE ARRAY TO STORE THE EXCESS DEATHS BY PROVINCE, PERIOD, RESAMPLING
nsim = 1000
excprovsim <- array(NA, dim=c(length(seqreg), length(labperiod1)+1, nsim+1),
                    dimnames=list(lab_reg, c("1gen20-31dec21",labperiod1),
                                  c("est",paste0("sim",seq(nsim)))))

# LOOP ACROSS PROVINCES
for(i in seq(seqreg)) {
  
  # PRINT
  cat(lab_reg[i],"")
  
  # RETRIEVE COEF/VCOV AND TOTAL DEATHS
  coef <- bluppost[[i]]$blup
  vcov <- bluppost[[i]]$vcov
  death <- subset(datamodel2, COD_REG21==seqreg[i] & data>=startdate)$decessi
  
  # COMPUTE ATTRIBUTABLE NUMBER (EXCESS), AND STORE THE SUM BY PERIOD
  an <- (1-exp(-bpost%*%coef))*death
  indcovid <- seqpost>=coviddate
  excprovsim[i,1,"est"] <- sum(an[indcovid], na.rm = T)
  excprovsim[i,-1,"est"] <- tapply(an, seqperiod, sum)
  
  # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
  set.seed(13041975)
  coefsim <- mvrnorm(nsim, coef, vcov)
  
  # LOOP ACROSS ITERATIONS AND DO AS ABOVE WITH RESAMPLES COEF
  for(s in seq(nsim)) {
    an <- (1-exp(-bpost%*%coefsim[s,]))*death
    excprovsim[i,1,s+1] <- sum(an[indcovid])
    excprovsim[i,-1,s+1] <- tapply(an, seqperiod, sum)
  }
}

# COLLAPSE BY REGION AND THEN FULL COUNTRY
excregsim <- apply(excprovsim, 2:3, tapply, seqreg, sum)
excitalysim <- apply(excregsim, 2:3, sum)

risultati_regione_iniz = cbind(excregsim[,c(2,3),"est"],lab_reg) %>%
  `colnames<-`(c("2020", "2021", "Regione")) %>%
  as.data.frame() %>%
  mutate(`2020` = round(as.numeric(`2020`), 0),
         `2021` = round(as.numeric(`2021`), 0)) %>%
  arrange(Regione) %>%
  tibble::rownames_to_column(., "COD_REG21")

IC_2020 <- apply(excregsim[,2,-1], 1, quantile, c(2.5,97.5)/100) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "COD_REG21")%>%
  mutate(`2.5%` = round(as.numeric(`2.5%`), 0),
         `97.5%` = round(as.numeric(`97.5%`), 0)) 
IC_2021 <- apply(excregsim[,3,-1], 1, quantile, c(2.5,97.5)/100)%>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "COD_REG21")%>%
  mutate(`2.5%` = round(as.numeric(`2.5%`), 0),
         `97.5%` = round(as.numeric(`97.5%`), 0))


risultati_regione <- left_join(risultati_regione_iniz, 
                               IC_2020,
                               by = "COD_REG21") %>%
  mutate(`2020` = paste0(`2020`, " (", `2.5%`, ", ", `97.5%`, ")")) %>%
  dplyr::select(, -c(5,6))%>%
  left_join(.,
            IC_2021,
            by = "COD_REG21") %>% 
  mutate(`2021` = paste0(`2021`, " (", `2.5%`, ", ", `97.5%`, ")")) %>%
  dplyr::select(, -c(5,6))


eccesso_rel_regioni <- risultati_regione_iniz %>%
  left_join(., mortality_data_regionale %>%
              mutate(anno = year(data)) %>%
              filter(anno>2019) %>%
              group_by(COD_REG21, anno)%>%
              summarise(decessi = sum(decessi)) %>%
              pivot_wider(names_from = anno, values_from = decessi) %>%
              `colnames<-`(c("COD_REG21", "2020_tot", "2021_tot")) %>%
              mutate(COD_REG21 = as.character(COD_REG21)),
            by = "COD_REG21") %>%
  group_by(COD_REG21) %>%
  mutate(`% 2020` = as.character(round(`2020`/(`2020_tot`-`2020`), 4)*100),
         `% 2021` = as.character(round(`2021`/(`2021_tot`-`2021`), 4)*100))%>%
  as.data.frame()

results_region_Scortichini <- left_join(risultati_regione %>%
                                          dplyr::select(c(`2020`, `2021`, Regione)), 
                                        eccesso_rel_regioni %>%
                                          dplyr::select(c(`% 2020`, `% 2021`, Regione)),
                                        by = "Regione") %>%
  tibble::column_to_rownames(var="Regione")
results_region_Scortichini



# Excess area and age class -----------------------------------------------
area = c("North Italy", "Center Italy", "South Italy and islands", "Italy")
age = c("0-49", "50-64", "65-79", "80+", "All ages")
seq_age_rip = 1:20
labage_rip = paste0(rep(area, each = 5),
                    " ",
                    rep(age, 4))
stage1list_age_rip = vector("list", length(seq_age_rip))
names(stage1list_age_rip) <- labage_rip




z = 1
for(i in area){
  for(j in age){
    dd = data.frame()
    
    if(i == "Italy" & j == "All ages"){
      dd <- mortality_data_regionale %>%
        filter(data>="2015-01-01") %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }else if(i == "Italy" & j != "All ages"){
      dd <- mortality_data_regionale%>%
        filter(classe_eta == j, data>="2015-01-01") %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }else if(i != "Italy" & j == "All ages"){
      dd <- mortality_data_regionale%>%
        filter(area == i, data>="2015-01-01") %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }else{
      dd <- mortality_data_regionale %>%
        filter(area == i, classe_eta == j, data>="2015-01-01")%>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }
    
    
    cat(labage_rip[z],"")
    
    
    # DEFINE BASIS FUNCTIONS FOR POST-PERIOD, SEASONALITY, TEMPERATURE, FLU
    # NB: USE onebasis TO SIMPLIFY PREDICTIONS AND PLOTTING
    kpost <- equalknots(dd$tspost, nkpost)
    bpost <- onebasis(dd$tspost, fun="bs", degree=2, knots=kpost)
    kseas <- equalknots(yday(dd$data), dfseas)
    bseas <- onebasis(yday(dd$data), fun="pbs", knots=kseas)
    cbtmean <- crossbasis(dd$temp, lag=lagtmean,
                          argvar=list(fun="bs", degree=2, knots=quantile(dd$temp, kpertmean/100)),
                          arglag=list(knots=logknots(lagtmean, nklagtmean)))
    #flu013 <- runMean(dd$flu, 0:13)
    
    # RUN THE MODEL
    # NB: PRESERVE MISSING TO COMPUTE RESIDUALS LATER
    mod <- glm(mformula, data=dd, family=quasipoisson, na.action="na.exclude")
    
    # SAVE THE RESULTS: COEF/VCOV, RESIDUALS, OVERDISPERSION
    loglik <- sum(dpois(mod$y,mod$fitted.values,log=TRUE))
    disp <- sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.res
    stage1list_age_rip[[z]] <- list(coef=coef(mod), vcov=vcov(mod), dispersion=disp,
                                    residuals=residuals(mod, type="deviance"))
    z = z + 1
  }
}



# MULTIVARIATE META-ANALYSIS OF COEFFICIENTS OF POST-PERIOD EXCESS
indpost <- grep("bpost", names(stage1list_age_rip[[1]]$coef))
coefpost <- t(sapply(stage1list_age_rip, function(x) x$coef[indpost]))
Scov <- lapply(stage1list_age_rip, function(x) x$vcov[indpost,indpost])
metapost <- mixmeta(coefpost, Scov)
bluppost <- blup(metapost, vcov=T)


# COMPUTE EXCESS MORTALITY

# REDEFINE BASIS 
kpost <- equalknots(datamodel2$tspost, nkpost)
bpost <- onebasis(unique(datamodel2$tspost), fun="bs", degree=2, knots=kpost)

# DEFINE ARRAY TO STORE THE EXCESS DEATHS BY PROVINCE, PERIOD, RESAMPLING
exc_area_age_sim <- array(NA, dim=c(length(seq_age_rip), length(labperiod1)+1, nsim+1),
                          dimnames=list(labage_rip, c("1gen20-31dec21",labperiod1),
                                        c("est",paste0("sim",seq(nsim)))))

sim_italy = matrix(NA, 730,1000)


z = 1
for(i in area){
  for(j in age){
    
    cat(labage_rip[z],"")
    # RETRIEVE COEF/VCOV AND TOTAL DEATHS
    coef <- bluppost[[z]]$blup
    vcov <- bluppost[[z]]$vcov
    if(i == "Italy" & j == "All ages"){
      death <- mortality_data_regionale %>%
        filter(data>="2015-01-01") %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi)) %>%
        filter(data>=startdate) %>%
        .$decessi
    }else if(i == "Italy" & j != "All ages"){
      death <- mortality_data_regionale%>%
        filter(data>="2015-01-01") %>%
        filter(classe_eta == j) %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi))%>%
        filter(data>=startdate) %>%
        .$decessi
    }else if(i != "Italy" & j == "All ages"){
      death <- mortality_data_regionale%>%
        filter(data>="2015-01-01") %>%
        filter(area == i) %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi)) %>%
        filter(data>=startdate) %>%
        .$decessi
    }else{
      death <- mortality_data_regionale %>%
        filter(area == i, classe_eta == j)%>%
        filter(data>="2015-01-01") %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi)) %>%
        filter(data>=startdate) %>%
        .$decessi
    }
    # COMPUTE ATTRIBUTABLE NUMBER (EXCESS), AND STORE THE SUM BY PERIOD
    an <- (1-exp(-bpost%*%coef))*death
    indcovid <- seqpost>=coviddate
    eccesso <- an[indcovid]
    exc_area_age_sim[z,1,"est"] <- sum(eccesso, na.rm = T)
    exc_area_age_sim[z,-1,"est"] <- tapply(an, seqperiod, sum)
    
    # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975)
    coefsim <- MASS::mvrnorm(nsim, coef, vcov)
    
    # LOOP ACROSS ITERATIONS AND DO AS ABOVE WITH RESAMPLES COEF
    for(s in seq(nsim)) {
      an <- (1-exp(-bpost%*%coefsim[s,]))*death
      exc_area_age_sim[z,1,s+1] <- sum(an[indcovid])
      exc_area_age_sim[z,-1,s+1] <- tapply(an, seqperiod, sum)
      
      if(z == 20){
        sim_italy[,s] = an
      }
      
    }
    z = z+1
  }
}

# COLLAPSE BY REGION AND THEN FULL COUNTRY
excregsim_areaage <- apply(exc_area_age_sim, 2:3, tapply, seq_age_rip, sum)
excitalysim <- apply(exc_area_age_sim, 2:3, sum)

risultati_rip_eta_iniz = cbind(excregsim_areaage[,c(2,3),"est"],labage_rip) %>%
  `colnames<-`(c("2020", "2021", "Rip")) %>%
  as.data.frame() %>%
  mutate(`2020` = round(as.numeric(`2020`), 0),
         `2021` = round(as.numeric(`2021`), 0))

IC_2020_rip <- apply(excregsim_areaage[,2,-1], 1, quantile, c(2.5,97.5)/100) %>%
  t() %>%
  as.data.frame() %>%
  mutate(`2.5%` = round(as.numeric(`2.5%`), 0),
         `97.5%` = round(as.numeric(`97.5%`), 0),
         Rip = labage_rip) 
IC_2021_rip <- apply(excregsim_areaage[,3,-1], 1, quantile, c(2.5,97.5)/100) %>%
  t() %>%
  as.data.frame() %>%
  mutate(`2.5%` = round(as.numeric(`2.5%`), 0),
         `97.5%` = round(as.numeric(`97.5%`), 0),
         Rip = labage_rip) 


risultati_rip_eta <- left_join(risultati_rip_eta_iniz, 
                               IC_2020_rip,
                               by = "Rip") %>%
  mutate(`2020` = paste0(`2020`, " (", `2.5%`, ", ", `97.5%`, ")")) %>%
  dplyr::select(, -c(4,5)) %>%
  left_join(.,
            IC_2021_rip,
            by = "Rip") %>% 
  mutate(`2021` = paste0(`2021`, " (", `2.5%`, ", ", `97.5%`, ")")) %>%
  dplyr::select(, -c(4,5))


eccesso_rel_rip_eta <-bind_rows(mortality_data_regionale %>%
                                  mutate(anno = year(data)) %>%
                                  filter(anno>2019) %>%
                                  group_by(area, classe_eta, anno)%>%
                                  summarise(decessi = sum(decessi)) %>%
                                  pivot_wider(names_from = anno, values_from = decessi),
                                mortality_data_regionale %>%
                                  mutate(anno = year(data)) %>%
                                  filter(anno>2019) %>%
                                  group_by(area, anno)%>%
                                  summarise(decessi = sum(decessi)) %>%
                                  pivot_wider(names_from = anno, values_from = decessi)%>%
                                  mutate(classe_eta = "Total") %>%
                                  relocate(.after = "area"),
                                mortality_data_regionale %>%
                                  mutate(anno = year(data)) %>%
                                  filter(anno>2019) %>%
                                  group_by(classe_eta, anno)%>%
                                  summarise(decessi = sum(decessi)) %>%
                                  pivot_wider(names_from = anno, values_from = decessi)%>%
                                  mutate(area = "Italy") %>%
                                  relocate(.before = "classe_eta"),
                                mortality_data_regionale %>%
                                  mutate(anno = year(data)) %>%
                                  filter(anno>2019) %>%
                                  group_by(anno)%>%
                                  summarise(decessi = sum(decessi)) %>%
                                  pivot_wider(names_from = anno, values_from = decessi)%>%
                                  mutate(area = "Italy",
                                         classe_eta = "Total") %>%
                                  relocate(classe_eta, .before = `2020`)%>%
                                  relocate(area, .before = `classe_eta`)) %>%
  mutate(area = ordered(area, levels = c("North Italy",
                                         "Center Italy",
                                         "South Italy and islands",
                                         "Italy")))%>%
  arrange(area,classe_eta) %>%
  group_by(area, classe_eta) %>%
  tibble::add_column(Rip = as.vector(labage_rip)) %>%
  `colnames<-`(c("area", "classe_eta", "2020_tot", "2021_tot", "Rip")) %>%
  left_join(., risultati_rip_eta_iniz,
            by = "Rip") %>%
  group_by(Rip) %>%
  mutate(`% 2020` = as.character(round(`2020`/(`2020_tot`-`2020`), 4)*100),
         `% 2021` = as.character(round(`2021`/(`2021_tot`-`2021`), 4)*100))%>%
  as.data.frame()

results_area_age_Scortichini <- left_join(risultati_rip_eta %>%
                                          dplyr::select(c(`2020`, `2021`, Rip)), 
                                          eccesso_rel_rip_eta %>%
                                          dplyr::select(c(`% 2020`, `% 2021`, Rip)),
                                        by = "Rip") %>%
  tibble::column_to_rownames(var="Rip")
results_area_age_Scortichini



eccesso_ITALIA = eccesso[1:730]
decessI_ITALIA_20_21 = mortality_data_regionale  %>%
  filter(data>=startdate) %>%
  group_by(data) %>%
  summarise(decessi = sum(decessi)) %>%
  .$decessi
baseline = decessI_ITALIA_20_21 - eccesso_ITALIA
LB_grafico = apply(matrix(decessI_ITALIA_20_21, 730, 1000)- sim_italy,
                   1,
                   quantile,
                   c(2.5,97.5)/100)[1,]
UB_grafico = apply(matrix(decessI_ITALIA_20_21, 730, 1000)- sim_italy,
                   1,
                   quantile,
                   c(2.5,97.5)/100)[2,]

# Plot excess 2020 2021 ---------------------------------------------------


library(ggplot2)
plot_baseline_Italy_Scortichini <- data.frame(data = seq(ymd("2020-01-01"),
                                                         ymd("2021-12-31"),
                                                         by = "1 day")[c(1:59,61:731)],
                                              baseline = baseline,
                                              decessi = decessI_ITALIA_20_21,
                                              LB = LB_grafico,
                                              UB = UB_grafico) %>%
  mutate(data = lubridate::floor_date(data, unit = "week", week_start = 3)) %>%
  group_by(data)%>%
  summarise(baseline = sum(baseline),
            decessi = sum(decessi),
            LB = sum(LB),
            UB = sum(UB)) %>%
  slice(-105) %>%
  ggplot(., 
         aes(x = data))+
  geom_path(aes(y = decessi, color = "Number of deaths"), lwd=.8)+
  geom_line(aes(y = baseline, color = "Expected number of deaths"), lwd=.8)+
  geom_ribbon(aes(ymin=LB,ymax=UB),alpha=0.3)+
  theme_epicurve()+
  scale_color_manual(name =" ",
                     values = c(`Number of deaths` = "Red",
                                `Expected number of deaths` = "Blue"))+
  ylab("Weekly number of deaths")+
  scale_x_date(date_breaks = "1 months", date_labels = "%b %Y",
               limits = c(ymd("2020-01-01"), ymd("2022-01-01")))
plot_baseline_Italy_Scortichini



# Change baseline to 2011-2019 area and age -------------------------------


stage1list_age_rip_1119 = vector("list", length(seq_age_rip))
names(stage1list_age_rip_1119) <- labage_rip


z = 1
for(i in area){
  for(j in age){
    dd = data.frame()
    
    if(i == "Italy" & j == "All ages"){
      dd <- mortality_data_regionale %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }else if(i == "Italy" & j != "All ages"){
      dd <- mortality_data_regionale%>%
        filter(classe_eta == j) %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }else if(i != "Italy" & j == "All ages"){
      dd <- mortality_data_regionale%>%
        filter(area == i) %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }else{
      dd <- mortality_data_regionale %>%
        filter(area == i, classe_eta == j)%>%
        group_by(data) %>%
        summarise(decessi = sum(decessi),
                  temp = mean(temp)) %>%
        mutate(tspost = pmax(as.numeric(.$data-startdate),0),
               y = decessi)
    }
    
    
    cat(labage_rip[z],"")
    
    
    # DEFINE BASIS FUNCTIONS FOR POST-PERIOD, SEASONALITY, TEMPERATURE, FLU
    # NB: USE onebasis TO SIMPLIFY PREDICTIONS AND PLOTTING
    kpost <- equalknots(dd$tspost, nkpost)
    bpost <- onebasis(dd$tspost, fun="bs", degree=2, knots=kpost)
    kseas <- equalknots(yday(dd$data), dfseas)
    bseas <- onebasis(yday(dd$data), fun="pbs", knots=kseas)
    cbtmean <- crossbasis(dd$temp, lag=lagtmean,
                          argvar=list(fun="bs", degree=2, knots=quantile(dd$temp, kpertmean/100)),
                          arglag=list(knots=logknots(lagtmean, nklagtmean)))
    #flu013 <- runMean(dd$flu, 0:13)
    
    # RUN THE MODEL
    # NB: PRESERVE MISSING TO COMPUTE RESIDUALS LATER
    mod <- glm(mformula, data=dd, family=quasipoisson, na.action="na.exclude")
    
    # SAVE THE RESULTS: COEF/VCOV, RESIDUALS, OVERDISPERSION
    loglik <- sum(dpois(mod$y,mod$fitted.values,log=TRUE))
    disp <- sum(residuals(mod,type="pearson")^2, na.rm=T)/mod$df.res
    stage1list_age_rip_1119[[z]] <- list(coef=coef(mod), vcov=vcov(mod), dispersion=disp,
                                         residuals=residuals(mod, type="deviance"))
    z = z + 1
  }
}

# MULTIVARIATE META-ANALYSIS OF COEFFICIENTS OF POST-PERIOD EXCESS
indpost <- grep("bpost", names(stage1list_age_rip_1119[[1]]$coef))
coefpost <- t(sapply(stage1list_age_rip_1119, function(x) x$coef[indpost]))
Scov <- lapply(stage1list_age_rip_1119, function(x) x$vcov[indpost,indpost])
metapost <- mixmeta(coefpost, Scov)
bluppost <- blup(metapost, vcov=T)


# COMPUTE EXCESS MORTALITY

# REDEFINE BASIS 
kpost <- equalknots(datamodel2$tspost, nkpost)
bpost <- onebasis(unique(datamodel2$tspost), fun="bs", degree=2, knots=kpost)

# DEFINE ARRAY TO STORE THE EXCESS DEATHS BY PROVINCE, PERIOD, RESAMPLING
exc_area_age_sim <- array(NA, dim=c(length(seq_age_rip), length(labperiod1)+1, nsim+1),
                          dimnames=list(labage_rip, c("1gen20-31dec21",labperiod1),
                                        c("est",paste0("sim",seq(nsim)))))



z = 1
for(i in area){
  for(j in age){
    
    cat(labage_rip[z],"")
    # RETRIEVE COEF/VCOV AND TOTAL DEATHS
    coef <- bluppost[[z]]$blup
    vcov <- bluppost[[z]]$vcov
    if(i == "Italy" & j == "All ages"){
      death <- mortality_data_regionale %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi)) %>%
        filter(data>=startdate) %>%
        .$decessi
    }else if(i == "Italy" & j != "All ages"){
      death <- mortality_data_regionale%>%
        filter(classe_eta == j) %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi))%>%
        filter(data>=startdate) %>%
        .$decessi
    }else if(i != "Italy" & j == "All ages"){
      death <- mortality_data_regionale%>%
        filter(area == i) %>%
        group_by(data) %>%
        summarise(decessi = sum(decessi)) %>%
        filter(data>=startdate) %>%
        .$decessi
    }else{
      death <- mortality_data_regionale %>%
        filter(area == i, classe_eta == j)%>%
        group_by(data) %>%
        summarise(decessi = sum(decessi)) %>%
        filter(data>=startdate) %>%
        .$decessi
    }
    # COMPUTE ATTRIBUTABLE NUMBER (EXCESS), AND STORE THE SUM BY PERIOD
    an <- (1-exp(-bpost%*%coef))*death
    indcovid <- seqpost>=coviddate
    exc_area_age_sim[z,1,"est"] <- sum(an[indcovid], na.rm = T)
    exc_area_age_sim[z,-1,"est"] <- tapply(an, seqperiod, sum)
    
    # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
    set.seed(13041975)
    coefsim <- mvrnorm(nsim, coef, vcov)
    
    # LOOP ACROSS ITERATIONS AND DO AS ABOVE WITH RESAMPLES COEF
    for(s in seq(nsim)) {
      an <- (1-exp(-bpost%*%coefsim[s,]))*death
      exc_area_age_sim[z,1,s+1] <- sum(an[indcovid])
      exc_area_age_sim[z,-1,s+1] <- tapply(an, seqperiod, sum)
    }
    z = z+1
  }
}

# COLLAPSE BY REGION AND THEN FULL COUNTRY
excregsim_areaage_1119 <- apply(exc_area_age_sim, 2:3, tapply, seq_age_rip, sum)
excitalysim <- apply(exc_area_age_sim, 2:3, sum)

risultati_rip_eta_iniz_11_19 = cbind(excregsim_areaage_1119[,c(2,3),"est"],labage_rip) %>%
  `colnames<-`(c("2020", "2021", "Rip")) %>%
  as.data.frame() %>%
  mutate(`2020` = round(as.numeric(`2020`), 0),
         `2021` = round(as.numeric(`2021`), 0))

IC_2020_rip_11_19 <- apply(excregsim_areaage_1119[,2,-1], 1, quantile, c(2.5,97.5)/100) %>%
  t() %>%
  as.data.frame() %>%
  mutate(`2.5%` = round(as.numeric(`2.5%`), 0),
         `97.5%` = round(as.numeric(`97.5%`), 0),
         Rip = labage_rip) 
IC_2021_rip_11_19 <- apply(excregsim_areaage_1119[,3,-1], 1, quantile, c(2.5,97.5)/100) %>%
  t() %>%
  as.data.frame() %>%
  mutate(`2.5%` = round(as.numeric(`2.5%`), 0),
         `97.5%` = round(as.numeric(`97.5%`), 0),
         Rip = labage_rip) 


risultati_rip_eta_11_19 <- left_join(risultati_rip_eta_iniz_11_19, 
                                     IC_2020_rip_11_19,
                                     by = "Rip") %>%
  mutate(`2020` = paste0(`2020`, " (", `2.5%`, ", ", `97.5%`, ")")) %>%
  dplyr::select(, -c(4,5)) %>%
  left_join(.,
            IC_2021_rip_11_19,
            by = "Rip") %>% 
  mutate(`2021` = paste0(`2021`, " (", `2.5%`, ", ", `97.5%`, ")")) %>%
  dplyr::select(, -c(4,5))
risultati_rip_eta_11_19
