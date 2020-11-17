#### Functions for RLE Scripts 
## EClausius November 2020 

#Calculate IDW weights 
IDW_indicators <- function(varname, xydist, p, dat){
  indicator_values <- pull(dat, varname)
  indicator_values <- na.omit(indicator_values)
  w <- 1/(xydist^p)
  indic <- matrix(indicator_values, nrow = 1)
  isreal <- which(!is.na(indic))
  nvals <- length(isreal)
  
  #unelegant way 
  # val <- matrix(rep(indic[-isnan], ncol(w)), ncol = ncol(w), byrow = F) *
  #   w[-isnan,]
  # val <- colSums(val)/colSums(w)
  #elegant way using matrix multiply
  val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
  as.numeric(val)
}

#test

#Get HEX values as a quantile, given a dataframe
gethex <- function(varname, dat, probs = seq(0, 1, by =0.05), pal = "RdBu"){
  #pal can be an RColorBrewer palette name or a vector of hex codes
  indicator_values <- pull(dat, varname)
  #leaflet::colorNumeric() 
  xcols <- leaflet::colorQuantile(pal, 
                                  indicator_values, reverse = TRUE,
                                  probs = probs)
  xcols(indicator_values)
}



# Fits a gam for a given location 
fitgam <- function(this_location, datin, varname, save_plots){
  # this_location = "Rottnest Island"
  # datin <- metrics_location
  # varname <- "CTI"
  form <- eval(parse(text = 
                       paste0(varname, " ~ s(year_month)")))
  datsub <- filter(datin, Location == this_location)
  datsub$SiteCode <- factor(datsub$SiteCode)
  xout <- NULL
  
  xout <- try(
    m1 <- gamm(form, method = "REML", 
               random = list(SiteCode= ~ 1),
               data = datsub), TRUE)
  
  if (class(xout) != "try-error"){
    png(filename = paste0("outputs/plots-gams/", varname, "/", this_location, ".png"))
    plot(xout$gam, main = paste0(this_location, varname))
    dev.off()
    
    #
    # Predictions 
    #
    
    #Regional trend
    preds <- with(datsub, data.frame(year_month = 
                                       seq(min(year_month), max(year_month), 
                                           length.out = 50), 
                                     SiteCode = NA))
    predmod <- predict(m1$gam, newdata = preds, se = TRUE)
    preds$fit <- predmod$fit
    preds$se <- predmod$se.fit
    preds$Location <- this_location
    preds$metric <- varname
    preds$fit <- ifelse(preds$fit < 0, 0, as.numeric(preds$fit))    ###change any values below 0 to 0 to make more ecological sense 
    
    #Site by site trend
    predsites <- with(datsub, expand.grid(year_month = 
                                            seq(min(Year), max(Year), 
                                                by = 1), 
                                          SiteCode = factor(unique(SiteCode))))
    #predict GAM component
    predmodsites <- predict(m1$gam, newdata = predsites, se = FALSE)
    #Extract random intercepts
    re <- coef(m1$lme)[ncol(coef(m1$lme))]
    # add the random intercepts back on to the global mean
    predsites$fit <- predmodsites + re[[1]][match(predsites$SiteCode,
                                                  gsub(".*/", "", rownames(re)) )]
    xtemp <- datsub %>% dplyr::select(SiteCode, year_month = Year) %>% distinct() %>%
      mutate(has_survey = "yes")
    ztemp <- datsub %>% dplyr::select(SiteCode, Lat, Lon) %>% distinct()
    
    spatial_pred <- predsites %>% 
      full_join(xtemp, by = c("year_month", "SiteCode")) %>% 
      full_join(ztemp, by = c("SiteCode")) %>%
      mutate(has_survey = ifelse(is.na(has_survey), "no", has_survey))
    spatial_pred$Location <- this_location
    spatial_pred$metric <- varname
    spatial_pred$fit <- ifelse(spatial_pred$fit <0 , 0, as.numeric(spatial_pred$fit))
    
    xsave <- list(model = xout, predictions = preds, site_predictions = spatial_pred)
    
    #
    # Plots 
    #
    if(save_plots){
      g1 <- ggplot() + 
        geom_point(data = datsub, 
                   aes_string(x = "year_month", y = varname), 
                   alpha = 0.15) + 
        geom_line(data = preds, aes(x = year_month, y = fit)) + 
        geom_ribbon(data = preds, aes(x = year_month, 
                                      ymin = fit -se, ymax = fit + se), 
                    alpha = 0.5, fill = "tomato") + 
        theme_bw() + 
        xlab("Year") + 
        ylab(varname)
      
      g2 <- ggplot() + 
        geom_line(data = preds, aes(x = year_month, y = fit)) + 
        geom_ribbon(data = preds, aes(x = year_month, 
                                      ymin = fit -se, ymax = fit + se), 
                    alpha = 0.5, fill = "tomato") + 
        theme_bw() + 
        xlab("Year") + 
        ylab(varname)
      
      residdf <- data.frame(fit = xout$gam$fitted, resid = xout$gam$residuals)
      
      g3 <- ggplot(residdf) + 
        geom_point(aes(x = fit, y = resid)) + 
        xlab("Fitted") + 
        ylab("Residuals")
      
      gall <- cowplot::plot_grid(g1, g2, g3, nrow = 1)
      
      g4 <- ggplot() + 
        geom_line(data = spatial_pred, aes(x = year_month, y = fit)) + 
        geom_point(data = spatial_pred, aes(x = year_month, 
                                            y = fit, color = has_survey), 
                   size = 1) + 
        geom_point(data = datsub, aes_string(x = "year_month", y = varname), 
                   alpha = 0.5, color = "black") + 
        facet_wrap(~SiteCode) + 
        theme_bw()
      
      g5 <- ggplot() + 
        geom_line(data = preds, aes(x = year_month, y = fit)) + 
        geom_ribbon(data = preds, aes(x = year_month, 
                                      ymin = fit - se, ymax = fit + se), 
                    alpha = 0.5, fill = "tomato") + 
        theme_bw() + 
        xlab("Year") + 
        ylab(varname) + 
        ylim(0, NA)
      
      
      ggsave(gall, width = 10, height = 3,
             filename = paste0("outputs/plots-gams/", varname, "/", this_location, ".png"))
      ggsave(g4, width = 10, height = 10, 
             filename = paste0("outputs/plots-gams/", varname, "/", this_location,
                               "_SiteRandom.png"))
      ggsave(g2, width = 10, height = 1, 
             filename = paste0("outputs/plots-gams/", varname, "/", this_location,
                               "_elongated_plot.png"))
      ggsave(g5, width = 10, height = 1, 
             filename = paste0("outputs/plots-gams/", varname, "/", this_location,
                               "_zero-vmax.png"))
      rm(g1, g2, g3, g4, g5, gall)
    }
    rm(datsub, residdf, preds, predsites, spatial_pred, xout)
    
    
  } else{
    
    xsave <- xout
    rm(datsub, xout)
  }
  
  xsave
}