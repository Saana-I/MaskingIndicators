  Sys.setenv(TZ='GMT')

##################################################################################
##### R tools

  set.seed(0)
  
  library(mgcv)
  
  ilogit <- function(x) {exp(x)/(exp(x)+1)}
  logit <- function(p) {log(p)-log(1-p)}
  R2 <- function(obs, fits) {
    SSE <- sum((obs-fits)^2)
    SST <- sum((obs-mean(obs))^2)
    return(1-SSE/SST)
  }
  
  myalpha <- 80
  greyt <- rgb(0, 0, 0, max = 255, alpha = myalpha)
  redt <- bluet <- rgb(255, 0, 0, max = 255, alpha = myalpha)
  bluet <- rgb(0, 0, 255, max = 255, alpha = myalpha)
  bluet2 <- rgb(100, 0, 255, max = 255, alpha = myalpha)

##################################################################################
##### Data 
  
  load("Supplement_data.Rd")

#################################################################################
##### Data transforms

  SELmax_base <- min(atab$SELmax, na.rm=T)-1
  atab$SELmax[is.na(atab$SELmax)] <- SELmax_base
  atab$SELmax <- atab$SELmax - SELmax_base
  
  atab$SELmax_CAS <- atab$SELmax
  atab$SELmax_CAS[atab$ExposureType=="MPS" | atab$ExposureType=="HPS"] <- 0
  
  atab$SELmax_PAS <- atab$SELmax
  atab$SELmax_PAS[atab$ExposureType=="CAS"] <- 0
  
  atab$SELmax2 <- atab$SELmax
  atab$SELmax2[atab$Angle > 90] <- 0
  
  atab$SELmax_CAS2 <- atab$SELmax_CAS
  atab$SELmax_CAS2[atab$Angle > 90] <- 0
  
  atab$SELmax_PAS2 <- atab$SELmax_PAS
  atab$SELmax_PAS2[atab$Angle > 90] <- 0
  
  atab$ID <- as.numeric(as.factor(atab$ind))
  
  atab$pitch_first2 <- (pi/2 - atab$pitch_first)/pi*180 
  atab$pitch_mean2 <- (pi/2 - atab$pitch_mean)/pi*180 
  
  atab$CR <- atab$ClickN/atab$duration
  atab$CR2 <- 1/atab$ICI_med
  atab$CR_log <- log(atab$CR2)
  atab$ICI_log <- log(atab$ICI_first)
  
  atab$duration_log <- log(atab$duration)
  atab$odba_log <- log(atab$odba_rms)
  
  atab$buzz <- as.numeric(atab$nextLabel=="BUZZ")
  
  atab$fluker <- atab$fluken/atab$duration
  atab$fluker_log <- log(atab$fluker)
  
  # Threshold for pauses from regular click trains
  tBool0 <- (atab$state=="2" | atab$state=="3" | atab$state=="4") & !is.na(atab$state) 
  tBool2 <- tBool0 & atab$label=="RC" & atab$ClickN > 5 & atab$CR <= 5
  tBool2[is.na(tBool2)] <- FALSE
  TH <- median(atab$ICI_first[tBool2], na.rm=T)
  
  atab$pause <- atab$pausedur > TH
  atab$pause[atab$pausedur < -TH] <- NA

  atab$ind_fac <- as.factor(atab$ind)
  
##################################################################################
##### Filters

  tBool0 <- (atab$state=="2" | atab$state=="3" | atab$state=="4") & !is.na(atab$state) 
  
  # Buzzes
  tBool1 <- tBool0 &  atab$label=="BUZZ"
  tBool1[is.na(tBool1)] <- FALSE
  sum(tBool1) # 1910
  
  bBool1 <- tBool1 & atab$ExposureType=="BP"
  eBool1 <- tBool1 & atab$SELmax > 0
  
  # Regular click trains
  tBool2 <- tBool0 & atab$label=="RC" & atab$ClickN > 5 & atab$CR <= 5
  tBool2[is.na(tBool2)] <- FALSE
  sum(tBool2) # 5963
  
  bBool2 <- tBool2 & atab$ExposureType=="BP"
  eBool2 <- tBool2 & atab$SELmax > 0
  
  excl2016 <- atab$GMTtime > as.POSIXct("2017-01-01 00:00:01 GMT")
  
  # Check sample sizes
  sum((atab$SELmax_CAS[tBool2]+SELmax_base)>=160) # 75
  sum((atab$SELmax_PAS[tBool2]+SELmax_base)>=160) # 27
  
  sum((atab$SELmax_CAS[tBool2]+SELmax_base)>=170) # 8
  sum((atab$SELmax_PAS[tBool2]+SELmax_base)>=170) # 2
  
  sum((atab$SeaState[tBool1])>=4, na.rm=T) # 129
  sum((atab$SELmax[tBool1]+SELmax_base)>=165) # 9
  
  sum((atab$SeaState[tBool2])>=4, na.rm=T) # 347
  sum((atab$SELmax[tBool2]+SELmax_base)>=165) # 37
  
  
##################################################################################
#####  STATISTICAL MODELLING
  
  # Function to show predictions during baseline
  pred_plots <- function(fit1, tBool, varname, log_var=FALSE) {
    
    mylims <- c(min(atab[tBool,varname], na.rm=T), max(atab[tBool,varname], na.rm=T))
    
    if(varname=="ICI_first" | varname=="duration") {  
      mylims <- c(min(atab[tBool,varname], na.rm=T), quantile(atab[tBool,varname], 0.97, na.rm=T))}
    
    if(varname=="duration_log" | varname=="odba_log" | varname=="CR_log") {  
      mylims <- c(min(exp(atab[tBool,varname]), na.rm=T), quantile(exp(atab[tBool,varname]), 0.97, na.rm=T))}
    
    par(mfrow=c(2,2), mar=c(5,4,1,1))
    
    preddata <- data.frame(SeaState=0:5)
    preddata$SELmax <- 0
    preddata$SELmax_CAS <- 0
    preddata$SELmax_PAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$Angle <- 90
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    preddata$pitch_first2 <- 90
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(atab$SeaState[tBool], y, 
         xlab="Sea state", ylab=varname, main="SeaState - Depth", col="darkgrey", ylim=mylims)
    for(j in seq(0,2000,100)) {
      preddata$depth_first <- j
      preddata$SeaState_TH <- preddata$SeaState
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(preddata$SeaState, y, col="blue")
    }
    lines(preddata$SeaState, y, col="cyan")
    
    preddata <- data.frame(SeaState=0:5)
    preddata$SELmax <- 0
    preddata$SELmax_CAS <- 0
    preddata$SELmax_PAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$Angle <- 90
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    preddata$pitch_first2 <- 90
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(atab$SeaState[tBool], y, 
         xlab="Sea state", ylab=varname, main="SeaState - V. Angle", col="darkgrey", ylim=mylims)
    for(j in seq(0,180,30)) {
      preddata$pitch_first2 <- j
      preddata$SeaState_TH <- preddata$SeaState
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(preddata$SeaState, y, col="blue")
    }
    lines(preddata$SeaState, y, col="cyan")
    
    preddata <- data.frame(depth_first=seq(0,2000,10))
    preddata$SELmax <- 0
    preddata$SELmax_CAS <- 0
    preddata$SELmax_PAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$Angle <- 90
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$Swell <- 0
    preddata$pitch_first2 <- 90
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(atab$depth_first[tBool], y, 
         xlab="Depth (m)", ylab=varname, main="Depth - SeaState", col="darkgrey", ylim=mylims)
    for(j in seq(0,5,0.5)) {
      preddata$SeaState <- j
      preddata$SeaState_TH <- preddata$SeaState
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(preddata$depth_first, y, col="blue")
    }
    lines(preddata$depth_first, y, col="cyan")
    
    preddata <- data.frame(pitch_first2=0:180)
    preddata$SELmax <- 0
    preddata$SELmax_CAS <- 0
    preddata$SELmax_PAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$Angle <- 90
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(atab$pitch_first2[tBool], y, 
         xlab="Vertical pitch (deg)", ylab=varname, main="V. Angle - SeaState", col="darkgrey", ylim=mylims)
    for(j in seq(0,5,0.5)) {
      preddata$SeaState <- j
      preddata$SeaState_TH <- preddata$SeaState
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(preddata$pitch_first2, y, col="blue")
    }
    lines(preddata$pitch_first2, y, col="cyan")
  }
  
  # Function to show predictions during exposures
  pred_plots2 <- function(fit1, tBool, varname, log_var=FALSE) {
    
    mylims <- c(min(atab[tBool,varname], na.rm=T), max(atab[tBool,varname], na.rm=T))
    
    if(varname=="ICI_first" | varname=="duration") {  
      mylims <- c(min(atab[tBool,varname], na.rm=T), quantile(atab[tBool,varname], 0.97, na.rm=T))}
    
    if(varname=="duration_log" | varname=="odba_log" | varname=="CR_log") {  
      mylims <- c(min(exp(atab[tBool,varname]), na.rm=T), quantile(exp(atab[tBool,varname]), 0.97, na.rm=T))}
    
    
    par(mfrow=c(2,2), mar=c(5,4,1,1))
    
    preddata <- data.frame(SELmax=0:107)
    preddata$SELmax_CAS <- preddata$SELmax
    preddata$SELmax_PAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$SeaState <- 3
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    preddata$pitch_first2 <- 45
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(SELmax_base + atab$SELmax[tBool], y, 
         xlab="SELmax CAS (dB)", ylab=varname, main="SELmax CAS - Angle", col="darkgrey", ylim=mylims)
    for(j in seq(0,180,10)) {
      preddata$Angle <- j
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(SELmax_base + preddata$SELmax, y, col="blue")
    }
    lines(SELmax_base + preddata$SELmax, y, col="cyan")
    
    preddata <- data.frame(SELmax=0:107)
    preddata$SELmax_PAS <- preddata$SELmax
    preddata$SELmax_CAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$SeaState <- 3
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    preddata$pitch_first2 <- 45
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(SELmax_base + atab$SELmax[tBool], y, 
         xlab="SELmax PAS", ylab=varname, main="SELmax PAS - Angle", col="darkgrey", ylim=mylims)
    for(j in seq(0,180,10)) {
      preddata$Angle <- j
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(SELmax_base + preddata$SELmax, y, col="blue")
    }
    lines(SELmax_base + preddata$SELmax, y, col="cyan")
    
    ##########
    
    preddata <- data.frame(SELmax=0:107)
    preddata$SELmax_CAS <- preddata$SELmax
    preddata$SELmax_PAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$Angle <- 30
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    preddata$pitch_first2 <- 45
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(SELmax_base + atab$SELmax[tBool], y, 
         xlab="SELmax CAS (dB)", ylab=varname, main="SELmax CAS - SeaState", col="darkgrey", ylim=mylims)
    for(j in seq(0,5,0.5)) {
      preddata$SeaState <- j
      preddata$SeaState_TH <- preddata$SeaState
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(SELmax_base + preddata$SELmax, y, col="blue")
    }
    lines(SELmax_base + preddata$SELmax, y, col="cyan")
    
    preddata <- data.frame(SELmax=0:107)
    preddata$SELmax_PAS <- preddata$SELmax
    preddata$SELmax_CAS <- 0
    preddata$SELmax2 <- preddata$SELmax
    preddata$SELmax_CAS2 <- preddata$SELmax_CAS
    preddata$SELmax_PAS2 <- preddata$SELmax_PAS
    preddata$Angle <- 30
    preddata$NS <- 0
    preddata$state <- as.character(3)
    preddata$depth_first <- 100
    preddata$Swell <- 0
    preddata$pitch_first2 <- 45
    y <- atab[tBool,varname]
    if(log_var) {y <- exp(y)}
    plot(SELmax_base + atab$SELmax[tBool], y, 
         xlab="SELmax PAS (dB)", ylab=varname, main="SELmax PAS - SeaState", col="darkgrey", ylim=mylims)
    for(j in seq(0,5,0.5)) {
      preddata$SeaState <- j
      preddata$SeaState_TH <- preddata$SeaState
      y <- predict(fit1$gam, newdata=preddata, type="response")
      if(log_var) {y <- exp(y)}
      lines(SELmax_base + preddata$SELmax, y, col="blue")
    }
    lines(SELmax_base + preddata$SELmax, y, col="cyan")
    
    
  }
  
  
###########################################################################################
##### RESPONE VARIABLES: Duration, Pause, ODBA_rms, Fluke rate, First ICI and Click rate  
  
###########################################################################################
##### NON-EXPOSURE MODELS (for model diagnostics)
  
  bfits <- list()
  
  ##################################################################################
  ##### Buzzes
  
  bfits[[1]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                         data=atab[bBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  bfits[[2]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                         data=atab[bBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  bfits[[3]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                         data=atab[bBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains
  
  bfits[[4]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first, 
                         data=atab[bBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  bfits[[5]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first, 
                         data=atab[bBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  bfits[[6]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first, 
                         data=atab[bBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  bfits[[7]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                         data=atab[bBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  
  # Estimated change in fluke rates with change in sea state
  diff(predict(bfits[[4]]$gam, newdata=data.frame(SeaState=c(1,2),pitch_first2=0.5, depth_first=200), type="response")*60)
  
  ##################################################################################
  ##### Refit without serial correlation
  
  bfits2 <- list()
  
  ##################################################################################
  ##### Buzzes
  
  bfits2[[1]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                          data=atab[bBool1,], family=gaussian(link = "identity"), random=list(ind=~1)),silent=T)
  
  bfits2[[2]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                          data=atab[bBool1,], family=binomial(link = "logit"), random=list(ind=~1)),silent=T)
  
  bfits2[[3]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                          data=atab[bBool1,], family=gaussian(link = "identity"), random=list(ind=~1)),silent=T)
  
  ##################################################################################
  ##### Regular click trains
  
  bfits2[[4]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first, 
                          data=atab[bBool2,], family=quasipoisson(link = "log"), random=list(ind=~1)),silent=T)
  
  bfits2[[5]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first, 
                          data=atab[bBool2,], family=Gamma(link = "log"), random=list(ind=~1)),silent=T)
  
  bfits2[[6]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first, 
                          data=atab[bBool2,], family=gaussian(link = "identity"), random=list(ind=~1)),silent=T)
  
  bfits2[[7]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first,
                          data=atab[bBool2,], family=binomial(link = "logit"), random=list(ind=~1)),silent=T)
  
  
  for(r in 1:7) {
    print(r)
    resp_var <- as.character(summary(bfits[[r]]$gam)$formula[2])
    if(class(bfits[[r]])[1]!="gamm") {print(bfits[[r]])} else {
      print(summary(bfits[[r]]$gam))
      
      
      tiff(filename=paste("3_model_diagnostics_", resp_var,".tiff",sep=""), 
           antialias="none", compression="zip",
           res = 600,
           width = 7*800, height = 6*800)
      
      par(mfrow=c(3,3), mar=c(4,4,4,3))
      #plot(1,1,bty="n", yaxt="n", xaxt="n",xlab="", ylab="", col=NA)
      #text(1,1,paste("Model ", r, ": ", resp_var, sep=""), cex=2)
      
      acf(residuals(bfits2[[r]]$lme,type="normalized"),main="Withour AR1")
      acf(residuals(bfits[[r]]$lme,type="normalized"),main="With AR1")
      
      gam.check(bfits[[r]]$gam, col=greyt)
      plot(bfits[[r]]$gam, shade=T,  scheme=3)
      vis.gam(bfits[[r]]$gam, view=c("pitch_first2", "SeaState"), 
              type="response", color="gray", plot.type="contour", main="")
      vis.gam(bfits[[r]]$gam, view=c("depth_first", "SeaState"), 
              type="response", color="gray", plot.type="contour", main="")
      #scan()
      options(graphics.record=FALSE) 
      dev.off()
    }
  }
  
  
##################################################################################
##### EXPOSURE MODELS (for statistical testing)
  
  efits <- list()
  
  ##################################################################################
  ##### Buzzes - Model structure 1 (SELsp + SELsp:sonar_angle)
  
  efits[[1]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax + SELmax:Angle,
                         data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[2]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax + SELmax:Angle,
                         data=atab[tBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[3]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax + SELmax:Angle,
                         data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains - Model structure 1 (SELsp + SELsp:sonar_angle)
  
  efits[[4]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax + SELmax:Angle,
                         data=atab[tBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[5]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax + SELmax:Angle,
                         data=atab[tBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[6]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax + SELmax:Angle,
                         data=atab[tBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Buzzes - Model structure 2 (SELsp, facing)
  
  efits[[7]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax2,
                         data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  
  efits[[8]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax2,
                         data=atab[tBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[9]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                         + NS + SELmax2,
                         data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains - Model structure 2 (SELsp, facing)
  
  efits[[10]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax2,
                          data=atab[tBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[11]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax2,
                          data=atab[tBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[12]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax2,
                          data=atab[tBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Buzzes - Model structure 3 (SELsp_CAS + SELsp_PAS + SELsp_CAS:sonar_angle + SELsp_PAS:sonar_angle)
  
  efits[[13]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[14]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[15]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains - Model structure 3 (SELsp_CAS + SELsp_PAS + SELsp_CAS:sonar_angle + SELsp_PAS:sonar_angle)
  
  efits[[16]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[17]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[18]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Buzzes - Model structure 4 (SELsp_CASfacing + SELsp_PASfacing)
  
  efits[[19]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  
  efits[[20]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[21]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains - Model structure 4 (SELsp_CASfacing + SELsp_PASfacing)
  
  efits[[22]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[23]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[24]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Buzzes - Model structure 5 (ti(SELsp) + ti(SELsp,sonar_angle))
  
  efits[[25]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[26]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[27]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains - Model structure 5 (ti(SELsp) + ti(SELsp,sonar_angle))
  
  efits[[28]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[29]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[30]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Buzzes - Model structure 6 (ti(SELsp_CAS) + ti(SELsp_CAS,sonar_angle) + ti(SELsp_PAS) + ti(SELsp_PAS,sonar_angle))
  
  efits[[31]] <- try(gamm(duration_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[32]] <- try(gamm(pause ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool1,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[33]] <- try(gamm(odba_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool1,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  ##################################################################################
  ##### Regular click trains - Model structure 6 (ti(SELsp_CAS) + ti(SELsp_CAS,sonar_angle) + ti(SELsp_PAS) + ti(SELsp_PAS,sonar_angle))
  
  efits[[34]] <- try(gamm(fluker ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=quasipoisson(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[35]] <- try(gamm(ICI_first ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=Gamma(link = "log"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[36]] <- try(gamm(CR_log ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=gaussian(link = "identity"), random=list(ind=~1), correlation=corAR1()),silent=T)


  ##################################################################################
  ##### P(Buzz) following regular click trains - fit to each 6 model structures
  
  efits[[37]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax + SELmax:Angle,
                          data=atab[tBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[38]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax2,
                          data=atab[tBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[39]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS + SELmax_CAS:Angle + SELmax_PAS + SELmax_PAS:Angle,
                          data=atab[tBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[40]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + SELmax_CAS2 + SELmax_PAS2,
                          data=atab[tBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[41]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax, bs="ts", k=3) + ti(SELmax, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  #efits[[42]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
  #                        + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
  #                        data=atab[tBool2,], family=binomial(link = "logit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  # message = iteration limit reached without convergence (10)
  
  #efits[[42]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
  #                        + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
  #                        data=atab[tBool2,], family=binomial(link = "probit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  efits[[42]] <- try(gamm(buzz ~ s(depth_first, bs="ts",k=5) + SeaState + SeaState:pitch_first2 + SeaState:depth_first
                          + NS + ti(SELmax_CAS, bs="ts", k=3) + ti(SELmax_CAS, Angle, bs="ts", k=3) + ti(SELmax_PAS, bs="ts", k=3) + ti(SELmax_PAS, Angle, bs="ts", k=3),
                          data=atab[tBool2,], family=binomial(link = "cauchit"), random=list(ind=~1), correlation=corAR1()),silent=T)
  
  
  # save(efits, file="2_model_fits.Rd")
  
  
  # Diagnostics & result checks
  for(r in 1:length(efits)) {
    print(r)
    resp_var <- as.character(summary(efits[[r]]$gam)$formula[2])
    if(class(efits[[r]])[1]!="gamm") {print(efits[[r]])} else {
      print(summary(efits[[r]]$gam))
      
      par(mfrow=c(3,4))
      plot(1,1,bty="n", yaxt="n", xaxt="n",xlab="", ylab="", col=NA)
      text(1,1,paste("Model ", r, ": ", resp_var, sep=""), cex=2)
      
      acf(residuals(efits[[r]]$lme,type="normalized"),main="standardized residual ACF")
      gam.check(efits[[r]]$gam, col=greyt)
      plot(efits[[r]]$gam, shade=T,  scheme=3)
      vis.gam(efits[[r]]$gam, view=c("pitch_first2", "SeaState"), 
              type="response", color="gray", plot.type="contour", main="")
      vis.gam(efits[[r]]$gam, view=c("depth_first", "SeaState"), 
              type="response", color="gray", plot.type="contour", main="")
      scan()
      
      print(summary(efits[[r]]$gam))
      
      tObj <- summary(efits[[r]]$gam)
      varname <- as.character(tObj$formula[2])
      log_var <- varname=="duration_log" | varname=="odba_log" | varname=="CR_log"
      if(!is.na(match(varname, c("duration_log","pause", "odba_log")))) {tBool <- tBool1 
      } else {tBool <- tBool2}
      
      pred_plots(efits[[r]], tBool, varname, log_var=log_var) # baseline plots
      scan()
      
      pred_plots2(efits[[r]], tBool, varname, log_var=log_var) # exposure plots
      scan()
    }
  }
  

  # Save results to a table
  rn <- 0
  for(r in 1:length(efits)) {
    if(class(efits[[r]])[1]=="gamm") {
      rn <- rn + 1
      tObj <- summary(efits[[r]]$gam)
      temp <- data.frame(rbind(tObj$p.table, tObj$s.table))
      temp <- cbind(ModelNo=rep(r, length(temp[,1])), Response=as.character(tObj$formula[2]), Covariates=row.names(temp), temp)
      if(rn==1) {tab <- temp} else {
        tab <- rbind(tab,temp)}
    }
  }
  View(tab)
  
  write.csv(tab, file="2_efits.csv") 
  
  
##################################################################################
##### Models for Sonar angle (Angle)
  
  atab$Angle_log <- logit(atab$Angle/180)
  
  sfits <- list()
  sfits[[1]] <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax, bs="ts",k=5),
                     data=atab[eBool2,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  
  sfits[[2]] <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax_PAS, bs="ts",k=5) + s(SELmax_CAS, bs="ts",k=5),
                     data=atab[eBool2,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  
  # refit without serial correlation
  sfits2 <- list()
  sfits2[[1]] <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax, bs="ts",k=5),
                      data=atab[eBool2,], family=gaussian, random=list(ind=~1))
  
  sfits2[[2]] <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax_PAS, bs="ts",k=5) + s(SELmax_CAS, bs="ts",k=5),
                      data=atab[eBool2,], family=gaussian, random=list(ind=~1))
  
  plot(sfits[[1]]$gam)
  acf(residuals(sfits[[1]]$lme,type="normalized"))
  summary(sfits[[1]]$gam)
  #  s(ExposureTime) 1.0447      4 4.399 6.57e-07 ***
  #  s(SELmax)       0.2494      4 0.081     0.22    
  
  
  # check if the relationships improve excluding poor mag data
  
  tBool <- atab$GMTtime > as.POSIXct("2017-01-01 00:00:01 GMT")
  temp1 <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax, bs="ts",k=5),
                data=atab[eBool2 & tBool,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  plot(temp1$gam)
  summary(temp1$gam)
  # s(ExposureTime) 9.213e-01      4 1.959 0.00324 **
  # s(SELmax)       3.522e-07      4 0.000 0.77617
  
  tBool <- atab$GMTtime > as.POSIXct("2017-01-01 00:00:01 GMT")
  temp2 <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax_PAS, bs="ts",k=5) + s(SELmax_CAS, bs="ts",k=5),
                data=atab[eBool2 & tBool,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  plot(temp2$gam)
  summary(temp2$gam)
  # s(ExposureTime) 9.067e-01      4 1.755 0.00431 **
  # s(SELmax_PAS)   1.815e-01      4 0.055 0.26442   
  # s(SELmax_CAS)   1.569e-06      4 0.000 0.33121 
  
  
  tiff(filename=paste("temp.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 7*800, height = 6*800)
  
  par(mfrow=c(3,3), mar=c(4,4,4,3))
  #plot(1,1,bty="n", yaxt="n", xaxt="n",xlab="", ylab="", col=NA)
  #text(1,1,paste("Model ", r, ": ", resp_var, sep=""), cex=2)
  
  acf(residuals(sfits2[[1]]$lme,type="normalized"),main="Without AR1")
  acf(residuals(sfits[[1]]$lme,type="normalized"),main="With AR1")
  
  gam.check(sfits[[1]]$gam, col=greyt)
  plot(sfits[[1]]$gam, shade=T,  scheme=3)
  
  vis.gam(sfits[[1]]$gam, view=c("ExposureTime", "SELmax"), 
          type="response", color="gray", plot.type="contour", main="")
  
  
  options(graphics.record=FALSE) 
  dev.off()
  
  
  ############### Other possible model structures
  
  sfits[[2]] <- gamm(Angle_log ~ ExposureTime + SELmax + ExposureTime:SELmax,
                     data=atab[eBool2,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  
  acf(residuals(sfits[[2]]$lme,type="normalized"))
  summary(sfits[[2]]$gam)
  # (Intercept)         -0.3226295  0.3340733  -0.966    0.335
  # ExposureTime        -0.0043959  0.0220589  -0.199    0.842
  # SELmax               0.0003580  0.0055379   0.065    0.948
  # ExposureTime:SELmax  0.0003213  0.0002633   1.220    0.223
  
  sfits[[3]] <- gamm(Angle_log ~ ti(ExposureTime, bs="ts",k=5) + ti(SELmax, bs="ts",k=5) + ti(ExposureTime, SELmax, bs="ts",k=5),
                     data=atab[eBool2,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  
  acf(residuals(sfits[[3]]$lme,type="normalized"))
  summary(sfits[[3]]$gam)
  # ti(ExposureTime)        1.032e+00      4 5.813 5.77e-07 ***
  # ti(SELmax)              8.819e-07      4 0.000    0.478    
  # ti(ExposureTime,SELmax) 2.578e+00     16 0.289    0.103  
  
  sfits[[4]] <- gamm(Angle_log ~ s(ExposureTime, bs="ts",k=5) + s(SELmax_CAS, bs="ts",k=5) + s(SELmax_PAS, bs="ts",k=5),
                     data=atab[eBool2,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  
  plot(sfits[[4]]$gam)
  acf(residuals(sfits[[4]]$lme,type="normalized"))
  summary(sfits[[4]]$gam)
  # s(ExposureTime) 1.092e+00      4 8.064  <2e-16 ***
  # s(SELmax_CAS)   1.922e-01      4 0.059   0.266    
  # s(SELmax_PAS)   8.372e-07      4 0.000   0.599
  

##################################################################################
##### Models for Surface angle (VAngle)   
  
  atab$VAngle_log <- logit(atab$pitch_first2/180)
  
  vfits <- list()
  vfits[[1]] <- gamm(VAngle_log ~ state + s(depth_first, bs="ts",k=5) + SeaState + SeaState:depth_first,
                     data=atab[bBool2,], family=gaussian, random=list(ind=~1), correlation=corAR1())
  
  # refit without serial correlation
  vfits2 <- list()
  vfits2[[1]] <- gamm(VAngle_log ~ s(depth_first, bs="ts",k=5) + + SeaState + SeaState:depth_first,
                      data=atab[bBool2,], family=gaussian, random=list(ind=~1))
  
  vfits2[[1]] <- gamm(VAngle_log ~ state + s(depth_first, bs="ts",k=5) + + SeaState + SeaState:depth_first,
                      data=atab[bBool2,], family=gaussian, random=list(ind=~1))
  
  plot(vfits[[1]]$gam)
  acf(residuals(vfits[[1]]$lme,type="normalized"))
  summary(vfits[[1]]$gam)
  # SeaState              6.953e-02  5.510e-02   1.262    0.207
  # SeaState:depth_first -1.644e-05  6.747e-05  -0.244    0.807
  
  
  tiff(filename=paste("temp.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 5*800, height = 6*800)
  
  par(mfrow=c(3,2), mar=c(4,4,4,3))
  #plot(1,1,bty="n", yaxt="n", xaxt="n",xlab="", ylab="", col=NA)
  #text(1,1,paste("Model ", r, ": ", resp_var, sep=""), cex=2)
  
  acf(residuals(vfits2[[1]]$lme,type="normalized"),main="Without AR1")
  acf(residuals(vfits[[1]]$lme,type="normalized"),main="With AR1")
  
  gam.check(vfits[[1]]$gam, col=greyt)
  #plot(vfits[[1]]$gam, shade=T,  scheme=3)
  
  #vis.gam(vfits[[1]]$gam, view=c("depth_first", "SeaState"), 
  #        type="response", color="gray", plot.type="contour", main="")
  
  
  options(graphics.record=FALSE) 
  dev.off()
  