  Sys.setenv(TZ='GMT')

##################################################################################
#####  R tools & colours

  ilogit <- function(x) {exp(x)/(exp(x)+1)}
  
  library(nimble)
  
  #library(R2jags)
  set.seed(0)
  library(lattice)
  library(coda)
  
  R2 <- function(obs, fits) {
    SSE <- sum((obs-fits)^2)
    SST <- sum((obs-mean(obs))^2)
    return(1-SSE/SST)
  }
  
  myPlots <- function(fit, par.name="tau[1]", titles="", diag=T) {

      tau.mcmc <- mcmc.list(mcmc(fit$samples$chain1[,par.name]),
                            mcmc(fit$samples$chain2[,par.name]),
                            mcmc(fit$samples$chain3[,par.name]))
      
      temp <- gelman.diag(tau.mcmc)
      par(mfrow=c(1,2))
      traceplot(tau.mcmc)
      title(par.name)
      gelman.plot(tau.mcmc, main=round(temp$psrf[[1]],3), auto.layout=F)
      
      par(mfrow=c(1,1))
      lattice.options()
      plotObj <- densityplot(tau.mcmc, auto.layout=F,
                             main=list(label=paste(par.name, "mean",
                                                   signif(summary(tau.mcmc)$statistics["Mean"]))))
      
      print(plotObj)
      
    
    return(tau.mcmc)
  }
  
  
  myalpha <- 80
  greyt <- rgb(0, 0, 0, max = 255, alpha = myalpha)
  redt <- bluet <- rgb(255, 0, 0, max = 255, alpha = myalpha)
  bluet <- rgb(0, 0, 255, max = 255, alpha = myalpha)
  bluet2 <- rgb(100, 0, 255, max = 255, alpha = myalpha)

##################################################################################
##### Data
  
  load("Supplement_data.Rd")
  
##################################################################################
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
  
##################################################################################
##### Filters
  
  tBool0 <- (atab$state=="2" | atab$state=="3" | atab$state=="4") & !is.na(atab$state) 
  
  # Buzzes
  tBool1 <- tBool0 &  atab$label=="BUZZ"
  tBool1[is.na(tBool1)] <- FALSE
  sum(tBool1) # 1910
  bBool1 <- tBool1 & atab$ExposureType=="BP"
  
  # Regular click trains
  tBool2 <- tBool0 & atab$label=="RC" & atab$ClickN > 5 & atab$CR <= 5
  tBool2[is.na(tBool2)] <- FALSE
  sum(tBool2) # 5963
  bBool2 <- tBool2 & atab$ExposureType=="BP"
  
##################################################################################
##### Set up data for the model

  tBool <- !is.na(atab$depth_first) & !is.na(atab$pitch_first) & !is.na(atab$zpeak0) & !is.na(atab$SeaState) 
  tBool <- tBool & !(is.na(atab$Angle) & (atab$SELmax_PAS+atab$SELmax_CAS)>0)

  tBool <- tBool2 & tBool
  
  N <- sum(tBool)
  NW <- max(atab$ID[tBool])
  
  N # 4045 (previous sample size, including NA for sea state: 5963)
  NW # 15
  
  # Constants
   sp.constants <- list(N=N, NW=NW, ID=atab$ID[tBool])

  # Data
    # Response variables
      sp.data <- list(zpeak0 = atab$zpeak0[tBool])
      sp.data$zpeak0[sp.data$zpeak0 >= 187.9] <- NA
      naBool <- is.na(sp.data$zpeak0)
      sp.data$isCensored <- as.numeric(naBool) # Right-censored: 1
      sp.data$censorLimitVec <- sp.data$zpeak0 - 0.1 # rep(floor(min(sp.data$zpeak0, na.rm=T)-1),N)
      sp.data$censorLimitVec[naBool] <- 188
    
      #summary(sp.data$zpeak0)
      #table(sp.data$isCensored, sp.data$censorLimitVec)
      
    # Covariates
      # Baseline
      sp.data$SeaState <- log10(atab$SeaState[tBool])
      sp.data$depth_first <- atab$depth_first[tBool]/100 # how many 100s meters away from surface?
      sp.data$pitch_first <- (pi/2 - atab$pitch_first[tBool])/pi*18 # how many tens of degrees away from surface?
      # Exposures
      sp.data$SELmax <- atab$SELmax[tBool]/10
      sp.data$Angle <- atab$Angle[tBool]/10
      sp.data$Angle[is.na(sp.data$Angle)] <- 18
      
      table(is.na(sp.data$zpeak0))
      length(sp.data$zpeak0) 
      sum(!is.na(sp.data$zpeak0))/length(sp.data$zpeak0) # 0.7364648
      sum(is.na(sp.data$zpeak0))/length(sp.data$zpeak0) # 0.2635352 
      tapply(is.na(sp.data$zpeak0), sp.constants$ID, mean)
      x <- tapply(is.na(sp.data$zpeak0), sp.constants$ID, mean)
      sum(x==0) # 5
      sum(x> 0 & x < 0.5) # 6
      sum(x >= 0.5) # 3
     
      
##################################################################################
##### MODEL SETUP
      
  # Initial value ranges
  
    b0_init <- c(min(sp.data$zpeak0, na.rm=T), max(sp.data$zpeak0, na.rm=T))
    tau_init <- c(0.01, 0.2)
    
  sp.inits <- function() { 
    
    zpeak0_init <- rep(NA, N)
    naBool <- is.na(sp.data$zpeak0)
    zpeak0_init[naBool] <- runif(sum(naBool), 188.1, 200) #min(sp.data$zpeak0, na.rm=T),  max(sp.data$zpeak0, na.rm=T))
    
    return(list(
      zpeak0=zpeak0_init,
      phi=runif(1, 0, 1),
      b0=runif(1, b0_init[1],  b0_init[2]),
      b0_tau=runif(1, tau_init[1],  tau_init[2]),
      b0_ID=runif(NW, b0_init[1], b0_init[2]),
      mp1=runif(1, 0, 5),  # Sea state
      mr1=runif(1, 0, 5),  # masking release f(depth)
      mr2=runif(1, 0, 5),  # masking release f(pitch)
      sp1=runif(1, 0, 5),  # CAS/PAS  
      sr1=runif(1, 0, 5),  # masking release f(Angle)
      tau=runif(NW, tau_init[1],  tau_init[2])
      ))
    }
  inits_list <- sp.inits()
  
  # Parameters
    sp.params <- c("tau", "b0_tau", "b0", "mp1", "mr1", "mr2", "sp1", "sr1", "y", "phi") #, "phi", "epsilon")
  
  # Model code
    source("3_Lombard_nimble_model.R")
    source("3_Lombard_nimble_model_null1.R")
    source("3_Lombard_nimble_model_null2.R")
    
  ## Final check
  print(str(sp.data))
  print(sp.params)
  sp.model
 
##################################################################################
##### MODEL FIT
  
  # Build & check model
  model <- nimbleModel(code=sp.model, constants=sp.constants, data=sp.data, 
                       inits=inits_list)
  
  model$initializeInfo() 
  # Missing values (NAs) or non-finite values were found
  # model$getNodeNames()

  # Check model-generated values for zpeak0
  set.seed(1)
  model$simulate("zpeak0")
  hist(model$zpeak0, col="lightgrey")
  hist(sp.data$zpeak0, add=T, col="darkgrey")
  
  # Fit model 
  fit <- nimbleMCMC(sp.model, constants=sp.constants, data=sp.data, 
                    inits=sp.inits, monitors=sp.params, 
                    summary=TRUE, samplesAsCodaMCMC = TRUE,
                    nchains=3, niter=24000, nburnin=12000, 
                    thin=max(1, floor(3 * (24000-12000) /1000)))
  
  # Null model for sea state
  fit_null1 <- nimbleMCMC(sp.model_null1, constants=sp.constants, data=sp.data, 
                    inits=sp.inits, monitors=sp.params, 
                    summary=TRUE, samplesAsCodaMCMC = TRUE,
                    nchains=3, niter=24000, nburnin=12000, 
                    thin=max(1, floor(3 * (24000-12000) /1000))) 
  
  # Null model for sonar
  fit_null2 <- nimbleMCMC(sp.model_null2, constants=sp.constants, data=sp.data, 
                    inits=sp.inits, monitors=sp.params, 
                    summary=TRUE, samplesAsCodaMCMC = TRUE,
                    nchains=3, niter=24000, nburnin=12000, 
                    thin=max(1, floor(3 * (24000-12000) /1000))) 
  
  save(sp.constants, sp.data, sp.inits, sp.params, sp.model, fit, fit_null1, fit_null2, 
       file="3_Lombard_nimble_model_objects.Rd")
  
  
##################################################################################
##### CONVERGENCE CHECKS
  
  temp <- myPlots(fit, "b0[1]"); plot(temp) # Intercept
  
  temp <- myPlots(fit, "mp1[1]"); plot(temp) # Sea state
  temp <- myPlots(fit_null1, "mp1[1]"); plot(temp) # Sea state
  temp <- myPlots(fit_null2, "mp1[1]"); plot(temp) # Sea state
  
  temp <- myPlots(fit, "mr1[1]"); plot(temp) # Masking release f(depth)
  temp <- myPlots(fit, "mr2[1]"); plot(temp) # Masking release f(surface angle)
  
  temp <- myPlots(fit, "sp1[1]"); plot(temp) # Sonar masking potential
  temp <- myPlots(fit_null1, "sp1[1]"); plot(temp) # Sea state
  temp <- myPlots(fit_null2, "sp1[1]"); plot(temp) # Sea state
  
  temp <- myPlots(fit, "sr1[1]"); plot(temp) # Sonar masking release f(Angle)
  temp <- myPlots(fit, "tau[1]"); plot(temp)
  temp <- myPlots(fit, "tau[2]"); plot(temp)
  temp <- myPlots(fit, "tau[15]"); plot(temp)
  temp <- myPlots(fit, "phi[1]"); plot(temp) # 
  
  
##################################################################################
##### GOODNESS-OF-FIT
  
  ##### Plot model fit
  
    par(mfrow=c(1,1), mar=c(4,4,3,3))
    
    parnames <- names(fit$summary$all.chains[,"Median"])
    tempI <- startsWith(parnames, "y")
    
    tcol <- rep(greyt, sp.constants$N)
    tcol[sp.data$SELmax>0] <- redt
    
    y <- fit$summary$all.chains[tempI,"Median"]
    plot(sp.data$zpeak0, y, xlab="Observed 0-to-peak", ylab="Predicted 0-to-peak", 
         col=tcol, pch=16, cex=0.5+sp.data$SELmax/15)
    abline(v=188, col="blue")
    abline(h=188, col="blue")
    abline(0,1, col="red")
    
  ##### Compare R-squared
    
    parnames <- names(fit$summary$all.chains[,"Median"])
    tempI <- startsWith(parnames, "y")
    y <- fit$summary$all.chains[tempI,"Median"]
    y_null1 <- fit_null1$summary$all.chains[tempI,"Median"]
    y_null2 <- fit_null2$summary$all.chains[tempI,"Median"]
    R2(sp.data$zpeak0[!is.na(sp.data$zpeak0)], y[!is.na(sp.data$zpeak0)]) # 0.08595962
    R2(sp.data$zpeak0[!is.na(sp.data$zpeak0)], y_null1[!is.na(sp.data$zpeak0)]) # 0.06048347
    R2(sp.data$zpeak0[!is.na(sp.data$zpeak0)], y_null2[!is.na(sp.data$zpeak0)]) # 0.088164
  
    
##################################################################################
##### MODEL ESTIMATES

  ##### Summary of posterior mean estimates
    
    parnames <- c("b0[1]", "b0_tau[1]", "phi[1]", "mp1[1]", "mr1[1]", "mr2[1]", "sp1[1]", "sr1[1]")
    stab <- fit$summary$all.chains[c(1:23),]
    stab1 <- fit_null1$summary$all.chains[c(1:23),]
    stab2 <- fit_null2$summary$all.chains[c(1:23),]
    
    stab <- cbind(stab, stab1, stab2)
    
    #write.csv(stab, file="3_Lombard_nimble_model_output.csv")
  
  ##### Figure 7
    
  get_preds <- function(preddata, intercept=FALSE) {
    
    if(intercept) {b0 <- fit$summary$all.chains["b0[1]","Mean"]} else {b0 <- 0}
    mp1 <- fit$summary$all.chains["mp1[1]","Median"]
    mr1 <- fit$summary$all.chains["mr1[1]","Median"]
    mr2 <- fit$summary$all.chains["mr2[1]","Median"]
    sp1 <- fit$summary$all.chains["sp1[1]","Median"]
    sr1 <- fit$summary$all.chains["sr1[1]","Median"]
    
    mp <- mp1[1]*log10(preddata$SeaState)
    mr <-  ilogit(-10 + mr1[1]*preddata$depth_first + mr2[1]*preddata$pitch_first)
    mt <- mp*(1-mr)
    
    sp <- (sp1[1]*preddata$SELmax)
    sr <-  ilogit(-10 + sr1[1]*preddata$Angle)
    st <- sp*(1-sr)
    
    y <- b0 + mt + st
    
    return(y)
  }
  get_CI <- function(preddata, intercept=FALSE) {
    
    A <- rbind(fit$samples$chain1, fit$samples$chain2, fit$samples$chain3)
    
    if(intercept) {b0 <- A[,"b0[1]"]} else {b0 <- 0}

    mp1 <- A[,"mp1[1]"]
    mr1 <- A[,"mr1[1]"]
    mr2 <- A[,"mr2[1]"]
    sp1 <- A[,"sp1[1]"]
    sr1 <- A[,"sr1[1]"]

    P <- matrix(NA, dim(preddata)[1], length(mp1))
    
    preddata$CI_lo <- NA
    preddata$CI_hi <- NA
      
    for(j in 1:(dim(P)[1])) {
    
      mp <- mp1*log10(preddata$SeaState[j])
      mr <-  ilogit(-10 + mr1*preddata$depth_first[j] + mr2*preddata$pitch_first[j])
      mt <- mp*(1-mr)
      
      sp <- (sp1*preddata$SELmax[j])
      sr <-  ilogit(-10 + sr1*preddata$Angle[j])
      st <- sp*(1-sr)
    
      y <- b0 + mt + st
      
      preddata$CI_lo[j] <- quantile(y, 0.025)
      preddata$CI_hi[j] <- quantile(y, 0.975)
      
    }
    
    return(list(CI_low=preddata$CI_lo, CI_hi=preddata$CI_hi))
  }
  

  tiff(filename=paste("Figure_7.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 6*800, height = 5.5*800)
  
  mymax <- 8
  par(mfrow=c(2,2), mar=c(4,4,3,2))
  
  for(d in c(0,90)) {
    
    preddata <- data.frame(depth_first=seq(0,12,0.2))
    preddata$Swell <- 0
    preddata$SeaState <- 5
    preddata$pitch_first <-  d/10
    preddata$SELmax <- 0
    preddata$Angle <- 0
    preddata$Buzz <- 0
    y <- get_preds(preddata)
    y_max <- mymax
    
    plot(preddata$depth_first*100, y, type="l", col=NA,
         xlab=c("Depth (m)"), ylab="dB change in 0-to-peak level", ylim=c(0, y_max),
         main=paste("Vertical angle ", d, " deg", sep=""))
    
    for(j in 1:5) {
      preddata$SeaState <- j
      y <- get_preds(preddata)
      lines(preddata$depth_first*100, y)
      text(10, y[1], j)
      
      tObj <- get_CI(preddata)
      y_lo <- tObj$CI_lo
      y_hi <- tObj$CI_hi
      lines(preddata$depth_first*100, y_lo, col="lightblue", lty=2)
      lines(preddata$depth_first*100, y_hi, col="salmon", lty=2)
    }
  }
  # 

    preddata <- data.frame(pitch_first = seq(0,180,1)/10)
    preddata$depth_first <- 1
    preddata$Swell <- 0
    preddata$SeaState <- 2
    preddata$Angle <-  0
    preddata$SELmax <- 0
    preddata$Buzz <- 0
    y <- get_preds(preddata)
    y_max <- mymax
    
    plot(preddata$pitch_first*10, y, type="l", col=NA,
         xlab=c("Vertical angle (deg)"), ylab="dB change in 0-to-peak level", ylim=c(0, y_max),
         main="Dive depth 100 m")
    
    for(j in 1:5) {
      preddata$SeaState <-j
      y <- get_preds(preddata)
      lines(preddata$pitch_first*10, y)
      text(5, y[1], j, col="darkgrey")
      
      tObj <- get_CI(preddata)
      y_lo <- tObj$CI_lo
      y_hi <- tObj$CI_hi
      lines(preddata$pitch_first*10, y_lo, col="lightblue",  lty=2)
      lines(preddata$pitch_first*10, y_hi, col="salmon",  lty=2)
    }
    
  
  for(d in c(1)) {
    
    SEL_pred <- c(90, 130, 170)
    
    preddata <- data.frame(Angle = seq(0,80,1)/10)
    preddata$depth_first <- 1
    preddata$Swell <- 0
    preddata$SeaState <- d
    preddata$pitch_first <-  d/10
    preddata$SELmax <- (max(SEL_pred)-SELmax_base)/10
    preddata$Buzz <- 0
    y <- get_preds(preddata)
    y_max <- mymax
    
    plot(preddata$Angle*10, y, type="l", col="grey",
         xlab=c("Angle to source (deg)"), ylab="dB change in 0-to-peak level", ylim=c(0, y_max),
         main="Sonar effect at Sea state 1")
    
    for(j in 1:length(SEL_pred)) {
      preddata$SELmax <- (SEL_pred[j]-SELmax_base)/10
      y <- get_preds(preddata)
      lines(preddata$Angle*10, y)
      text(5, y[1], SEL_pred[j], col="darkgrey")
    
      tObj <- get_CI(preddata)
      y_lo <- tObj$CI_lo
      y_hi <- tObj$CI_hi
      lines(preddata$Angle*10, y_lo, col="lightblue",  lty=2)
      lines(preddata$Angle*10, y_hi, col="salmon",  lty=2)
      }
    
  }
  
  options(graphics.record=FALSE) 
  dev.off()
  
  max(y) # 4.390377
  y_lo[which.max(y)] # 0.5863372
  y_hi[which.max(y)] # 8.477364
  