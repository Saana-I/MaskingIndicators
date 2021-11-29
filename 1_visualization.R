  Sys.setenv(TZ='GMT')

##################################################################################
##### R tools & colours
  
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
#####  Data transforms

  atab$CR <- atab$ClickN/atab$duration
  atab$CR2 <- 1/atab$ICI_med
  
  SELmax_base <- min(atab$SELmax, na.rm=T)-1
  atab$SELmax[is.na(atab$SELmax)] <- SELmax_base
  atab$SELmax <- atab$SELmax - SELmax_base
  
  SELcum_base <- min(atab$SELcum, na.rm=T)-1
  atab$SELcum[is.na(atab$SELcum)] <- SELcum_base
  atab$SELcum <- atab$SELcum - SELcum_base
  
  atab$SELmax_CAS <- atab$SELmax
  atab$SELmax_CAS[atab$ExposureType=="MPS" | atab$ExposureType=="HPS"] <- 0
  
  atab$SELmax_PAS <- atab$SELmax
  atab$SELmax_PAS[atab$ExposureType=="CAS"] <- 0
  
  atab$SELmax_CAS_fac <- round(atab$SELmax_CAS/20)*20
  atab$SELmax_PAS_fac <- round(atab$SELmax_PAS/20)*20
  
  atab$pitch_first2 <- (pi/2 - atab$pitch_first)/pi*180 
  atab$pitch_mean2 <- (pi/2 - atab$pitch_mean)/pi*180 

  atab$fluker <- atab$fluken/atab$duration
  atab$fluker_log <- log(atab$fluker)
  
##################################################################################
##### Filters

  # Selecting regular click trains
  
  tBool0 <- (atab$state=="2" | atab$state=="3" | atab$state=="4") & !is.na(atab$state) 
  tBool0 <- tBool0 & atab$label=="RC" & atab$ClickN > 5 & atab$CR <= 5
  tBool0[is.na(tBool0)] <- FALSE
  
  # Seecting buzz trains
  
  tBool0B <- (atab$state=="2" | atab$state=="3" | atab$state=="4") & !is.na(atab$state)
  tBool0B <- tBool0B & atab$label=="BUZZ"
  tBool0B[is.na(tBool0B)] <- FALSE
  
  
##################################################################################
##### Figure 1
  
  whales <- unique(atab$ind)
  
  tcol <- rep(NA, length(atab$ind))
  tcol[atab$label=="RC"] <- greyt
  tcol[atab$label=="BUZZ"] <- bluet
  
  plot_temp <- function(Exposure="CAS") {
    plot(-40:40, -40:40, col="NA", xlab="Time (min)", ylab="Depth (m)", ylim=c(-2000,0), las=1)
    
    for(w in 1:length(whales)) {
      
      baseI <- which(etab$ind==whales[w] & etab$Session=="Baseline")
      expI <- which(etab$ind==whales[w] & etab$Session==Exposure)
      
      if(length(expI)>0) {
        expI <- expI[1]
        
        ST <- etab$Start_UTC[expI]-40*60
        ET <- etab$End_UTC[expI]
        
        tBool1 <- atab$ind==whales[w] & (atab$GMTtime >= ST  & atab$GMTtime < ET)
        tBool2 <- atab$ind==whales[w] & ((atab$GMTtime+atab$duration) >= ST  & (atab$GMTtime+atab$duration) < ET)
        tBool <- tBool1 | tBool2 
        
        tvec <- (as.numeric(atab$GMTtime[tBool])-as.numeric(etab$Start_UTC[expI]))/60
        tvec2 <- tvec+atab$duration[tBool]/60
        
        segments(x0=tvec, x1=tvec2,y0=-atab$depth_first[tBool], y1=-atab$depth_last[tBool], 
                 col=tcol[tBool], lwd=4*(1-atab$Angle[tBool]/180))
        
      }
    }
  }
  
  tiff(filename=paste("Figure_1.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 6.5*800, height = 4.5*800)
  
  par(mfrow=c(2,2), mar=c(4,4,2,1))
  
  plot_temp("NS")
  mtext(side=3, at=-30, "No-sonar")
  plot_temp("CAS")
  mtext(side=3, at=-30, "CAS")
  legend(0,-1000, title="Sonar angle (deg)", legend=c(30, 60, 90), lwd=4*(1-c(30, 60, 90)/180), col=greyt)
  plot_temp("MPAS")
  mtext(side=3, at=-30, "MPAS")
  plot_temp("HPAS")
  mtext(side=3, at=-30, "HPAS")
  
  options(graphics.record=FALSE) 
  dev.off()
  
  
##################################################################################
##### Figure 2
  
  tcol <- rep(greyt, length(tBool0))
  tcol[atab$state=="3"] <- bluet
  
  x <- atab$SELmax + SELmax_base
  tcex <- x/80
  
  
  tiff(filename=paste("Figure_2.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 6*800, height = 5.5*800)
  
  par(mfrow=c(2,2),mar=c(4,4,3,1))
  
  SS_min <- c(0, 2, 3, 4)
  SS_max <- c(2, 3, 4, 5)
  
  for(j in 1:length(SS_min)) {
    
    tBool <- tBool0 & atab$ExposureType=="BP" & atab$SeaState >= SS_min[j] & atab$SeaState < SS_max[j]
    if(j==length(SS_min)) {tBool <- tBool0 & atab$ExposureType=="BP" & atab$SeaState >= SS_min[j] & atab$SeaState <= SS_max[j]}
    
    plot(atab$Angle[tBool], atab$pitch_first2[tBool], cex=tcex[tBool], 
         pch=16, col=greyt, ylim=c(0,180), xlim=c(0,180), cex.axis=1.2,
         xlab="", ylab="")
    
    tBool2 <- tBool0 & atab$ExposureType!="BP" & atab$Focal==TRUE & atab$SeaState >= SS_min[j] & atab$SeaState < SS_max[j]
    
    points(atab$Angle[tBool2], atab$pitch_first2[tBool2], cex=tcex[tBool2], 
           col=atab$col[tBool2], pch=atab$pch[tBool2])
    
    mtext("Sonar angle (deg)", side=1, line=2.6)
    mtext("Surface angle (deg)", side=2, line=2.6)
    mtext(paste("Sea state [", SS_min[j], " - ", SS_max[j], ")", sep=""), side=3, line=0.5, font=2)
    if(j==length(SS_min)) {mtext(paste("Sea state [", SS_min[j], " - ", SS_max[j], "]", sep=""), side=3, line=0.5, font=2)}
    
    if(j==1) {
      tvec <- c(-pi/2, 0, pi/2)
      legend(-7, 185,  seq(90,170,40), title=expression(paste(plain("SEL (dB re 1"), mu,plain(Pa)^2, "s)", sep="")), pch=15, col=bluet, 
             pt.cex= seq(90,170,40)/80, cex=1, bty="n")
    }
  }
  
  options(graphics.record=FALSE) 
  dev.off()
  
  
  
##################################################################################
##### Figure 4 
  
  
  tiff(filename=paste("Figure_4.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 5.5*800, height = 6.5*800)
  
  tcol <- rep(greyt, length(tBool0))
  tcol[atab$state=="3"] <- bluet
  
  col_var <- (atab$pitch_first+pi/2)/(pi)

  x <- atab$SELmax + SELmax_base
  mycex <- (atab$SELmax + SELmax_base)/80
  
  par(mfrow=c(3,2),mar=c(4,4,3,1), bg="white")
  
  SS_min <- c(0,1, 2,   2.5,  3,4)
  SS_max <- c(1,2, 2.5, 3,    4,6)
  
  for(j in 1:length(SS_min)) {
    
    tBool <- tBool0 & atab$SeaState >= SS_min[j] & atab$SeaState < SS_max[j]
    tBool[is.na(tBool)] <- FALSE
    tBool[is.na(col_var)] <- FALSE
    
    eBool <- tBool & atab$ExposureType=="BP"
    plot(atab$depth_first[eBool], atab$ICI_first[eBool], cex=0.2+mycex[eBool],
         pch=1, col=greyt, ylim=c(0,2), xlim=c(0,2000), cex.axis=1.2,
         xlab="", ylab="")
    points(atab$depth_first[tBool], atab$ICI_first[tBool], col=gray(1-col_var[tBool]), pch=16, cex=mycex[tBool])
    eBool <- tBool & atab$ExposureType=="CAS"
    points(atab$depth_first[eBool], atab$ICI_first[eBool], col=bluet, pch=1, cex=0.1+mycex[eBool])
    eBool <- tBool & (atab$ExposureType=="HPS" | atab$ExposureType=="MPS")
    points(atab$depth_first[eBool], atab$ICI_first[eBool], col=redt, pch=1, cex=0.1+mycex[eBool])
    abline(h=c(0.02, 0.5, 1, 1.5), lty=2, col=greyt)
    if(j==1) {
      tvec <- c(-pi/2, 0, pi/2)
      legend(1300, 1.05, c("-45","0","45"), title="Pitch (deg)", pch=16, 
             col=gray(1-((c(-pi/4,0,pi/4)+pi/2)/(pi))), bg="white", cex=1.5, pt.cex=1)#, bty="n")
      legend(1300, 2, seq(80,180,30)[c(2,4)], title="SEL (dB)", pch=16, 
             col=greyt, cex=1.5, pt.cex=seq(80,180,30)[c(2,4)]/80, bg="white")#, bty="n")
      legend(740, 2, c("CAS", "PAS"), cex=1.5, col=c(bluet, redt), pch=1)
    }
    
    mtext("Depth (m)", side=1, line=2.6)
    mtext("First ICI (s)", side=2, line=2.6)
    mtext(paste("Sea state [", SS_min[j], ", ", SS_max[j], ")", sep=""), side=3, line=0.5, font=2)
    
  }
  
  options(graphics.record=FALSE) 
  dev.off()

  
##################################################################################
##### Figure 5 
  
  tiff(filename=paste("Figure_5.tiff",sep=""), 
       antialias="none", compression="zip",
       res = 600,
       width = 5.5*800, height = 6.5*800)
  
  tcol <- rep(greyt, length(tBool0))
  tcol[atab$state=="3"] <- bluet
  
  col_var <- (atab$pitch_mean +pi/2)/(pi)

  mycex <- (atab$SELmax + SELmax_base)/80
  
  par(mfrow=c(3,2),mar=c(4,4,3,1), bg="white")
  
  SS_min <- c(0, 1, 2,   2.5,  3, 4)
  SS_max <- c(1, 2, 2.5, 3,    4, 6)
  
  for(j in 1:length(SS_min)) {
    
    tBool <- tBool0 & atab$SeaState >= SS_min[j] & atab$SeaState < SS_max[j]
    tBool[is.na(tBool)] <- FALSE
    tBool[is.na(col_var)] <- FALSE
    
    eBool <- tBool & atab$ExposureType=="BP"
    plot(atab$depth_first[eBool], atab$CR2[eBool], cex=0.2+mycex[eBool],
         pch=1, col=greyt, ylim=c(0,5), xlim=c(0,2000), cex.axis=1.2,
         xlab="", ylab="")
    points(atab$depth_first[tBool], atab$CR2[tBool], col=gray(1-col_var[tBool]), pch=16, cex=mycex[tBool])
    eBool <- tBool & atab$ExposureType=="CAS"
    points(atab$depth_first[eBool], atab$CR2[eBool], col=bluet, pch=1, cex=0.1+mycex[eBool])
    eBool <- tBool & (atab$ExposureType=="HPS" | atab$ExposureType=="MPS")
    points(atab$depth_first[eBool], atab$CR2[eBool], col=redt, pch=1, cex=0.1+mycex[eBool])
    abline(h=0:5, lty=2, col=greyt)
    if(j==1) {
      tvec <- c(-pi/2, 0, pi/2)
      legend(1300, 2.6, c("-45","0","45"), title="Pitch (deg)", pch=16, 
             col=gray(1-((c(-pi/4,0,pi/4)+pi/2)/(pi))), bg="white", cex=1.5, pt.cex=1)#, bty="n")
      legend(1300, 5, seq(80,180,30)[c(2,4)], title="SEL (dB)", pch=16, 
             col=greyt, cex=1.5, pt.cex=seq(80,180,30)[c(2,4)]/80, bg="white")#, bty="n")
      legend(740, 5, c("CAS", "PAS"), cex=1.5, col=c(bluet, redt), pch=1)
    }
    mtext("Depth (m)", side=1, line=2.6)
    mtext("Click rate (/s)", side=2, line=2.6)
    mtext(paste("Sea state [", SS_min[j], ", ", SS_max[j], ")", sep=""), side=3, line=0.5, font=2)
    
  }
  
  options(graphics.record=FALSE) 
  dev.off()
  
  
  
##################################################################################
##### Figure 6  

  atab_orig <- atab
  atab$ID <- as.numeric(as.factor(atab$ind))
  atab <- atab[tBool0 & !is.na(atab$depth_first),]
  
  varname <- c("fluker", "zpeak0", "ICI_first", "CR", "CR2")
  ylabs <- c("Stroke rate (/s)", "Zero-to-peak level (dB)", "First ICI (s)", "Click rate (s-1)", "Click rate (s-1)")
  
  for(v in 1:2) { # 1:length(varname)) {
    
    atab$y <- atab[,varname[v]]
    atab2 <- atab[!is.na(atab$y),]# atab[!is.na(atab$y) & atab$depth_first < 250,]
    
    
    tiff(filename=paste("Figure_6_", varname[v],".tiff",sep=""), 
         antialias="none", compression="zip",
         res = 600,
         width = 7*800, height = 6*800)
    
    par(mfrow=c(3,3), mar=c(4,4,3,2))
    
    deg_add <- 45
    
    for(d in c(0,45,90)) {
      #d <- 0
      
      
      mylims <- quantile(atab2$y, c(0.05, 0.95))
      if(varname[v]=="ICI_first") {mylims <- c(0,1.5)}
      if(varname[v]=="CR") {mylims <- c(0.5,3.5)}
      if(varname[v]=="CR2") {mylims <- c(0.5,3.5)}
      
      ### Baseline
      sBool <- (atab2$pitch_first2 >= d & atab2$pitch_first2 < (d+ deg_add)) & atab2$depth_first < 200
      sBool <- sBool & atab2$SELmax==0
      sBool[is.na(sBool)] <- FALSE
      whales <- unique(atab2$ID)
      
      plot(atab2$SeaState[sBool], atab2$y[sBool], col=greyt, xlab="Sea state (Beaufort)", ylab=ylabs[v], pch=16, cex=0.5,
           main=paste("<200 m, surface angle: ", d, "-", d+deg_add," degrees", sep=""), 
           ylim=mylims, xlim=c(0,5), las=1)
      
      for(w in whales) {
        
        wBool <- sBool & atab2$ID==w
        if(sum(wBool, na.rm=T)>0) {
          y <- tapply(atab2$y[wBool], round(atab2$SeaState[wBool]), median, na.rm=T)
          lines(as.numeric(names(y)), y, lwd=1.2, type="b", pch=16)
        }
      }
      
      ### Sonar exposures - shallow
      sBool <- (atab2$Angle >= d & atab2$Angle < (d+ deg_add)) & atab2$depth_first < 200
      sBool <- sBool & atab2$SELmax>0
      
      plot(SELmax_base+atab2$SELmax_PAS[sBool], atab2$y[sBool], col=redt, 
           xlab=expression(paste(plain("SEL (dB re 1"), mu,plain(Pa)^2, "s)", sep="")),  ylab=ylabs[v], pch=16, cex=0.5,
           main=paste("<200 m, sonar angle: ", d, "-", d+deg_add," degrees", sep=""),
           ylim=mylims, xlim=c(80,170), las=1)
      points(SELmax_base+atab2$SELmax_CAS[sBool], atab2$y[sBool], col=bluet, pch=16, cex=0.5)
      
      # Add PAS
      for(w in whales) {
        
        wBool <- sBool & (atab2$SELmax_PAS > 0) & atab2$ID==w
        if(sum(wBool, na.rm=T)>0) {
          #lines(SELmax_base+atab2$SELmax_PAS[wBool]*10, atab2$y[wBool], lwd=1.2, type="b", pch=16, col="red")
          y <- tapply(atab2$y[wBool], atab2$SELmax_PAS_fac[wBool], median, na.rm=T)
          lines(SELmax_base + as.numeric(names(y)), y, lwd=1.2, type="b", pch=16, col="red")
        }
      }
      
      # Add CAS
      
      for(w in whales) {
        
        wBool <- sBool & (atab2$SELmax_CAS > 0) & atab2$ID==w
        if(sum(wBool, na.rm=T)>0) {
          #lines(SELmax_base+atab2$SELmax_CAS[wBool]*10, atab2$y[wBool], lwd=1.2, type="b", pch=16, col="blue")
          y <- tapply(atab2$y[wBool], atab2$SELmax_CAS_fac[wBool], median, na.rm=T)
          lines(SELmax_base + as.numeric(names(y)), y, lwd=1.2, type="b", pch=16, col="blue")
        }
      }
      
      
      ### Sonar exposures - deep
      sBool <- (atab2$Angle >= d & atab2$Angle < (d+ deg_add)) & atab2$depth_first >= 200
      sBool <- sBool & atab2$SELmax>0
      
      plot(SELmax_base+atab2$SELmax_PAS[sBool], atab2$y[sBool], col=redt, 
           xlab=expression(paste(plain("SEL (dB re 1"), mu,plain(Pa)^2, "s)", sep="")),  ylab=ylabs[v], pch=16, cex=0.5,
           main=paste(">200 m, sonar angle: ", d, "-", d+deg_add," degrees", sep=""),
           ylim=mylims, xlim=c(80,170), las=1)
      points(SELmax_base+atab2$SELmax_CAS[sBool], atab2$y[sBool], col=bluet, pch=16, cex=0.5)
      
      # Add PAS
      for(w in whales) {
        
        wBool <- sBool & (atab2$SELmax_PAS > 0) & atab2$ID==w
        if(sum(wBool, na.rm=T)>0) {
          #lines(SELmax_base+atab2$SELmax_PAS[wBool]*10, atab2$y[wBool], lwd=1.2, type="b", pch=16, col="red")
          y <- tapply(atab2$y[wBool], atab2$SELmax_PAS_fac[wBool], median, na.rm=T)
          lines(SELmax_base + as.numeric(names(y)), y, lwd=1.2, type="b", pch=16, col="red")
        }
      }
      
      # Add CAS
      
      for(w in whales) {
        
        wBool <- sBool & (atab2$SELmax_CAS > 0) & atab2$ID==w
        if(sum(wBool, na.rm=T)>0) {
          #lines(SELmax_base+atab2$SELmax_CAS[wBool]*10, atab2$y[wBool], lwd=1.2, type="b", pch=16, col="blue")
          y <- tapply(atab2$y[wBool], atab2$SELmax_CAS_fac[wBool], median, na.rm=T)
          lines(SELmax_base + as.numeric(names(y)), y, lwd=1.2, type="b", pch=16, col="blue")
        }
      }
      
    }
    
    options(graphics.record=FALSE) 
    dev.off()
    
  }
  