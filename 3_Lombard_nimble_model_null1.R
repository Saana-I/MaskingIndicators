

sp.model_null1 <- nimbleCode({
    
    #model{
    
    ## Likelihood
    
    for(j in 2:N) { # loop across data

      # Observation model
      isCensored[j] ~ dinterval(zpeak0[j], censorLimitVec[j])
      zpeak0[j] ~ dnorm(y[j] + phi[1]*epsilon[j-1], tau.cor[ID[j]])
      
      # Serial correlation in predictor
      epsilon[j] <- (zpeak0[j] - y[j]) #- phi[1]*epsilon[j-1]
      
      # Linear predictor
      y[j] <- b0_ID[ID[j]] + mp[j]*(1-mr[j]) + sp[j]*(1-sr[j])
      
      # Masking potential
      mp[j] <- -mp1[1]*SeaState[j]
      
      # Masking release (proportion)
      logit(mr[j]) <- -10 + mr1[1]*depth_first[j] + mr2[1]*pitch_first[j]

      # Masking potential from sonar
      sp[j] <- sp1[1]*SELmax[j]
      
      # Masking release from sonar
      logit(sr[j]) <- -10 + sr1[1]*Angle[j]

    }
  
    ## 1st time step
    
    for(j in 1) {
      
      # Observation model
      isCensored[j] ~ dinterval(zpeak0[j], censorLimitVec[j])
      zpeak0[j] ~ dnorm(y[j], tau.cor[ID[j]])
      
      # Serial correlation in predictor
      epsilon[j] <- (zpeak0[j] - y[j])
      
      # Linear predictor
      y[j] <- b0_ID[ID[j]] + mp[j]*(1-mr[j]) + sp[j]*(1-sr[j])
      
      # Masking potential
      mp[j] <- -mp1[1]*SeaState[j]
      
      # Masking release (proportion)
      logit(mr[j]) <- -10 + mr1[1]*depth_first[j] + mr2[1]*pitch_first[j]
      
      # Masking potential from sonar
      sp[j] <- sp1[1]*SELmax[j]
      
      # Masking release from sonar
      logit(sr[j]) <- -10 + sr1[1]*Angle[j]
      
    }
    
    ## Priors

    for(k in 1:NW) {
      tau.cor[k] <- tau[k]#/(1-phi[1]*phi[1])
      tau[k] ~ dgamma(1, 1)
      b0_ID[k] ~ dnorm(b0[1], b0_tau[1])
    }
    
    phi[1] ~ dunif(-1,1)
    b0_tau[1] ~ dgamma(1, 1)
    b0[1] ~ dunif(100, 200)
    
    mp1[1] ~ dgamma(1, 1)   # Sea state

    mr1[1] ~ dgamma(1, 1)    # masking release f(depth)
    mr2[1] ~ dgamma(1, 1)    # masking release f(pitch)
    
    sp1[1] ~ dgamma(1, 1)   # CAS/PAS  
    
    sr1[1] ~ dgamma(1, 1)    # masking release f(Angle)
    
    })

