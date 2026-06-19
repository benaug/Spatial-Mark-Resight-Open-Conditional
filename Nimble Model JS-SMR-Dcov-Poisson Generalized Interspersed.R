NimModel <- nimbleCode({
  #Density covariates
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  D.beta1 ~ dnorm(0,sd=10)
  # D.intercept <- exp(D.beta0)*cellArea
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells]) #separate this component so s's do not depend on D.intercept
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  
  ##Abundance##
  lambda.y1 <- D.intercept*pi.denom #Expected starting population size
  N[1] ~ dpois(lambda.y1) #Realized starting population size
  for(g in 2:n.primary){
    N[g] <- N.survive[g-1] + N.recruit[g-1] #yearly abundance
    #N.recruit and N.survive information also contained in z/z.start + z.stop
    #N.recruit has distributions assigned below, but survival distributions defined on z
  }
  N.super <- N[1] + sum(N.recruit[1:(n.primary-1)]) #size of superpopulation
  
  #Recruitment
  gamma.fixed ~ dunif(0,2)
  for(g in 1:(n.primary-1)){
    # gamma[g] ~ dunif(0,2) # yearly recruitment priors
    gamma[g] <- gamma.fixed
    ER[g] <- N[g]*gamma[g] #yearly expected recruits
    N.recruit[g] ~ dpois(ER[g]) #yearly realized recruits
  }
  
  #Individual covariates
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]]) #categorical likelihood for this cell, equivalent to zero's trick
  }
  
  #Survival (phi must have M x n.primary - 1 dimension for custom updates to work)
  #without individual or year effects, use for loop to plug into phi[i,g]
  phi.fixed ~ dunif(0,1)
  for(i in 1:M){
    for(g in 1:(n.primary-1)){ #plugging same individual phi's into each year for custom update
      phi[i,g] <- phi.fixed
    }
    #survival likelihood (bernoulli) that only sums from z.start to z.stop
    z[i,1:n.primary] ~ dSurvival(phi=phi[i,1:(n.primary-1)],z.start=z.start[i],z.stop=z.stop[i],z.super=z.super[i])
    #telemetry survival likelihood
    #fixes z states, 1=known alive, 0=known dead, NA=unknown
    #currently assume censoring is uninformative
    tel.z.states[i,1:n.primary] ~ dSurvivalTel(z=z[i,1:n.primary],z.super=z.super[i])
  }
  
  ##Observation Model##
  #sample type observation model priors (Dirichlet), fixed across years
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])
  sigma.fixed ~ dunif(0,10)
  for(g in 1:n.primary){ #sigma informed by data except in years with no capture effort and no telemetry
    # sigma[g] ~ dunif(0,10) #sigma varies by year, shared across methods
    sigma[g] <- sigma.fixed #sigma fixed across years, shared across methods
  }
  #Marking process
  for(g in 1:n.mark.years){
    p0[g] ~ dunif(0,1)
    for(i in 1:M){
      pd[i,mark.years[g],1:J.mark[mark.years[g]]] <- GetDetectionProb(
        s=s[i,1:2],
        X=X.mark[mark.years[g],1:J.mark[mark.years[g]],1:2],
        J=J.mark[mark.years[g]],
        sigma=sigma[mark.years[g]],
        p0=p0[g],
        z=z[i,mark.years[g]],
        z.super=z.super[i])
      
      y.mark[i,mark.years[g],1:J.mark[mark.years[g]]] ~ dBinomialVector(
        pd[i,mark.years[g],1:J.mark[mark.years[g]]],
        K1D=K1D.mark[mark.years[g],1:J.mark[mark.years[g]]],
        z=z[i,mark.years[g]],
        z.super=z.super[i])
    }
  }
  
  #Sighting process
  for(g in 1:n.sight.years){
    lam0[g] ~ dunif(0,15)
    for(i in 1:M){
      lam[i,sight.years[g],1:J.sight[sight.years[g]]] <- GetDetectionRate(
        s=s[i,1:2],
        X=X.sight[sight.years[g],1:J.sight[sight.years[g]],1:2],
        J=J.sight[sight.years[g]],
        sigma=sigma[sight.years[g]],
        lam0=lam0[g],
        z=z[i,sight.years[g]],
        z.super=z.super[i])
      
      for(k in 1:K.sight[sight.years[g]]){
        y.sight[i,sight.years[g],1:J.sight[sight.years[g]],k] ~ dPoissonVector(
          lam[i,sight.years[g],1:J.sight[sight.years[g]]] *
            K2D.sight[sight.years[g],1:J.sight[sight.years[g]],k],
          z=z[i,sight.years[g]],
          z.super=z.super[i])
        
        y.event[i,sight.years[g],1:J.sight[sight.years[g]],k,1:3] ~ dmultiOpen(
          y.sight=y.sight[i,sight.years[g],1:J.sight[sight.years[g]],k],
          mark.states=mark.states[i,sight.years[g],k],
          theta.marked=theta.marked[1:3],
          theta.unmarked=theta.unmarked[1:3],
          capcounts=capcounts[sight.years[g],i])
      }
    }
    
    capcounts[sight.years[g],1:M] <- Getcapcounts(
      ID=ID[sight.years[g],1:n.samples[sight.years[g]]],
      capcounts.ID=capcounts.ID[sight.years[g],1:M])
    
    n.cap[sight.years[g]] <- Getncap(capcounts=capcounts[sight.years[g],1:M])
  }
  
  #Telemetry
  for(i in 1:n.tel.inds){
    for(g in 1:n.tel.sessions[i]){
      locs[i,g,1:max.n.tel.locs,1:2] ~ dNormVector(
        s=s[tel.ID[i],1:2],
        sigma=sigma[tel.session[i,g]],
        n.locs.ind=n.locs.ind[i,g],
        max.n.tel.locs=max.n.tel.locs)
    }
  }
})

#custom updates:
#1) for marked individuals: update z.start, then update z.stop
#2) for unmarked individuals: update entire z vectors
#3) N.super/z.super update 
