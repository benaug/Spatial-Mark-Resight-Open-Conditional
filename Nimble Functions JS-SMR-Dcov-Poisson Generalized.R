#telemetry survival vector distribution
dSurvivalTel <- nimbleFunction(
  run = function(x = double(1), z = double(1), z.super = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      for(g in 1:length(z)){
        if(x[g]!=2){ #2 codes NA (not yet collared or uninformatively censored), 1 is alive, 0 is dead
          #z invalid if conflicts with tel.
          if(x[g] != z[g]){
            logProb <- -Inf
          }
        }
      }
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rSurvivalTel <- nimbleFunction(
  run = function(n = integer(0), z = double(1),z.super = double(0)) {
    returnType(double(1))
    return(rep(0,length(z)))
  }
)

#telemetry location vector distribution
dNormVector <- nimbleFunction(
  run = function(x = double(2), s = double(1), sigma = double(0), n.locs.ind = double(0),
                 max.n.tel.locs = double(0), log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(n.locs.ind>0){
      for(i in 1:n.locs.ind){
        logProb <- logProb + dnorm(x[i,1], mean = s[1], sd = sigma, log = TRUE)
        logProb <- logProb + dnorm(x[i,2], mean = s[2], sd = sigma, log = TRUE)
      }
    }
    return(logProb) 
  }
)

rNormVector <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sigma = double(0), n.locs.ind = double(0),
                 max.n.tel.locs = double(0)) {
    returnType(double(2))
    out <- matrix(0,nrow=max.n.tel.locs,ncol=2)
    return(out)
  }
)

dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
    returnType(double(0))
    return(0)
  }
)

#this is used to restrict likelihood evaluation to only the years relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.primary <- length(phi)+1
      #extract first and last survival event years
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.primary){ #if surv.start beyond last year, no survival events, logProb=0
        if(surv.stop > n.primary){ #but can't survive past n.primary
          surv.stop <- n.primary 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event years
          logProb <- logProb + dbinom(x[g], size = 1, p = phi[g-1], log = TRUE)
        }
      }
    }
    return(logProb)
  }
)

#make dummy random vector generator to make nimble happy
rSurvival <- nimbleFunction(
  run = function(n = integer(0),phi = double(1), z.start = double(0), z.stop = double(0),z.super = double(0)) {
    returnType(double(1))
    n.primary <- length(phi)
    return(rep(0,n.primary))
  }
)

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this year
    }
    if(z==1){ #otherwise calculate
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super*z==0){
      if(sum(x)>0){
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      return(sum(dpois(x,lambda=lam,log=TRUE)))
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0), lam = double(1), z = double(0), z.super = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

# Function to calculate detection rate, but skip when z=0
GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J))
    }else{
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)
#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K1D = double(1), z = double(0),z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==0 | z.super==1&z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, size = K1D, prob = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0),pd = double(1), K1D = double(1), z = double(0),z.super = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#custom multinomial distribution for open-pop SMR with year-specific marked status
#x is y.event[i,g,1:J,1:3]
#y.sight is total detections y.sight[i,g,1:J]
#mark.states = 1 if individual is marked in this year, 0 otherwise
dmultiOpen <- nimbleFunction(
  run=function(x=double(2), y.sight=double(1), mark.states=double(0),
               theta.marked=double(1), theta.unmarked=double(1),
               capcounts=double(0), log=integer(0)){
    returnType(double(0))
    J <- nimDim(y.sight)[1]
    logProb <- 0
    if(capcounts==0){
      return(0)
    }else{
      for(j in 1:J){
        if(y.sight[j]>0){
          if(mark.states==1){
            logProb <- logProb + dmulti(x[j,1:3], size=y.sight[j], prob=theta.marked[1:3], log=TRUE)
          }else{
            logProb <- logProb + dmulti(x[j,2:3], size=y.sight[j], prob=theta.unmarked[2:3], log=TRUE)
          }
        }
      }
      return(logProb)
    }
  }
)

#dummy random generator
rmultiOpen <- nimbleFunction(
  run=function(n=integer(0), y.sight=double(1), mark.states=double(0),
               theta.marked=double(1), theta.unmarked=double(1),
               capcounts=double(0)){
    returnType(double(2))
    J <- nimDim(y.sight)[1]
    out <- matrix(0,nrow=J,ncol=3)
    return(out)
  }
)

#calculates how many samples each individual is currently allocated.
Getcapcounts <- nimbleFunction(
  run = function(ID=double(1),capcounts.ID=double(1)){
    returnType(double(1))
    n.samples <- nimDim(ID)[1]
    capcounts <- capcounts.ID
    for(l in 1:n.samples){
      capcounts[ID[l]] <- capcounts[ID[l]] + 1
    }
    return(capcounts)
  }
)

#calculate number of captured individuals
Getncap <- nimbleFunction(
  run = function(capcounts=double(1)){
    returnType(double(0))
    M <- nimDim(capcounts)[1]
    nstate <- numeric(M, value = 0)
    for(i in 1:M){
      if(capcounts[i]>0){
        nstate[i] <- 1
      }
    }
    n.cap <- sum(nstate)
    return(n.cap)
  }
)

IDSamplerOpen <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    M <- control$M
    J.sight <- control$J.sight
    n.primary <- control$n.primary
    n.sight.years <- control$n.sight.years
    sight.years <- control$sight.years
    mark.states <- control$mark.states
    event.type <- control$event.type
    match <- control$match
    n.samples <- control$n.samples
    this.j <- control$this.j
    calcNodes <- model$getDependencies(c("y.sight","y.event","ID"))
  },
  run = function() {
    z <- model$z
    y.sight <- model$y.sight
    y.event <- model$y.event
    ID.curr <- model$ID
    
    for(g in 1:n.sight.years){
      gg <- sight.years[g]
      if(n.samples[gg] > 0){
        Jg <- J.sight[gg]
        #precalculate log likelihoods for this year
        ll.y <- matrix(0,nrow=M,ncol=Jg)
        ll.y.event <- matrix(0,nrow=M,ncol=Jg)
        for(i in 1:M){
          for(j in 1:Jg){
            if(z[i,gg]==1){
              ll.y[i,j] <- dpois(y.sight[i,gg,j],model$K1D.sight[gg,j]*model$lam[i,gg,j],log=TRUE)
              if(y.sight[i,gg,j]>0){
                if(mark.states[i,gg]==1){
                  ll.y.event[i,j] <- dmulti(y.event[i,gg,j,1:3],y.sight[i,gg,j],model$theta.marked[1:3],log=TRUE)
                }else{
                  if(y.event[i,gg,j,1]>0){
                    ll.y.event[i,j] <- -Inf
                  }else{
                    ll.y.event[i,j] <- dmulti(y.event[i,gg,j,1:3],y.sight[i,gg,j],model$theta.unmarked[1:3],log=TRUE)
                  }
                }
              }
            }else{
              if(y.sight[i,gg,j]>0){
                ll.y[i,j] <- -Inf
              }else{
                ll.y[i,j] <- 0
              }
              ll.y.event[i,j] <- 0
            }
          }
        }
        ll.y.cand <- ll.y
        ll.y.event.cand <- ll.y.event
        ID.cand <- ID.curr
        y.sight.cand <- y.sight
        y.event.cand <- y.event
        
        ###update IDs within year gg
        for(l in 1:n.samples[gg]){
          j.curr <- this.j[gg,l]
          e.curr <- event.type[gg,l]
          propprobs <- model$lam[1:M,gg,j.curr]
          for(ii in 1:M){
            if(!match[gg,l,ii] | z[ii,gg]==0){
              propprobs[ii] <- 0
            }
          }
          propprobs <- propprobs/sum(propprobs)
          ID.cand[gg,l] <- rcat(1,prob=propprobs)
          if(ID.cand[gg,l] != ID.curr[gg,l]){
            swapped <- c(ID.curr[gg,l],ID.cand[gg,l])
            forprob <- propprobs[swapped[2]]
            backprob <- propprobs[swapped[1]]
            #move one sample from old ID to proposed ID
            y.event.cand[swapped[1],gg,j.curr,e.curr] <- y.event[swapped[1],gg,j.curr,e.curr] - 1
            y.event.cand[swapped[2],gg,j.curr,e.curr] <- y.event[swapped[2],gg,j.curr,e.curr] + 1
            y.sight.cand[swapped[1],gg,j.curr] <- y.sight[swapped[1],gg,j.curr] - 1
            y.sight.cand[swapped[2],gg,j.curr] <- y.sight[swapped[2],gg,j.curr] + 1
            #update Poisson likelihoods for the two affected individuals at this detector/year
            ll.y.cand[swapped[1],j.curr] <- dpois(y.sight.cand[swapped[1],gg,j.curr],
                                                  model$lam[swapped[1],gg,j.curr]*model$K1D.sight[gg,j.curr],log=TRUE)
            ll.y.cand[swapped[2],j.curr] <- dpois(y.sight.cand[swapped[2],gg,j.curr],
                                                  model$lam[swapped[2],gg,j.curr]*model$K1D.sight[gg,j.curr],log=TRUE)
            #update event likelihood for old ID
            if(y.sight.cand[swapped[1],gg,j.curr]==0){
              ll.y.event.cand[swapped[1],j.curr] <- 0
            }else{
              if(mark.states[swapped[1],gg]==1){
                ll.y.event.cand[swapped[1],j.curr] <- dmulti(y.event.cand[swapped[1],gg,j.curr,1:3],
                                                             y.sight.cand[swapped[1],gg,j.curr],
                                                             model$theta.marked[1:3],log=TRUE)
              }else{
                if(y.event.cand[swapped[1],gg,j.curr,1]>0){
                  ll.y.event.cand[swapped[1],j.curr] <- -Inf
                }else{
                  ll.y.event.cand[swapped[1],j.curr] <- dmulti(y.event.cand[swapped[1],gg,j.curr,1:3],
                                                               y.sight.cand[swapped[1],gg,j.curr],
                                                               model$theta.unmarked[1:3],log=TRUE)
                }
              }
            }
            #update event likelihood for new ID
            if(y.sight.cand[swapped[2],gg,j.curr]==0){
              ll.y.event.cand[swapped[2],j.curr] <- 0
            }else{
              if(mark.states[swapped[2],gg]==1){
                ll.y.event.cand[swapped[2],j.curr] <- dmulti(y.event.cand[swapped[2],gg,j.curr,1:3],
                                                             y.sight.cand[swapped[2],gg,j.curr],
                                                             model$theta.marked[1:3],log=TRUE)
              }else{
                if(y.event.cand[swapped[2],gg,j.curr,1]>0){
                  ll.y.event.cand[swapped[2],j.curr] <- -Inf
                }else{
                  ll.y.event.cand[swapped[2],j.curr] <- dmulti(y.event.cand[swapped[2],gg,j.curr,1:3],
                                                               y.sight.cand[swapped[2],gg,j.curr],
                                                               model$theta.unmarked[1:3],log=TRUE)
                }
              }
            }
            #select-sample proposal correction
            focalprob <- y.event[swapped[1],gg,j.curr,e.curr]/n.samples[gg]
            focalbackprob <- y.event.cand[swapped[2],gg,j.curr,e.curr]/n.samples[gg]
            lp.initial <- sum(ll.y[swapped,j.curr]) + sum(ll.y.event[swapped,j.curr])
            lp.proposed <- sum(ll.y.cand[swapped,j.curr]) + sum(ll.y.event.cand[swapped,j.curr])
            log_MH_ratio <- (lp.proposed + log(backprob) + log(focalbackprob)) -
              (lp.initial + log(forprob) + log(focalprob))
            accept <- decide(log_MH_ratio)
            if(accept){
              y.event[swapped[1],gg,j.curr,e.curr] <- y.event.cand[swapped[1],gg,j.curr,e.curr]
              y.event[swapped[2],gg,j.curr,e.curr] <- y.event.cand[swapped[2],gg,j.curr,e.curr]
              y.sight[swapped[1],gg,j.curr] <- y.sight.cand[swapped[1],gg,j.curr]
              y.sight[swapped[2],gg,j.curr] <- y.sight.cand[swapped[2],gg,j.curr]
              ll.y[swapped[1],j.curr] <- ll.y.cand[swapped[1],j.curr]
              ll.y[swapped[2],j.curr] <- ll.y.cand[swapped[2],j.curr]
              ll.y.event[swapped[1],j.curr] <- ll.y.event.cand[swapped[1],j.curr]
              ll.y.event[swapped[2],j.curr] <- ll.y.event.cand[swapped[2],j.curr]
              ID.curr[gg,l] <- ID.cand[gg,l]
            }else{
              y.event.cand[swapped[1],gg,j.curr,e.curr] <- y.event[swapped[1],gg,j.curr,e.curr]
              y.event.cand[swapped[2],gg,j.curr,e.curr] <- y.event[swapped[2],gg,j.curr,e.curr]
              y.sight.cand[swapped[1],gg,j.curr] <- y.sight[swapped[1],gg,j.curr]
              y.sight.cand[swapped[2],gg,j.curr] <- y.sight[swapped[2],gg,j.curr]
              ID.cand[gg,l] <- ID.curr[gg,l]
            }
          }
        }
      }
    }
    model$y.sight <<- y.sight
    model$y.event <<- y.event
    model$ID <<- ID.curr
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods = list(reset=function(){})
)

#all z updates live here
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    M <- control$M
    # n.marked.all <- control$n.marked.all
    n.cap.all <- control$n.cap.all
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    mark.states <- control$mark.states
    tel.z.states <- control$tel.z.states
    y2D <- control$y2D
    mark.years <- control$mark.years
    sight.years <- control$sight.years
    n.mark.years <- control$n.mark.years
    n.sight.years <- control$n.sight.years
    z.super.ups <- control$z.super.ups
    n.primary <- control$n.primary
    z.nodes <- control$z.nodes
    tel.z.states.nodes <- control$tel.z.states.nodes
    y.mark.nodes <- control$y.mark.nodes
    y.sight.nodes <- control$y.sight.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    N.nodes <- control$N.nodes
    ER.nodes <- control$ER.nodes
    N.survive.nodes <- control$N.survive.nodes
    N.recruit.nodes <- control$N.recruit.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    #precompute entry counts
    entry.counts.curr <- rep(0,n.primary+1)
    for(g in 1:n.primary){
      entry.counts.curr[g] <- sum(model$z.start==g & model$z.super==1)
    }
    entry.counts.curr[n.primary + 1] <- sum(model$z.super==0)
    #precompute y.req - these are z's that cannot be turned off because they currently have allocated samples
    y.req <- y2D
    for(i in 1:M){
      if(model$z.super[i]==1){
        for(g in 1:n.sight.years){
          gg <- sight.years[g]
          # if(sum(model$y.sight[i,gg,1:J.sight[gg]]) > 0){
          #   y.req[i,gg] <- 1
          # }
          #faster to use capcounts, already summed to individual by year
          if(model$capcounts[gg,i] > 0){
            y.req[i,gg] <- 1
          }
        }
      }
    }
    
    #1) Detected guy updates: z.start, z.stop
    # Detected guys are dynamically updated, these updates happen if sum(y.req[i,]>0)
    # 1a) z start update (z.stop update below): Gibbs, compute full conditional
    for(i in 1:M){
      if(model$z.super[i]==1 & sum(y.req[i,])>0 & y.req[i,1]==0){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y.req[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.primary)
        i.idx.mark <- seq(i,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
        i.idx.sight <- seq(i,M*n.sight.years,M) #used to reference correct sighting process nodes
        
        for(g in 1:first.det){ #must be recruited in year with first detection or before
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[g:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          
          #update N, N.recruit, N.survive. These individuals always in superpopulation
          #1) Update N
          model$N <<- N.curr - z.curr + z.prop
          #2) Update N.recruit
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes when a z changes
          
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.mark.nodes[i.idx.mark]) + model$calculate(y.sight.nodes[i.idx.sight])
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          
          # Add the full multinomial coefficient prior log-prob for this proposed configuration
          entry.counts.prop <- entry.counts.curr
          #z.super always 1 for detected guys
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          lp.prior <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.tel.z.states + lp.prior
        }
        maxlp <- max(lp.start) #deal with overflow
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr #set back to original
        
        if(model$z.start[i]!=z.start.prop){#if proposal is same as current, no need to replace anything
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[model$z.start[i]:first.det] <- 1 #must be alive until first detection
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary] #fill in remaining current z values, keeping death event the same
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr #set back to original first
          if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
            model$N.recruit[z.start.curr-1] <<- N.recruit.curr[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
            model$N.recruit[z.start.prop-1] <<- N.recruit.curr[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes
          
          #update these logProbs
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.sight.nodes[i.idx.sight])
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            for(j in 1:J.sight[gg]){
              mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
            }
          }
          #recompute entry counts
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            for(j in 1:J.sight[gg]){
              model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.sight.nodes[i.idx.sight])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    
    #1b) z stop update (z.start update above): Gibbs, compute full conditional
    # Detected guys are dynamically updated, these updates happen if sum(y.req[i,]>0)
    for(i in 1:M){
      if(model$z.super[i]==1 & sum(y.req[i,])>0 & y.req[i,n.primary]==0){
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y.req[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.primary)
        i.idx.mark <- seq(i,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
        i.idx.sight <- seq(i,M*n.sight.years,M) #used to reference correct sighting process nodes (y.um, y.unk and lam nodes
        for(g in (last.det):n.primary){ #can't die on or before year of last detection
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          #update N, number of recruits does not change going backwards
          model$N <<- N.curr - z.curr + z.prop
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes when a z changes
          #get these logProbs
          lp.N1 <- model$calculate(N.nodes[1])
          lp.N.recruit <- model$calculate(N.recruit.nodes)
          lp.y <- model$calculate(y.mark.nodes[i.idx.mark]) +  model$calculate(y.sight.nodes[i.idx.sight])
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          #no prior term, z.stop update does not change it
          lp.stop[g] <- lp.N1 + lp.N.recruit + lp.y + lp.surv + lp.tel.z.states
        }
        maxlp <- max(lp.stop) #deal with overflow
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr #set back to original
        if(model$z.stop[i]!=z.stop.prop){#if proposal differs from current
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.primary)
          z.prop[last.det:model$z.stop[i]] <- 1 #must be alive between last detection and this z.stop
          z.prop[1:(last.det)] <- z.curr[1:(last.det)] #fill in remaining current z values, keeping death event the same
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
          model$calculate(lam.nodes[i.idx.sight]) #update lam nodes
          #update these logProbs
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.sight.nodes[i.idx.sight])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            for(j in 1:J.sight[gg]){
              mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            for(j in 1:J.sight[gg]){
              model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
            }
          }
          #set these logProbs back
          model$calculate(N.nodes[1])
          model$calculate(N.recruit.nodes)
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.sight.nodes[i.idx.sight])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    #2) undetected guy update. Only if in the superpopulation.
    # Undetected guys are dynamically updated, these updates happen if sum(y.req[i,]==0)
    # Metropolis-Hastings, Propose z vectors from priors
    #entry counts current after z.start update
    for(i in (n.cap.all+1):M){
      if(model$z.super[i]==1 & sum(y.req[i,])==0){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        i.idx.mark <- seq(i,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
        i.idx.sight <- seq(i,M*n.sight.years,M) #used to reference correct sighting process nodes (y.um, y.unk and lam nodes
        #get forwards recruitment probabilities
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        #get initial logProbs
        lp.initial.entry <- model$getLogProb(N.nodes[1])
        lp.initial.entry <- lp.initial.entry + model$getLogProb(N.recruit.nodes)
        lp.initial.y.mark <- model$getLogProb(y.mark.nodes[i.idx.mark])
        lp.initial.y.sight <- model$getLogProb(y.sight.nodes[i.idx.sight])
        lp.initial.surv <- model$getLogProb(z.nodes[i])
        lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[i])
        log.prior.curr <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr + 1)))

        #track proposal probs - survival is symmetric, but not recruitment and detection
        log.prop.for <- log.prop.back <- 0

        #simulate recruitment
        z.start.prop <- rcat(1,recruit.probs.for)
        z.prop <- rep(0,n.primary)
        z.prop[z.start.prop] <- 1
        log.prop.for <- log.prop.for + log(recruit.probs.for[z.start.prop])

        #simulate survival
        if(z.start.prop < n.primary){#if you don't recruit in final year
          for(g in (z.start.prop+1):n.primary){
            z.prop[g] <- rbinom(1,1,model$phi[i,g-1]*z.prop[g-1])
            log.prop.for <- log.prop.for + dbinom(z.prop[g],1,model$phi[i,g-1]*z.prop[g-1],log=TRUE)
          }
        }
        z.on.prop <- which(z.prop==1)
        z.stop.prop <- max(z.on.prop)
        model$z[i,] <<- z.prop
        model$z.start[i] <<- z.start.prop
        model$z.stop[i] <<- z.stop.prop

        #update N, N.recruit, N.survive only if individual is in superpopulation
        #1) Update N
        model$N <<- model$N - z.curr + z.prop
        #2) Update N.recruit
        if(z.start.curr > 1){ #if wasn't in pop in year 1 in current, remove recruit event
          model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
        }
        if(z.start.prop > 1){ #if wasn't in pop in year 1 in proposal, add recruit event
          model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
        }
        #3) Update N.survive
        model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g

        model$calculate(ER.nodes) #update ER when N updated
        model$calculate(pd.nodes[i.idx.mark]) #update pd nodes when a z changes
        model$calculate(lam.nodes[i.idx.sight]) #update lam nodes after z changes
        #get proposed logProbs
        lp.proposed.entry <- model$calculate(N.nodes[1])
        lp.proposed.entry <- lp.proposed.entry + model$calculate(N.recruit.nodes)
        lp.proposed.y.mark <- model$calculate(y.mark.nodes[i.idx.mark])
        lp.proposed.y.sight <- model$calculate(y.sight.nodes[i.idx.sight])
        lp.proposed.surv <- model$calculate(z.nodes[i])
        lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[i])

        # Full multinomial coefficient prior for proposed configuration
        entry.counts.prop <- entry.counts.curr
        entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
        entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
        log.prior.prop <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop + 1)))

        #get backwards proposal probs
        recruit.probs.back <- c(model$lambda.y1,model$ER)
        recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
        log.prop.back <- log.prop.back + log(recruit.probs.back[z.start.curr])
        if(z.start.curr < n.primary){#if you don't recruit in final year
          for(g in (z.start.curr+1):n.primary){
            log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[i,g-1]*z.curr[g-1],log=TRUE)
          }
        }
        lp.initial.total <- lp.initial.entry + lp.initial.y.mark + lp.initial.y.sight +
          lp.initial.surv + lp.initial.tel.z.states + log.prior.curr
        lp.proposed.total <- lp.proposed.entry + lp.proposed.y.mark + lp.proposed.y.sight +
          lp.proposed.surv + lp.proposed.tel.z.states + log.prior.prop

        #MH step
        log_MH_ratio <- (lp.proposed.total + log.prop.back) - (lp.initial.total + log.prop.for)
        # log_MH_ratio <- (lp.proposed) - (lp.initial)
        accept <- decide(log_MH_ratio)

        if(accept){
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            for(j in 1:J.sight[gg]){
              mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
            }
          }
          entry.counts.curr <- entry.counts.prop
        }else{
          model[["z.start"]][i] <<- mvSaved["z.start",1][i]
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["N.recruit"]] <<- mvSaved["N.recruit",1]
          model[["ER"]] <<- mvSaved["ER",1]
          for(g in 1:n.mark.years){
            gg <- mark.years[g]
            for(j in 1:J.mark[gg]){
              model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
            }
          }
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            for(j in 1:J.sight[gg]){
              model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
            }
          }
          #set these logProbs back
          model$calculate(N.recruit.nodes)
          model$calculate(N.nodes[1])
          model$calculate(y.mark.nodes[i.idx.mark])
          model$calculate(y.sight.nodes[i.idx.sight])
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    #3) update z.super: Metropolis-Hastings. only involves unmarked individuals
    #entry counts current coming out of undetected ind update
    for(up in 1:z.super.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z.super==1)
        non.init <- length(z.on)
        pick <- rcat(1,rep(1/non.init,non.init))
        pick <- z.on[pick]
        # prereject turning off known captured individuals or any individuals currently allocated sighting samples
        if(pick <= n.cap.all){
          reject <- TRUE
        }
        if(!reject){
          for(g in 1:n.sight.years){
            gg <- sight.years[g]
            if(model$capcounts[gg,pick] > 0){
              reject <- TRUE
            }
          }
        }
        if(!reject){
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]

          #p select off guy
          log.p.select.for <- log(1/non.init)
          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))
          pick.idx.mark <- seq(pick,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
          pick.idx.sight <- seq(pick,M*n.sight.years,M)

          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes[pick.idx.mark])
          lp.initial.y.sight <- model$getLogProb(y.sight.nodes[pick.idx.sight])
          lp.initial.surv <- model$getLogProb(z.nodes[pick])
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])

          # propose new N.super/z.super/z.start/z.stop
          model$N.super <<-  model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.primary)

          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N - z.curr
          #2) Update N.recruit
          if(z.start.curr > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary]-model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[pick.idx.mark]) #turn pd off
          model$calculate(lam.nodes[pick.idx.sight]) #turn lam off
        
          #Reverse proposal probs
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          if(z.start.curr < n.primary){
            for(g in (z.start.curr+1):n.primary){
              log.prop.back <- log.prop.back + dbinom(z.curr[g],1,model$phi[pick,g-1]*z.curr[g-1],log=TRUE)
            }
          }

          #get proposed logProbs for N, N.recruit, and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- model$calculate(y.mark.nodes[pick.idx.mark]) #will always be 0
          lp.proposed.y.sight <- model$calculate(y.sight.nodes[pick.idx.sight]) #will always be 0
          lp.proposed.surv <- model$calculate(z.nodes[pick]) #will always be 0
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick]) #will always be 0

          lp.initial.total <- lp.initial.N + lp.initial.y.mark + lp.initial.y.sight + 
            lp.initial.N.recruit + lp.initial.surv + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.y.mark + lp.proposed.y.sight +
            lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.tel.z.states

          #backwards prior and select probs
          #move from class z.start.curr in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.primary + 1] <- entry.counts.prop[n.primary + 1] + 1

          #p select on guy
          noff.back <- sum(model$z.super == 0)
          log.p.select.back <- log(1/noff.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.for <- 0
          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)

          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              for(j in 1:J.sight[gg]){
                model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
              }
            }
            #set these logProbs back
            model$calculate(y.mark.nodes[pick.idx.mark])
            model$calculate(y.sight.nodes[pick.idx.sight])
            model$calculate(z.nodes[pick])
            model$calculate(tel.z.states.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
          }
        }
      }else{#add
        if(model$N.super[1] < M){ #cannot update if z.super maxed out. Need to raise M
          z.off <- which(model$z.super==0)
          noff.init <- length(z.off)
          pick <- rcat(1,rep(1/noff.init,noff.init)) #select one of these individuals
          pick <- z.off[pick]
          pick.idx.mark <- seq(pick,M*n.mark.years,M) #used to reference correct marking process nodes (y.mark and pd nodes)
          pick.idx.sight <- seq(pick,M*n.sight.years,M)

          non.init <- sum(model$z.super == 1)

          #p select off guy
          log.p.select.for <- log(1/noff.init)

          #log multinomial coefficient prior
          log.z.prior.for <- - (lgamma(M+1) - sum(lgamma(entry.counts.curr+1)))

          #get initial logProbs (survival logProb does not change)
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes[pick.idx.mark]) #will always be 0
          lp.initial.y.sight <- model$getLogProb(y.sight.nodes[pick.idx.sight]) #will always be 0
          lp.initial.surv <- model$getLogProb(z.nodes[pick]) #will always be 0
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick]) #will always be 0

          # Propose new z.start for the new on individual
          recruit.probs.for <- c(model$lambda.y1,model$ER)
          recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
          z.start.prop <- rcat(1,recruit.probs.for)  # propose entry cohort
          log.prop.for <- log(recruit.probs.for[z.start.prop])
          model$z.start[pick] <<- z.start.prop

          # Simulate survival path
          model$z[pick,] <<- 0 # initialize to 0
          model$z[pick,z.start.prop] <<- 1
          if(z.start.prop < n.primary){
            for(g in (z.start.prop+1):n.primary){
              model$z[pick,g] <<- rbinom(1,1,model$phi[pick,g-1]*model$z[pick,g-1])
              log.prop.for <- log.prop.for + dbinom(model$z[pick,g],1,model$phi[pick,g-1]*model$z[pick,g-1],log=TRUE)
            }
          }
          # Update z.stop
          z.on.prop <- which(model$z[pick,] == 1)
          z.stop.prop <- max(z.on.prop)
          model$z.stop[pick] <<- z.stop.prop

          #propose new N/z
          model$N.super <<-  model$N.super + 1
          model$z.super[pick] <<- 1

          #update N, N.recruit, N.survive
          #1) Update N
          model$N <<- model$N + model$z[pick,]
          #2) Update N.recruit
          if(model$z.start[pick] > 1){ #if wasn't in pop in year 1
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          #3) Update N.survive
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit #survivors are guys alive in year g-1 minus recruits in this year g
          model$calculate(ER.nodes) #update ER when N updated
          model$calculate(pd.nodes[pick.idx.mark]) #turn pd on
          model$calculate(lam.nodes[pick.idx.sight]) #turn lam on
         
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- model$calculate(y.mark.nodes[pick.idx.mark])
          lp.proposed.y.sight <- model$calculate(y.sight.nodes[pick.idx.sight])
          lp.proposed.surv <- model$calculate(z.nodes[pick])
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])

          lp.initial.total <- lp.initial.N + lp.initial.y.mark + lp.initial.y.sight  +
            lp.initial.N.recruit + lp.initial.surv + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.y.mark + lp.proposed.y.sight  +
            lp.proposed.N.recruit + lp.proposed.surv + lp.proposed.tel.z.states

          #backwards prior and select probs
          #move from class g in z.super==0 to class g in z.super==1
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.prop[n.primary + 1] <- entry.counts.prop[n.primary + 1] - 1

          #p select on guy
          non.back <- sum(model$z.super == 1)
          log.p.select.back <- log(1/non.back)
          #log multinomial coefficient prior
          log.z.prior.back <- - (lgamma(M+1) - sum(lgamma(entry.counts.prop+1)))
          log.prop.back <- 0

          #MH step
          log_MH_ratio <- (lp.proposed.total + log.z.prior.back + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.z.prior.for + log.p.select.for + log.prop.for)

          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1] <<- model[["z.super"]]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]] <<- mvSaved["z.super",1]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            for(g in 1:n.mark.years){
              gg <- mark.years[g]
              for(j in 1:J.mark[gg]){
                model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
              }
            }
            for(g in 1:n.sight.years){
              gg <- sight.years[g]
              for(j in 1:J.sight[gg]){
                model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
              }
            }
            #set these logProbs back
            model$calculate(y.mark.nodes[pick.idx.mark])
            model$calculate(y.sight.nodes[pick.idx.sight])
            model$calculate(z.nodes[pick])
            model$calculate(tel.z.states.nodes[pick])
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
          }
        }
      }
    }
    
    #copy back to mySaved to update logProbs.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
