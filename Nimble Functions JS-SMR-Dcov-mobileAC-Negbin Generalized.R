dHab1 <- nimbleFunction(
  run = function(x = double(1),pi.cell = double(1),cells = double(2),res = double(0),
                 dSS = double(2),xlim = double(1),ylim = double(1),z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super==1){
      # uniform over state space cell
      logProb.unif <- - 2*log(res)
      # cell likelihood
      cell <- cells[trunc(x[1]/res)+1, trunc(x[2]/res)+1]
      logProb.cell <- log(pi.cell[cell])
      logProb <- logProb.unif + logProb.cell
    }else{
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rHab1 <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(1),cells = double(2),res = double(0),
                 dSS = double(2),xlim = double(1),ylim = double(1),z.super = double(0)){
    returnType(double(1))
    if(z.super==1){
      #simulate cell from pi.cell then uniform within cell
      n.cells <- nimDim(pi.cell)[1]
      s.cell <- rcat(1, pi.cell[1:n.cells])
      x.sim <- runif(1, dSS[s.cell,1] - res/2, dSS[s.cell,1] + res/2)
      y.sim <- runif(1, dSS[s.cell,2] - res/2, dSS[s.cell,2] + res/2)
      return(c(x.sim, y.sim))
    }else{
      return(c(0,0))
    }
  }
)

dHabMove <- nimbleFunction(
  run = function(x = double(1),s.prev = double(1),use.dist = double(1),dSS = double(2),
                 cells = double(2),res = double(0),sigma.move = double(0),z.super = double(0),
                 log = integer(0)){
    returnType(double(0))
    if(z.super==1){
      #find cell with lookup table
      cell <- cells[trunc(x[1]/res)+1,trunc(x[2]/res)+1]
      #get cell bounds
      x.min <- dSS[cell,1] - res/2
      x.max <- dSS[cell,1] + res/2
      y.min <- dSS[cell,2] - res/2
      y.max <- dSS[cell,2] + res/2
      # cell selection logProb
      logProb.cell <- log(use.dist[cell])
      # continuous within-cell location logProb
      logProb.x <- dnorm(x[1],s.prev[1],sigma.move, log=TRUE) -
        log(pnorm(x.max,s.prev[1],sigma.move) - pnorm(x.min,s.prev[1],sigma.move))
      logProb.y <- dnorm(x[2],s.prev[2],sigma.move,log=TRUE) -
        log(pnorm(y.max,s.prev[2],sigma.move) - pnorm(y.min,s.prev[2],sigma.move))
      logProb <- logProb.cell + logProb.x + logProb.y
    }else{
      logProb <- 0
    }
    if(log){
      return(logProb)
    }else{
      return(exp(logProb))
    }
  }
)

rHabMove <- nimbleFunction(
  run = function(n = integer(0),s.prev = double(1),use.dist = double(1),dSS = double(2),
                 cells = double(2),res = double(0),sigma.move = double(0),z.super = double(0)) {
    returnType(double(1))
    if(z.super==1){
      n.cells <- nimDim(dSS)[1]
      s.cell <- rcat(1,use.dist[1:n.cells])
      x.min <- dSS[s.cell,1] - res/2
      x.max <- dSS[s.cell,1] + res/2
      y.min <- dSS[s.cell,2] - res/2
      y.max <- dSS[s.cell,2] + res/2
      F.a <- pnorm(x.min,s.prev[1],sigma.move)
      F.b <- pnorm(x.max,s.prev[1],sigma.move)
      x.sim <- qnorm(runif(1,F.a,F.b),s.prev[1],sigma.move)
      F.a <- pnorm(y.min,s.prev[2],sigma.move)
      F.b <- pnorm(y.max,s.prev[2],sigma.move)
      y.sim <- qnorm(runif(1,F.a,F.b),s.prev[2],sigma.move)
      return(c(x.sim,y.sim))
    }else{
      return(c(0,0))
    }
  }
)


getAvail <- nimbleFunction(
  run = function(s = double(1),sigma=double(0),res=double(0),x.vals=double(1),y.vals=double(1),
                 n.cells.x=integer(0),n.cells.y=integer(0),z.super=double(0)) {
    returnType(double(1))
    if(z.super==1){
      avail.dist.x <- rep(0,n.cells.x)
      avail.dist.y <- rep(0,n.cells.y)
      delta <- 1e-8 #this sets the degree of trimming used to get individual availability distributions
      x.limits <- qnorm(c(delta,1-delta),mean=s[1],sd=sigma)
      y.limits <- qnorm(c(delta,1-delta),mean=s[2],sd=sigma)
      #convert to grid edges instead of centroids
      x.vals.edges <- c(x.vals - res/2, x.vals[n.cells.x]+0.5*res)
      y.vals.edges <- c(y.vals - res/2, y.vals[n.cells.y]+0.5*res)
      #trim in x and y direction
      if(x.vals.edges[1]<x.limits[1]){
        x.start <- floor((x.limits[1] - x.vals.edges[1]) / res) + 1
      }else{
        x.start <- 1
      }
      if(x.vals.edges[n.cells.x]>x.limits[2]){
        x.stop <- ceiling((x.limits[2] - x.vals.edges[1]) / res)
      }else{
        x.stop <- n.cells.x
      }
      if(y.vals.edges[1]<y.limits[1]){
        y.start <- floor((y.limits[1] - y.vals.edges[1]) / res) + 1
      }else{
        y.start <- 1
      }
      if(y.vals.edges[n.cells.y]>y.limits[2]){
        y.stop <- ceiling((y.limits[2] - y.vals.edges[1]) / res)
      }else{
        y.stop <- n.cells.y
      }
      pnorm.x <- rep(0,n.cells.x+1)
      pnorm.y <- rep(0,n.cells.y+1)
      #get pnorms
      for(l in x.start:(x.stop+1)){
        pnorm.x[l] <- pnorm(x.vals.edges[l],mean=s[1],sd=sigma)
      }
      for(l in y.start:(y.stop+1)){
        pnorm.y[l] <- pnorm(y.vals.edges[l],mean=s[2],sd=sigma)
      }
      for(l in (x.start):(x.stop)){
        avail.dist.x[l] <- pnorm.x[l+1] - pnorm.x[l]
      }
      for(l in (y.start):(y.stop)){
        avail.dist.y[l] <- pnorm.y[l+1] - pnorm.y[l]
      }
      avail.dist.tmp <- matrix(0,n.cells.x,n.cells.y)
      sum.dist <- 0
      for(i in x.start:x.stop){
        for(j in y.start:y.stop){
          avail.dist.tmp[i,j] <- avail.dist.x[i]*avail.dist.y[j]
          sum.dist <- sum.dist + avail.dist.tmp[i,j]
        }
      }
      avail.dist <- c(avail.dist.tmp)
      #if any probability mass is outside state space, normalize
      if(sum.dist<1){
        avail.dist <- avail.dist/sum.dist
      }
    }else{
      n.cells <- n.cells.x*n.cells.y
      avail.dist <- rep(0,n.cells)
    }
    return(avail.dist)
  }
)

getUse <- nimbleFunction(
  run = function(rsf = double(1),avail.dist=double(1),z.super = double(0)){
    returnType(double(1))
    if(z.super==1){
      use.dist <- rsf*avail.dist
      use.dist <- use.dist/sum(use.dist)
    }else{
      n.cells <- nimDim(rsf)[1]
      use.dist <- rep(0,n.cells)
    }
    return(use.dist)
  }
)

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

#this is used to restrict likelihood evaluation to only the primary occasions relevant for survival for each individual
dSurvival <- nimbleFunction(
  run = function(x = double(1), phi = double(1), z.start = double(0), z.stop = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    logProb <- 0
    if(z.super==1){
      n.primary <- length(phi)+1
      #extract first and last survival event primary occasions
      surv.start <- z.start+1
      surv.stop <- z.stop+1 #count death events, first z[i,]=0
      if(surv.start <= n.primary){ #if surv.start beyond last primary occasion, no survival events, logProb=0
        if(surv.stop > n.primary){ #but can't survive past n.primary
          surv.stop <- n.primary 
        }
        for(g in surv.start:surv.stop){ #sum logprob over survival event primary occasions
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
    n.primary <- length(phi)+1
    return(rep(0,n.primary))
  }
)

GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0), z.super=double(0)){ 
    returnType(double(1))
    if(z.super==0 | z.super==1&z==0){
      return(rep(0,J)) #skip calculation if not is superpop, or in superpop, but not alive in this primary occasion
    }
    if(z==1){ #otherwise calculate
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dNBVector <- nimbleFunction(
  run = function(x = double(1), mu = double(1), theta.d = double(0), K1D = double(1), z = double(0), z.super = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z.super*z==0){
      if(sum(x)>0){
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      p <- theta.d/(theta.d+mu)
      logProb <- sum(dnbinom(x,prob=p,size=theta.d*K1D,log=TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rNBVector <- nimbleFunction(
  run = function(n = integer(0), mu = double(1), theta.d = double(0), K1D = double(1), z = double(0), z.super = double(0)) {
    returnType(double(1))
    J <- nimDim(mu)[1]
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

#custom multinomial distribution for open-pop SMR with primary occasion-specific marked status
#x is y.event[i,g,1:J,1:3]
#y.sight is total detections y.sight[i,g,1:J]
#mark.states = 1 if individual is marked in this primary occasion, 0 otherwise
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
    n.sight.g <- control$n.sight.g
    sight.g <- control$sight.g
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
    for(g in 1:n.sight.g){
      gg <- sight.g[g]
      if(n.samples[gg] > 0){
        Jg <- J.sight[gg]
        #precalculate log likelihoods for this primary occasion
        ll.y <- matrix(0,nrow=M,ncol=Jg)
        ll.y.event <- matrix(0,nrow=M,ncol=Jg)
        for(i in 1:M){
          for(j in 1:Jg){
            if(z[i,gg]==1){
              p <- model$theta.d[g]/(model$theta.d[g]+model$lam[i,gg,j])
              ll.y[i,j] <- dnbinom(y.sight[i,gg,j],prob=p,size=model$theta.d[g]*model$K1D.sight[gg,j],log=TRUE)
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
        ###update IDs within primary occasion gg
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
            #update Negbin likelihoods for the two affected individuals at this detector/primary occasion
            p <- model$theta.d[g]/(model$theta.d[g]+model$lam[swapped[1],gg,j.curr])
            ll.y.cand[swapped[1],j.curr] <- dnbinom(y.sight.cand[swapped[1],gg,j.curr],
                                                    prob=p,size=model$theta.d[g]*model$K1D.sight[gg,j.curr],log=TRUE)
            
            p <- model$theta.d[g]/(model$theta.d[g]+model$lam[swapped[2],gg,j.curr])
            ll.y.cand[swapped[2],j.curr] <- dnbinom(y.sight.cand[swapped[2],gg,j.curr],
                                                    prob=p,size=model$theta.d[g]*model$K1D.sight[gg,j.curr],log=TRUE)
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
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    y2D.fixed <- control$y2D
    mark.g <- control$mark.g
    sight.g <- control$sight.g
    n.mark.g <- control$n.mark.g
    n.sight.g <- control$n.sight.g
    z.super.ups <- control$z.super.ups
    n.primary <- control$n.primary
    xlim <- control$xlim
    ylim <- control$ylim
    x.vals <- control$x.vals
    y.vals <- control$y.vals
    cells <- control$cells
    dSS <- control$dSS
    n.cells <- control$n.cells
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    res <- control$res
    s.nodes <- control$s.nodes
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
    #Conditional SMR: observed-individual status is dynamic because latent-ID sightings
    #can move among individuals. y2D.fixed contains only fixed evidence (marking,
    #known-ID sightings, telemetry-alive states). Add the CURRENT ID allocations via
    #capcounts, which is recalculated by IDSamplerOpen before this sampler runs.
    y.req <- y2D.fixed
    z.obs <- rep(0,M)
    for(i in 1:M){
      for(g2 in 1:n.sight.g){
        gg <- sight.g[g2]
        if(model$capcounts[gg,i] > 0){
          y.req[i,gg] <- 1
        }
      }
      if(sum(y.req[i,]) > 0){
        z.obs[i] <- 1
      }
    }
    
    #precompute entry counts: initial cohort, recruitment cohorts, never-entered class
    entry.counts.curr <- rep(0,n.primary+1)
    for(i in 1:M){
      if(model$z.super[i]==1){
        z.start.i <- model$z.start[i]
        entry.counts.curr[z.start.i] <- entry.counts.curr[z.start.i] + 1
      }else{
        entry.counts.curr[n.primary+1] <- entry.counts.curr[n.primary+1] + 1
      }
    }
    
    #1) Detected-individual updates: z.start, then z.stop.
    #Detection status and the constraining occasions come from current y.req.
    
    #1a) z.start Gibbs update
    for(i in 1:M){
      if(z.obs[i]==1 & model$z.super[i]==1 & y.req[i,1]==0){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        N.curr <- model$N
        N.recruit.curr <- model$N.recruit
        dets <- which(y.req[i,]>0)
        first.det <- min(dets)
        lp.start <- rep(-Inf,n.primary)
        mark.before <- which(mark.g < first.det)
        sight.before <- which(sight.g < first.det)
        
        #Remove focal individual once from its current entry class. Conditional on the
        #remaining allocation, the candidate-specific inverse-multinomial term is n_g+1.
        entry.counts.minus <- entry.counts.curr
        entry.counts.minus[z.start.curr] <- entry.counts.minus[z.start.curr] - 1
        lp.N1.not1 <- 0
        
        for(g in 1:first.det){
          z.start.prop <- g
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[g:first.det] <- 1
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary]
          }
          model$z[i,] <<- z.prop
          
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr
          if(z.start.curr > 1){
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          
          #Only occasions before first detection can differ across z.start candidates.
          for(h in 1:(first.det-1)){
            model$calculate(ER.nodes[h])
          }
          if(length(mark.before)>0){
            for(h in 1:length(mark.before)){
              idx.g <- mark.before[h]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(h in 1:length(sight.before)){
              idx.g <- sight.before[h]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
          }
          
          #There are only two possible N[1] values: entry at 1 versus entry after 1.
          if(g==1){
            lp.N1 <- model$calculate(N.nodes[1])
          }else{
            if(g==2){
              lp.N1.not1 <- model$calculate(N.nodes[1])
            }
            lp.N1 <- lp.N1.not1
          }
          lp.N.recruit <- 0
          for(h in 1:(first.det-1)){
            lp.N.recruit <- lp.N.recruit + model$calculate(N.recruit.nodes[h])
          }
          lp.y.mark <- 0
          lp.y.sight <- 0
          if(length(mark.before)>0){
            for(h in 1:length(mark.before)){
              idx.g <- mark.before[h]
              lp.y.mark <- lp.y.mark + model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            #Conditional formulation: only this individual's y.sight nodes change.
            for(h in 1:length(sight.before)){
              idx.g <- sight.before[h]
              lp.y.sight <- lp.y.sight + model$calculate(y.sight.nodes[i+(idx.g-1)*M])
            }
          }
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          lp.prior <- log(entry.counts.minus[g]+1)
          lp.start[g] <- lp.N1 + lp.N.recruit + lp.y.mark + lp.y.sight +
            lp.surv + lp.tel.z.states + lp.prior
        }
        
        maxlp <- max(lp.start)
        prop.probs <- exp(lp.start-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.start.prop <- rcat(1,prop.probs)
        model$z.start[i] <<- z.start.curr
        
        if(z.start.prop != z.start.curr){
          model$z.start[i] <<- z.start.prop
          z.prop <- rep(0,n.primary)
          z.prop[z.start.prop:first.det] <- 1
          if(first.det < n.primary){
            z.prop[(first.det+1):n.primary] <- z.curr[(first.det+1):n.primary]
          }
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.recruit <<- N.recruit.curr
          if(z.start.curr > 1){
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          
          for(h in 1:(first.det-1)){
            model$calculate(ER.nodes[h])
          }
          if(length(mark.before)>0){
            for(h in 1:length(mark.before)){
              idx.g <- mark.before[h]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
            for(h in 1:length(mark.before)){
              idx.g <- mark.before[h]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(h in 1:length(sight.before)){
              idx.g <- sight.before[h]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
            for(h in 1:length(sight.before)){
              idx.g <- sight.before[h]
              model$calculate(y.sight.nodes[i+(idx.g-1)*M])
            }
          }
          model$calculate(N.nodes[1])
          for(h in 1:(first.det-1)){
            model$calculate(N.recruit.nodes[h])
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          
          mvSaved["z.start",1][i] <<- model[["z.start"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["N.recruit",1] <<- model[["N.recruit"]]
          mvSaved["ER",1] <<- model[["ER"]]
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              gg <- mark.g[mark.before[g2]]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
              }
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              gg <- sight.g[sight.before[g2]]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
              }
            }
          }
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
          if(length(mark.before)>0){
            for(g2 in 1:length(mark.before)){
              gg <- mark.g[mark.before[g2]]
              for(j in 1:J.mark[gg]){
                model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
              }
            }
          }
          if(length(sight.before)>0){
            for(g2 in 1:length(sight.before)){
              gg <- sight.g[sight.before[g2]]
              for(j in 1:J.sight[gg]){
                model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
              }
            }
          }
          model$calculate(N.nodes[1])
          for(h in 1:(first.det-1)){
            model$calculate(N.recruit.nodes[h])
          }
          if(length(mark.before)>0){
            for(h in 1:length(mark.before)){
              idx.g <- mark.before[h]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.before)>0){
            for(h in 1:length(sight.before)){
              idx.g <- sight.before[h]
              model$calculate(y.sight.nodes[i+(idx.g-1)*M])
            }
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    
    #1b) z.stop Gibbs update
    for(i in 1:M){
      if(z.obs[i]==1 & model$z.super[i]==1 & y.req[i,n.primary]==0){
        z.curr <- model$z[i,]
        z.stop.curr <- model$z.stop[i]
        N.curr <- model$N
        dets <- which(y.req[i,]>0)
        last.det <- max(dets)
        lp.stop <- rep(-Inf,n.primary)
        mark.after <- which(mark.g > last.det)
        sight.after <- which(sight.g > last.det)
        
        for(g in last.det:n.primary){
          model$z.stop[i] <<- g
          z.prop <- rep(0,n.primary)
          z.prop[last.det:g] <- 1
          z.prop[1:last.det] <- z.curr[1:last.det]
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          
          if(last.det < n.primary-1){
            for(h in (last.det+1):(n.primary-1)){
              model$calculate(ER.nodes[h])
            }
          }
          if(length(mark.after)>0){
            for(h in 1:length(mark.after)){
              idx.g <- mark.after[h]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(h in 1:length(sight.after)){
              idx.g <- sight.after[h]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
          }
          
          if(last.det < n.primary-1){
            lp.N.recruit <- 0
            for(h in (last.det+1):(n.primary-1)){
              lp.N.recruit <- lp.N.recruit + model$calculate(N.recruit.nodes[h])
            }
          }else{
            lp.N.recruit <- 0
          }
          lp.y.mark <- 0
          lp.y.sight <- 0
          if(length(mark.after)>0){
            for(h in 1:length(mark.after)){
              idx.g <- mark.after[h]
              lp.y.mark <- lp.y.mark + model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(h in 1:length(sight.after)){
              idx.g <- sight.after[h]
              lp.y.sight <- lp.y.sight + model$calculate(y.sight.nodes[i+(idx.g-1)*M])
            }
          }
          lp.surv <- model$calculate(z.nodes[i])
          lp.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          lp.stop[g] <- lp.N.recruit + lp.y.mark + lp.y.sight + lp.surv + lp.tel.z.states
        }
        
        maxlp <- max(lp.stop)
        prop.probs <- exp(lp.stop-maxlp)
        prop.probs <- prop.probs/sum(prop.probs)
        z.stop.prop <- rcat(1,prop.probs)
        model$z.stop[i] <<- z.stop.curr
        
        if(z.stop.prop != z.stop.curr){
          model$z.stop[i] <<- z.stop.prop
          z.prop <- rep(0,n.primary)
          z.prop[last.det:z.stop.prop] <- 1
          z.prop[1:last.det] <- z.curr[1:last.det]
          model$z[i,] <<- z.prop
          model$N <<- N.curr - z.curr + z.prop
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          
          if(last.det < n.primary-1){
            for(h in (last.det+1):(n.primary-1)){
              model$calculate(ER.nodes[h])
            }
            for(h in (last.det+1):(n.primary-1)){
              model$calculate(N.recruit.nodes[h])
            }
          }
          if(length(mark.after)>0){
            for(h in 1:length(mark.after)){
              idx.g <- mark.after[h]
              model$calculate(pd.nodes[i+(idx.g-1)*M])
            }
            for(h in 1:length(mark.after)){
              idx.g <- mark.after[h]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(h in 1:length(sight.after)){
              idx.g <- sight.after[h]
              model$calculate(lam.nodes[i+(idx.g-1)*M])
            }
            for(h in 1:length(sight.after)){
              idx.g <- sight.after[h]
              model$calculate(y.sight.nodes[i+(idx.g-1)*M])
            }
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
          
          mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
          mvSaved["z",1][i,] <<- model[["z"]][i,]
          mvSaved["N",1] <<- model[["N"]]
          mvSaved["N.survive",1] <<- model[["N.survive"]]
          mvSaved["ER",1] <<- model[["ER"]]
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              gg <- mark.g[mark.after[g2]]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
              }
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              gg <- sight.g[sight.after[g2]]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
              }
            }
          }
        }else{
          model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
          model[["z"]][i,] <<- mvSaved["z",1][i,]
          model[["N"]] <<- mvSaved["N",1]
          model[["N.survive"]] <<- mvSaved["N.survive",1]
          model[["ER"]] <<- mvSaved["ER",1]
          if(length(mark.after)>0){
            for(g2 in 1:length(mark.after)){
              gg <- mark.g[mark.after[g2]]
              for(j in 1:J.mark[gg]){
                model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
              }
            }
          }
          if(length(sight.after)>0){
            for(g2 in 1:length(sight.after)){
              gg <- sight.g[sight.after[g2]]
              for(j in 1:J.sight[gg]){
                model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
              }
            }
          }
          if(last.det < n.primary-1){
            for(h in (last.det+1):(n.primary-1)){
              model$calculate(N.recruit.nodes[h])
            }
          }
          if(length(mark.after)>0){
            for(h in 1:length(mark.after)){
              idx.g <- mark.after[h]
              model$calculate(y.mark.nodes[i+(idx.g-1)*M])
            }
          }
          if(length(sight.after)>0){
            for(h in 1:length(sight.after)){
              idx.g <- sight.after[h]
              model$calculate(y.sight.nodes[i+(idx.g-1)*M])
            }
          }
          model$calculate(z.nodes[i])
          model$calculate(tel.z.states.nodes[i])
        }
      }
    }
    
    #2) Undetected individuals currently in the superpopulation.
    #Joint MH proposal for entry occasion and complete survival history. The survival
    #history is proposed from the demographic survival model, so its proposal probability
    #cancels exactly with the dSurvival target term.
    for(i in 1:M){
      if(z.obs[i]==0 & model$z.super[i]==1){
        z.curr <- model$z[i,]
        z.start.curr <- model$z.start[i]
        z.stop.curr <- model$z.stop[i]
        
        recruit.probs.for <- c(model$lambda.y1,model$ER)
        recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
        z.start.prop <- rcat(1,recruit.probs.for)
        log.prop.for <- log(recruit.probs.for[z.start.prop])
        
        z.prop <- rep(0,n.primary)
        z.prop[z.start.prop] <- 1
        z.stop.prop <- z.start.prop
        if(z.start.prop < n.primary){
          for(g in (z.start.prop+1):n.primary){
            if(z.prop[g-1]==1){
              z.prop[g] <- rbinom(1,1,model$phi[i,g-1])
              if(z.prop[g]==1){
                z.stop.prop <- g
              }
            }
          }
        }
        
        if(z.start.prop != z.start.curr | z.stop.prop != z.stop.curr){
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- 0
          for(h in 1:n.mark.g){
            lp.initial.y.mark <- lp.initial.y.mark + model$getLogProb(y.mark.nodes[i+(h-1)*M])
          }
          lp.initial.y.sight <- 0
          for(h in 1:n.sight.g){
            lp.initial.y.sight <- lp.initial.y.sight + model$getLogProb(y.sight.nodes[i+(h-1)*M])
          }
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[i])
          
          model$z[i,] <<- z.prop
          model$z.start[i] <<- z.start.prop
          model$z.stop[i] <<- z.stop.prop
          model$N <<- model$N - z.curr + z.prop
          if(z.start.curr > 1){
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          if(z.start.prop > 1){
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          model$calculate(ER.nodes)
          for(h in 1:n.mark.g){
            model$calculate(pd.nodes[i+(h-1)*M])
          }
          for(h in 1:n.sight.g){
            model$calculate(lam.nodes[i+(h-1)*M])
          }
          
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- 0
          for(h in 1:n.mark.g){
            lp.proposed.y.mark <- lp.proposed.y.mark + model$calculate(y.mark.nodes[i+(h-1)*M])
          }
          lp.proposed.y.sight <- 0
          for(h in 1:n.sight.g){
            lp.proposed.y.sight <- lp.proposed.y.sight + model$calculate(y.sight.nodes[i+(h-1)*M])
          }
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[i])
          
          if(z.start.prop != z.start.curr){
            log.z.prior.ratio <- log(entry.counts.curr[z.start.prop]+1) -
              log(entry.counts.curr[z.start.curr])
          }else{
            log.z.prior.ratio <- 0
          }
          
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          
          lp.initial.total <- lp.initial.N + lp.initial.N.recruit + lp.initial.y.mark +
            lp.initial.y.sight + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.N.recruit + lp.proposed.y.mark +
            lp.proposed.y.sight + lp.proposed.tel.z.states
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.prop.back) -
            (lp.initial.total + log.prop.for)
          accept <- decide(log_MH_ratio)
          
          if(accept){
            #Synchronize the survival logProb only after acceptance; it was omitted from
            #the MH calculation because it cancels with the survival proposal.
            model$calculate(z.nodes[i])
            mvSaved["z.start",1][i] <<- model[["z.start"]][i]
            mvSaved["z.stop",1][i] <<- model[["z.stop"]][i]
            mvSaved["z",1][i,] <<- model[["z"]][i,]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["ER",1] <<- model[["ER"]]
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][i,gg,j] <<- model[["pd"]][i,gg,j]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][i,gg,j] <<- model[["lam"]][i,gg,j]
              }
            }
            entry.counts.prop <- entry.counts.curr
            entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
            entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
            entry.counts.curr <- entry.counts.prop
          }else{
            model[["z.start"]][i] <<- mvSaved["z.start",1][i]
            model[["z.stop"]][i] <<- mvSaved["z.stop",1][i]
            model[["z"]][i,] <<- mvSaved["z",1][i,]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["ER"]] <<- mvSaved["ER",1]
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              for(j in 1:J.mark[gg]){
                model[["pd"]][i,gg,j] <<- mvSaved["pd",1][i,gg,j]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              for(j in 1:J.sight[gg]){
                model[["lam"]][i,gg,j] <<- mvSaved["lam",1][i,gg,j]
              }
            }
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            for(h in 1:n.mark.g){
              model$calculate(y.mark.nodes[i+(h-1)*M])
            }
            for(h in 1:n.sight.g){
              model$calculate(y.sight.nodes[i+(h-1)*M])
            }
            model$calculate(tel.z.states.nodes[i])
            #z survival logProb was never recalculated for the proposal.
          }
        }
      }
    }
    
    #3) Superpopulation-size update. Candidate sets are the CURRENT undetected on/off
    #individuals. Because z.obs was rebuilt after the ID update, an individual that gains
    #or loses its last latent-ID sighting automatically changes eligibility here.
    #For mobile ACs, add moves propose the entire s trajectory from its exact movement
    #prior. The movement target/proposal terms therefore cancel exactly, just like the
    #survival target/proposal terms. We still update and synchronize s/avail.dist/use.dist.
    z.on <- rep(0,M)
    z.off <- rep(0,M)
    non.curr <- 0
    noff.curr <- 0
    for(i in 1:M){
      if(z.obs[i]==0){
        if(model$z.super[i]==1){
          non.curr <- non.curr + 1
          z.on[non.curr] <- i
        }else{
          noff.curr <- noff.curr + 1
          z.off[noff.curr] <- i
        }
      }
    }
    
    for(up in 1:z.super.ups){
      updown <- rbinom(1,1,0.5)
      
      if(updown==0){
        #Remove one currently entered, undetected individual.
        non.init <- non.curr
        if(non.init>0){
          pick.pos <- rcat(1,rep(1/non.init,non.init))
          pick <- z.on[pick.pos]
          z.start.curr <- model$z.start[pick]
          z.curr <- model$z[pick,]
          log.p.select.for <- log(1/non.init)
          
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- 0
          for(h in 1:n.mark.g){
            lp.initial.y.mark <- lp.initial.y.mark + model$getLogProb(y.mark.nodes[pick+(h-1)*M])
          }
          lp.initial.y.sight <- 0
          for(h in 1:n.sight.g){
            lp.initial.y.sight <- lp.initial.y.sight + model$getLogProb(y.sight.nodes[pick+(h-1)*M])
          }
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])
          
          model$N.super <<- model$N.super - 1
          model$z.super[pick] <<- 0
          model$z.start[pick] <<- 0
          model$z.stop[pick] <<- 0
          model$z[pick,] <<- rep(0,n.primary)
          model$N <<- model$N - z.curr
          if(z.start.curr > 1){
            model$N.recruit[z.start.curr-1] <<- model$N.recruit[z.start.curr-1] - 1
          }
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          model$calculate(ER.nodes)
          
          #The spatial process is gated by z.super. Keep the off-state trajectory at zero.
          for(g in 1:n.primary){
            model$s[pick,g,1:2] <<- c(0,0)
          }
          for(g in 1:(n.primary-1)){
            model$avail.dist[pick,g,1:n.cells] <<- rep(0,n.cells)
            model$use.dist[pick,g,1:n.cells] <<- rep(0,n.cells)
          }
          for(g in 1:n.primary){
            model$calculate(s.nodes[pick+(g-1)*M])
          }
          for(h in 1:n.mark.g){
            model$calculate(pd.nodes[pick+(h-1)*M])
          }
          for(h in 1:n.sight.g){
            model$calculate(lam.nodes[pick+(h-1)*M])
          }
          
          #The reverse add move proposes the removed individual's old entry cohort.
          #Its old survival history and old s trajectory are proposed from their exact
          #priors, so those target/proposal terms cancel.
          recruit.probs.back <- c(model$lambda.y1,model$ER)
          recruit.probs.back <- recruit.probs.back/sum(recruit.probs.back)
          log.prop.back <- log(recruit.probs.back[z.start.curr])
          
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- 0
          for(h in 1:n.mark.g){
            lp.proposed.y.mark <- lp.proposed.y.mark + model$calculate(y.mark.nodes[pick+(h-1)*M])
          }
          lp.proposed.y.sight <- 0
          for(h in 1:n.sight.g){
            lp.proposed.y.sight <- lp.proposed.y.sight + model$calculate(y.sight.nodes[pick+(h-1)*M])
          }
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])
          
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.curr] <- entry.counts.prop[z.start.curr] - 1
          entry.counts.prop[n.primary+1] <- entry.counts.prop[n.primary+1] + 1
          noff.back <- noff.curr + 1
          log.p.select.back <- log(1/noff.back)
          log.z.prior.ratio <- log(entry.counts.curr[n.primary+1]+1) -
            log(entry.counts.curr[z.start.curr])
          log.prop.for <- 0
          
          lp.initial.total <- lp.initial.N + lp.initial.N.recruit + lp.initial.y.mark +
            lp.initial.y.sight + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.N.recruit + lp.proposed.y.mark +
            lp.proposed.y.sight + lp.proposed.tel.z.states
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.p.select.for + log.prop.for)
          accept <- decide(log_MH_ratio)
          
          if(accept){
            model$calculate(z.nodes[pick])
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            mvSaved["s",1][pick,1:n.primary,1:2] <<- model[["s"]][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              mvSaved["avail.dist",1][pick,g,1:n.cells] <<- model[["avail.dist"]][pick,g,1:n.cells]
              mvSaved["use.dist",1][pick,g,1:n.cells] <<- model[["use.dist"]][pick,g,1:n.cells]
            }
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
            z.on[pick.pos] <- z.on[non.curr]
            z.on[non.curr] <- 0
            non.curr <- non.curr - 1
            noff.curr <- noff.curr + 1
            z.off[noff.curr] <- pick
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            model[["s"]][pick,1:n.primary,1:2] <<- mvSaved["s",1][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              model[["avail.dist"]][pick,g,1:n.cells] <<- mvSaved["avail.dist",1][pick,g,1:n.cells]
              model[["use.dist"]][pick,g,1:n.cells] <<- mvSaved["use.dist",1][pick,g,1:n.cells]
            }
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              for(j in 1:J.mark[gg]){
                model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              for(j in 1:J.sight[gg]){
                model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
              }
            }
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            for(h in 1:n.mark.g){
              model$calculate(y.mark.nodes[pick+(h-1)*M])
            }
            for(h in 1:n.sight.g){
              model$calculate(y.sight.nodes[pick+(h-1)*M])
            }
            for(g in 1:n.primary){
              model$calculate(s.nodes[pick+(g-1)*M])
            }
            model$calculate(tel.z.states.nodes[pick])
          }
        }
        
      }else{
        #Add one currently never-entered, undetected individual.
        noff.init <- noff.curr
        if(noff.init>0){
          pick.pos <- rcat(1,rep(1/noff.init,noff.init))
          pick <- z.off[pick.pos]
          log.p.select.for <- log(1/noff.init)
          
          lp.initial.N <- model$getLogProb(N.nodes[1])
          lp.initial.N.recruit <- model$getLogProb(N.recruit.nodes)
          lp.initial.y.mark <- 0
          for(h in 1:n.mark.g){
            lp.initial.y.mark <- lp.initial.y.mark + model$getLogProb(y.mark.nodes[pick+(h-1)*M])
          }
          lp.initial.y.sight <- 0
          for(h in 1:n.sight.g){
            lp.initial.y.sight <- lp.initial.y.sight + model$getLogProb(y.sight.nodes[pick+(h-1)*M])
          }
          lp.initial.tel.z.states <- model$getLogProb(tel.z.states.nodes[pick])
          
          recruit.probs.for <- c(model$lambda.y1,model$ER)
          recruit.probs.for <- recruit.probs.for/sum(recruit.probs.for)
          z.start.prop <- rcat(1,recruit.probs.for)
          log.prop.for <- log(recruit.probs.for[z.start.prop])
          model$z.start[pick] <<- z.start.prop
          
          #Propose survival history from its exact survival model.
          model$z[pick,] <<- 0
          model$z[pick,z.start.prop] <<- 1
          z.stop.prop <- z.start.prop
          if(z.start.prop < n.primary){
            for(g in (z.start.prop+1):n.primary){
              if(model$z[pick,g-1]==1){
                model$z[pick,g] <<- rbinom(1,1,model$phi[pick,g-1])
                if(model$z[pick,g]==1){
                  z.stop.prop <- g
                }
              }
            }
          }
          model$z.stop[pick] <<- z.stop.prop
          model$N.super <<- model$N.super + 1
          model$z.super[pick] <<- 1
          model$N <<- model$N + model$z[pick,]
          if(z.start.prop > 1){
            model$N.recruit[z.start.prop-1] <<- model$N.recruit[z.start.prop-1] + 1
          }
          model$N.survive <<- model$N[2:n.primary] - model$N.recruit
          model$calculate(ER.nodes)
          
          #Propose the complete mobile-AC trajectory from its exact spatial/movement prior.
          model$s[pick,1,1:2] <<- rHab1(1,pi.cell=model$pi.cell[1:n.cells],
                                        cells=cells[1:n.cells.x,1:n.cells.y],
                                        dSS=dSS[1:n.cells,1:2],res=res,
                                        xlim=xlim,ylim=ylim,z.super=1)
          for(g in 2:n.primary){
            model$avail.dist[pick,g-1,1:n.cells] <<- getAvail(
              s=model$s[pick,g-1,1:2],sigma=model$sigma.move[1],res=res,
              x.vals=x.vals,y.vals=y.vals,n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
            model$use.dist[pick,g-1,1:n.cells] <<- getUse(
              rsf=model$rsf[1:n.cells],avail.dist=model$avail.dist[pick,g-1,1:n.cells],z.super=1)
            model$s[pick,g,1:2] <<- rHabMove(1,s.prev=model$s[pick,g-1,1:2],
                                             use.dist=model$use.dist[pick,g-1,1:n.cells],
                                             dSS=dSS[1:n.cells,1:2],cells=cells[1:n.cells.x,1:n.cells.y],
                                             res=res,sigma.move=model$sigma.move[1],z.super=1)
          }
          for(g in 1:n.primary){
            model$calculate(s.nodes[pick+(g-1)*M])
          }
          for(h in 1:n.mark.g){
            model$calculate(pd.nodes[pick+(h-1)*M])
          }
          for(h in 1:n.sight.g){
            model$calculate(lam.nodes[pick+(h-1)*M])
          }
          
          lp.proposed.N <- model$calculate(N.nodes[1])
          lp.proposed.N.recruit <- model$calculate(N.recruit.nodes)
          lp.proposed.y.mark <- 0
          for(h in 1:n.mark.g){
            lp.proposed.y.mark <- lp.proposed.y.mark + model$calculate(y.mark.nodes[pick+(h-1)*M])
          }
          lp.proposed.y.sight <- 0
          for(h in 1:n.sight.g){
            lp.proposed.y.sight <- lp.proposed.y.sight + model$calculate(y.sight.nodes[pick+(h-1)*M])
          }
          lp.proposed.tel.z.states <- model$calculate(tel.z.states.nodes[pick])
          
          entry.counts.prop <- entry.counts.curr
          entry.counts.prop[z.start.prop] <- entry.counts.prop[z.start.prop] + 1
          entry.counts.prop[n.primary+1] <- entry.counts.prop[n.primary+1] - 1
          non.back <- non.curr + 1
          log.p.select.back <- log(1/non.back)
          log.z.prior.ratio <- log(entry.counts.curr[z.start.prop]+1) -
            log(entry.counts.curr[n.primary+1])
          log.prop.back <- 0
          
          #Survival and movement target/proposal terms cancel exactly.
          lp.initial.total <- lp.initial.N + lp.initial.N.recruit + lp.initial.y.mark +
            lp.initial.y.sight + lp.initial.tel.z.states
          lp.proposed.total <- lp.proposed.N + lp.proposed.N.recruit + lp.proposed.y.mark +
            lp.proposed.y.sight + lp.proposed.tel.z.states
          log_MH_ratio <- (lp.proposed.total + log.z.prior.ratio + log.p.select.back + log.prop.back) -
            (lp.initial.total + log.p.select.for + log.prop.for)
          accept <- decide(log_MH_ratio)
          
          if(accept){
            model$calculate(z.nodes[pick])
            mvSaved["z.start",1][pick] <<- model[["z.start"]][pick]
            mvSaved["z.stop",1][pick] <<- model[["z.stop"]][pick]
            mvSaved["z",1][pick,] <<- model[["z"]][pick,]
            mvSaved["z.super",1][pick] <<- model[["z.super"]][pick]
            mvSaved["N",1] <<- model[["N"]]
            mvSaved["N.survive",1] <<- model[["N.survive"]]
            mvSaved["N.recruit",1] <<- model[["N.recruit"]]
            mvSaved["N.super",1][1] <<- model[["N.super"]]
            mvSaved["ER",1] <<- model[["ER"]]
            mvSaved["s",1][pick,1:n.primary,1:2] <<- model[["s"]][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              mvSaved["avail.dist",1][pick,g,1:n.cells] <<- model[["avail.dist"]][pick,g,1:n.cells]
              mvSaved["use.dist",1][pick,g,1:n.cells] <<- model[["use.dist"]][pick,g,1:n.cells]
            }
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              for(j in 1:J.mark[gg]){
                mvSaved["pd",1][pick,gg,j] <<- model[["pd"]][pick,gg,j]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              for(j in 1:J.sight[gg]){
                mvSaved["lam",1][pick,gg,j] <<- model[["lam"]][pick,gg,j]
              }
            }
            entry.counts.curr <- entry.counts.prop
            z.off[pick.pos] <- z.off[noff.curr]
            z.off[noff.curr] <- 0
            noff.curr <- noff.curr - 1
            non.curr <- non.curr + 1
            z.on[non.curr] <- pick
          }else{
            model[["z.start"]][pick] <<- mvSaved["z.start",1][pick]
            model[["z.stop"]][pick] <<- mvSaved["z.stop",1][pick]
            model[["z"]][pick,] <<- mvSaved["z",1][pick,]
            model[["z.super"]][pick] <<- mvSaved["z.super",1][pick]
            model[["N"]] <<- mvSaved["N",1]
            model[["N.survive"]] <<- mvSaved["N.survive",1]
            model[["N.recruit"]] <<- mvSaved["N.recruit",1]
            model[["N.super"]] <<- mvSaved["N.super",1][1]
            model[["ER"]] <<- mvSaved["ER",1]
            model[["s"]][pick,1:n.primary,1:2] <<- mvSaved["s",1][pick,1:n.primary,1:2]
            for(g in 1:(n.primary-1)){
              model[["avail.dist"]][pick,g,1:n.cells] <<- mvSaved["avail.dist",1][pick,g,1:n.cells]
              model[["use.dist"]][pick,g,1:n.cells] <<- mvSaved["use.dist",1][pick,g,1:n.cells]
            }
            for(g2 in 1:n.mark.g){
              gg <- mark.g[g2]
              for(j in 1:J.mark[gg]){
                model[["pd"]][pick,gg,j] <<- mvSaved["pd",1][pick,gg,j]
              }
            }
            for(g2 in 1:n.sight.g){
              gg <- sight.g[g2]
              for(j in 1:J.sight[gg]){
                model[["lam"]][pick,gg,j] <<- mvSaved["lam",1][pick,gg,j]
              }
            }
            model$calculate(N.nodes[1])
            model$calculate(N.recruit.nodes)
            for(h in 1:n.mark.g){
              model$calculate(y.mark.nodes[pick+(h-1)*M])
            }
            for(h in 1:n.sight.g){
              model$calculate(y.sight.nodes[pick+(h-1)*M])
            }
            for(g in 1:n.primary){
              model$calculate(s.nodes[pick+(g-1)*M])
            }
            model$calculate(tel.z.states.nodes[pick])
          }
        }
      }
    }
    
    #Update stored log probabilities after all custom state changes.
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset=function(){})
)

truncGammaPoisSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model,mvSaved,target,control){
    calcNodes <- model$getDependencies(target)
    upper <- model$getBound(target,"upper")
    if(target=="gamma"){
      is.fixed.gamma <- TRUE
      g <- 1
      recruitNodes <- grep("^N.recruit\\[",model$getNodeNames(stochOnly=TRUE),value=TRUE)
      n.recruit <- length(recruitNodes)
    }else{
      is.fixed.gamma <- FALSE
      g <- as.integer(gsub("[^0-9]","",target))
      n.recruit <- 1
    }
  },
  run = function(){
    if(is.fixed.gamma){
      count <- 0
      rate <- 0
      for(j in 1:n.recruit){
        count <- count+model$N.recruit[j]
        rate <- rate+model$N[j]
      }
    }else{
      count <- model$N.recruit[g]
      rate <- model$N[g]
    }
    if(rate>0){
      shape <- count+1
      p.upper <- pgamma(upper,shape=shape,rate=rate)
      model[[target]] <<- qgamma(runif(1,0,p.upper),shape=shape,rate=rate)
    }else{
      model[[target]] <<- runif(1,0,upper)
    }
    model$calculate(calcNodes)
    copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
  },
  methods=list(reset=function(){})
)