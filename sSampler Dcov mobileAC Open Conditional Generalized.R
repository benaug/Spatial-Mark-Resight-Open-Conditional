sSampler1 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    s.nodes <- control$s.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    y.mark.nodes <- control$y.mark.nodes
    y.sight.nodes <- control$y.sight.nodes
    calcNodes <- control$calcNodes
    adaptive <- extractControlElement(control,'adaptive',TRUE)
    adaptInterval <- extractControlElement(control,'adaptInterval',200)
    adaptFactorExponent <- extractControlElement(control,'adaptFactorExponent',0.8)
    scale <- extractControlElement(control,'scale',1)
    scaleOriginal <- scale
    timesRan <- 0
    timesAccepted <- 0
    timesAdapted <- 0
    scaleHistory <- c(0,0)
    acceptanceHistory <- c(0,0)
    if(nimbleOptions('MCMCsaveHistory')){
      saveMCMChistory <- TRUE
    }else{
      saveMCMChistory <- FALSE
    }
    optimalAR <- 0.44
    gamma1 <- 0
    if(adaptFactorExponent < 0) stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0) stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function(){
    z.super <- model$z.super[i]
    z <- model$z[i,g]
    if(z.super==1 & z==1){
      s.cand <- c(rnorm(1,model$s[i,g,1],scale),rnorm(1,model$s[i,g,2],scale))
      inbox <- s.cand[1] < xlim[2] & s.cand[1] > xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        lp.initial.s <- model$getLogProb(s.nodes)
        lp.initial.y.mark <- model$getLogProb(y.mark.nodes)
        lp.initial.y.sight <- model$getLogProb(y.sight.nodes)
        model$s[i,g,1:2] <<- s.cand
        lp.proposed.s <- model$calculate(s.nodes)
        if(J.mark[g] > 0){
          model$calculate(pd.nodes)
        }
        if(J.sight[g] > 0){
          model$calculate(lam.nodes)
        }
        lp.proposed.y.mark <- model$calculate(y.mark.nodes)
        lp.proposed.y.sight <- model$calculate(y.sight.nodes)
        lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.sight
        lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.sight
        log_MH_ratio <- lp.proposed - lp.initial
        accept <- decide(log_MH_ratio)
        if(accept){
          copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
        }else{
          copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
        }
        if(adaptive){
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump) timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory,timesAdapted)
          scaleHistory[timesAdapted] <<- scale
          setSize(acceptanceHistory,timesAdapted)
          acceptanceHistory[timesAdapted] <<- acceptanceRate
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {
      returnType(double(1))
      if(saveMCMChistory){
        return(scaleHistory)
      }else{
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1,0))
      }
    },
    getAcceptanceHistory = function() {
      returnType(double(1))
      if(saveMCMChistory){
        return(acceptanceHistory)
      }else{
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1,0))
      }
    },
    reset = function() {
      scale <<- scaleOriginal
      timesRan <<- 0
      timesAccepted <<- 0
      timesAdapted <<- 0
      if(saveMCMChistory){
        scaleHistory <<- c(0,0)
        acceptanceHistory <<- c(0,0)
      }
      gamma1 <<- 0
    }
  )
)

sSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    calcNodes <- control$calcNodes
    adaptive <- extractControlElement(control,'adaptive',TRUE)
    adaptInterval <- extractControlElement(control,'adaptInterval',200)
    adaptFactorExponent <- extractControlElement(control,'adaptFactorExponent',0.8)
    scale <- extractControlElement(control,'scale',1)
    scaleOriginal <- scale
    timesRan <- 0
    timesAccepted <- 0
    timesAdapted <- 0
    scaleHistory <- c(0,0)
    acceptanceHistory <- c(0,0)
    if(nimbleOptions('MCMCsaveHistory')){
      saveMCMChistory <- TRUE
    }else{
      saveMCMChistory <- FALSE
    }
    optimalAR <- 0.44
    gamma1 <- 0
    if(adaptFactorExponent < 0) stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0) stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function() {
    z.super <- model$z.super[i]
    z <- model$z[i,g]
    if(z.super==1 & z==0){
      s.cand <- c(rnorm(1,model$s[i,g,1],scale),rnorm(1,model$s[i,g,2],scale))
      inbox <- s.cand[1] < xlim[2] & s.cand[1] > xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        model_lp_initial <- model$getLogProb(calcNodes)
        model$s[i,g,1:2] <<- s.cand
        model_lp_proposed <- model$calculate(calcNodes)
        log_MH_ratio <- model_lp_proposed - model_lp_initial
        accept <- decide(log_MH_ratio)
        if(accept){
          copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
        }else{
          copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
        }
        if(adaptive){
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump) timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory,timesAdapted)
          scaleHistory[timesAdapted] <<- scale
          setSize(acceptanceHistory,timesAdapted)
          acceptanceHistory[timesAdapted] <<- acceptanceRate
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {
      returnType(double(1))
      if(saveMCMChistory){
        return(scaleHistory)
      }else{
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1,0))
      }
    },
    getAcceptanceHistory = function() {
      returnType(double(1))
      if(saveMCMChistory){
        return(acceptanceHistory)
      }else{
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1,0))
      }
    },
    reset = function() {
      scale <<- scaleOriginal
      timesRan <<- 0
      timesAccepted <<- 0
      timesAdapted <<- 0
      if(saveMCMChistory){
        scaleHistory <<- c(0,0)
        acceptanceHistory <<- c(0,0)
      }
      gamma1 <<- 0
    }
  )
)

sSampler3 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    i <- control$i
    xlim <- control$xlim
    ylim <- control$ylim
    jump.multiplier <- control$jump.multiplier
    sig.move.fixed <- control$sig.move.fixed
    calcNodes <- control$calcNodes
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    s.nodes <- control$s.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    y.mark.nodes <- control$y.mark.nodes
    y.sight.nodes <- control$y.sight.nodes
  },
  run = function(){
    if(model$z.super[i]==1){
      if(sig.move.fixed==TRUE){
        scale.jump <- jump.multiplier*model$sigma.move[1]
      }else{
        scale.jump <- jump.multiplier*model$sigma.move[i]
      }
      s.cand <- c(rnorm(1,model$s[i,g,1],scale.jump),
                  rnorm(1,model$s[i,g,2],scale.jump))
      inbox <- s.cand[1] < xlim[2] & s.cand[1] > xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        if(model$z[i,g]==0){
          model_lp_initial <- model$getLogProb(calcNodes)
          model$s[i,g,1:2] <<- s.cand
          model_lp_proposed <- model$calculate(calcNodes)
          log_MH_ratio <- model_lp_proposed - model_lp_initial
          accept <- decide(log_MH_ratio)
          if(accept){
            copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
          }else{
            copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
          }
        }else{
          lp.initial.s <- model$getLogProb(s.nodes)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes)
          lp.initial.y.sight <- model$getLogProb(y.sight.nodes)
          model$s[i,g,1:2] <<- s.cand
          lp.proposed.s <- model$calculate(s.nodes)
          if(J.mark[g]>0){
            model$calculate(pd.nodes)
          }
          if(J.sight[g]>0){
            model$calculate(lam.nodes)
          }
          lp.proposed.y.mark <- model$calculate(y.mark.nodes)
          lp.proposed.y.sight <- model$calculate(y.sight.nodes)
          lp.initial <- lp.initial.s + lp.initial.y.mark + lp.initial.y.sight
          lp.proposed <- lp.proposed.s + lp.proposed.y.mark + lp.proposed.y.sight
          log_MH_ratio <- lp.proposed - lp.initial
          accept <- decide(log_MH_ratio)
          if(accept){
            copy(from=model,to=mvSaved,row=1,nodes=calcNodes,logProb=TRUE)
          }else{
            copy(from=mvSaved,to=model,row=1,nodes=calcNodes,logProb=TRUE)
          }
        }
      }
    }
  },
  methods = list(
    reset = function() {}
  )
)
