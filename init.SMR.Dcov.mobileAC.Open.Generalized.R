e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.mobileAC.Open.Generalized <- function(data,inits=NA,M=NA,obsmod=NA){
  if(M < (data$n.cap.all)+1) stop("M must be larger than the number of captured individuals plus at least one unmarked individual.")
  library(abind)
  
  n.marked <- data$n.marked
  n.marked.all <- data$n.marked.all
  n.cap.all <- data$n.cap.all
  n.primary <- data$n.primary
  
  mark.states <- matrix(0,M,n.primary)
  mark.states[1:n.marked.all,] <- data$mark.states
  
  #augment tel.z.states, code NA as 2
  tel.z.states <- matrix(2,M,n.primary)
  tel.z.states[1:n.marked.all,] <- data$tel.z.states
  tel.z.states[is.na(tel.z.states)] <- 2
  
  J.mark <- unlist(lapply(data$X.mark,nrow))
  J.sight <- unlist(lapply(data$X.sight,nrow))
  J.mark.max <- max(J.mark)
  J.sight.max <- max(J.sight)
  K.mark <- data$K.mark
  K.sight <- data$K.sight
  K.mark.max <- max(K.mark)
  K.sight.max <- max(K.sight)
  locs <- data$locs
  
  #augment marking data and pull out sighting data
  y.mark <- array(0,dim=c(M,n.primary,J.mark.max))
  y.mark[1:n.cap.all,,] <- data$y.mark
  y.mID <- data$y.mID
  y.mnoID <- data$y.mnoID
  y.um <- data$y.um
  y.unk <- data$y.unk
  
  #reformat trap/effort arrays
  ID.marked <- matrix(0,max(n.marked),n.primary)
  X.mark <- array(0,dim=c(n.primary,J.mark.max,2))
  K1D.mark <- matrix(0,n.primary,J.mark.max)
  X.sight <- array(0,dim=c(n.primary,J.sight.max,2))
  K1D.sight <- matrix(0,n.primary,J.sight.max)
  
  for(g in 1:n.primary){
    if(n.marked[g]>0){
      ID.marked[1:n.marked[g],g] <- data$ID.marked[[g]]
    }
    if(J.mark[g]>0){
      X.mark[g,1:J.mark[g],1:2] <- data$X.mark[[g]]
      K1D.mark[g,1:J.mark[g]] <- data$K1D.mark[[g]]
    }
    if(J.sight[g]>0){
      X.sight[g,1:J.sight[g],1:2] <- data$X.sight[[g]]
      K1D.sight[g,1:J.sight[g]] <- data$K1D.sight[[g]]
    }
  }
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  p0 <- inits$p0
  lam0 <- inits$lam0
  sigma <- inits$sigma
  sigma.move.init <- inits$sigma.move
  rsf.beta.init <- inits$rsf.beta
  D.beta1.init <- inits$D.beta1
  n.cells <- nrow(data$dSS)
  n.cells.x <- length(data$x.vals)
  n.cells.y <- length(data$y.vals)
  lambda.cell <- data$InSS*exp(D.beta1.init*data$D.cov)
  pi.cell <- lambda.cell/sum(lambda.cell)
  rsf <- data$InSS*exp(rsf.beta.init*data$D.cov)
  if(!(obsmod %in% c("poisson","negbin"))) stop("obsmod must be 'poisson' or 'negbin'.")
  if(obsmod=="negbin"){
    theta.d <- inits$theta.d
    if(length(theta.d)!=n.primary) stop("inits$theta.d must have length n.primary for obsmod='negbin'.")
  }
  
  #initialize known/captured/telemetry ACs using known evidence before latent allocation
  # has.mark <- apply(y.mark[1:n.cap.all,,],1,sum)>0
  # has.mID <- rep(FALSE,n.cap.all)
  # has.mID[1:n.marked.all] <- rowSums(y.mID)>0
  # has.tel <- (1:n.cap.all)%in%data$tel.ID
  # idx <- which(has.mark|has.mID|has.tel)

  #build conditional sighting histories
  y.sight <- array(0,dim=c(M,n.primary,J.sight.max))
  y.event <- array(0,dim=c(M,n.primary,J.sight.max,3))
  
  #known marked-ID sightings are category 1
  y.sight[1:n.marked.all,,] <- y.mID
  y.event[1:n.marked.all,,,1] <- y.mID
  
  #provisional mobile activity centers for assigning latent detections
  z.super.pre <- rep(1,M)
  s.pre <- array(0,dim=c(M,n.primary,2))
  avail.dist.pre <- use.dist.pre <- array(0,dim=c(M,n.primary-1,n.cells))
  on.inds <- which(z.super.pre==1)
  
  for(i in on.inds){
    obs2D <- rep(0,n.primary)
    for(g in 1:n.primary){
      obs2D[g] <- sum(y.mark[i,g,]) + sum(y.sight[i,g,])
    }
    if(!any(is.na(data$tel.ID))){
      tel.idx <- which(data$tel.ID==i)
      if(length(tel.idx)>0){
        for(tt in 1:data$n.tel.sessions[tel.idx]){
          gg <- data$tel.session[tel.idx,tt]
          if(!is.na(gg) && data$n.locs.ind[tel.idx,tt]>0){
            obs2D[gg] <- obs2D[gg] + data$n.locs.ind[tel.idx,tt]
          }
        }
      }
    }
    dets <- which(obs2D>0)
    if(length(dets)>0){
      first.det <- min(dets)
      last.det <- max(dets)
      for(g in dets){
        locs.g <- matrix(numeric(0),nrow=0,ncol=2)
        trapcaps <- which(y.mark[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.mark[[g]][trapcaps,,drop=FALSE])
        }
        trapcaps <- which(y.sight[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.sight[[g]][trapcaps,,drop=FALSE])
        }
        if(!any(is.na(data$tel.ID))){
          tel.idx <- which(data$tel.ID==i)
          if(length(tel.idx)>0){
            tel.g.idx <- which(data$tel.session[tel.idx,]==g)
            if(length(tel.g.idx)>0){
              nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
              if(nloc>0){
                locs.g <- rbind(locs.g,cbind(data$locs[tel.idx,tel.g.idx,1:nloc,1],data$locs[tel.idx,tel.g.idx,1:nloc,2]))
              }
            }
          }
        }
        mean.loc <- c(mean(locs.g[,1]),mean(locs.g[,2]))
        mean.loc[1] <- min(max(mean.loc[1],xlim[1]+1e-5),xlim[2]-1e-5)
        mean.loc[2] <- min(max(mean.loc[2],ylim[1]+1e-5),ylim[2]-1e-5)
        mean.cell.x <- trunc(mean.loc[1]/data$res) + 1
        mean.cell.y <- trunc(mean.loc[2]/data$res) + 1
        mean.cell <- data$cells[mean.cell.x,mean.cell.y]
        if(data$InSS[mean.cell]==1){
          s.pre[i,g,] <- mean.loc
        }else{
          dists <- sqrt((data$dSS[,1] - mean.loc[1])^2 + (data$dSS[,2] - mean.loc[2])^2)
          dists[data$InSS==0] <- Inf
          pick <- which.min(dists)
          s.pre[i,g,1] <- runif(1,data$dSS[pick,1] - data$res/2,data$dSS[pick,1] + data$res/2)
          s.pre[i,g,2] <- runif(1,data$dSS[pick,2] - data$res/2,data$dSS[pick,2] + data$res/2)
        }
      }
      if(first.det>1){
        for(g in (first.det-1):1){
          avail <- getAvail(s=s.pre[i,g+1,],sigma=sigma.move.init,res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.pre[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.pre[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      if(last.det<n.primary){
        for(g in (last.det+1):n.primary){
          avail <- getAvail(s=s.pre[i,g-1,],sigma=sigma.move.init,res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.pre[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.pre[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      if(last.det>first.det){
        for(g in first.det:(last.det-1)){
          if(!(g+1)%in%dets){
            avail <- getAvail(s=s.pre[i,g,],sigma=sigma.move.init,res=data$res,
                              x.vals=data$x.vals,y.vals=data$y.vals,
                              n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
            use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
            s.cell <- sample(n.cells,1,prob=use)
            s.pre[i,g+1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
            s.pre[i,g+1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
          }
        }
      }
    }else{
      s.cell <- sample(n.cells,1,prob=pi.cell)
      s.pre[i,1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
      s.pre[i,1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      for(g in 2:n.primary){
        avail <- getAvail(s=s.pre[i,g-1,],sigma=sigma.move.init,res=data$res,
                          x.vals=data$x.vals,y.vals=data$y.vals,
                          n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
        use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
        s.cell <- sample(n.cells,1,prob=use)
        s.pre[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
        s.pre[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      }
    }
    for(g in 2:n.primary){
      avail.dist.pre[i,g-1,] <- getAvail(s=s.pre[i,g-1,1:2],sigma=sigma.move.init,res=data$res,
                                         x.vals=data$x.vals,y.vals=data$y.vals,
                                         n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
      use.dist.pre[i,g-1,] <- getUse(rsf=rsf,avail.dist=avail.dist.pre[i,g-1,],z.super=1)
    }
  }
  
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      D.sight <- e2dist(s.pre[,g,],X.sight[g,1:J.sight[g],1:2])
      lamd <- lam0[g]*exp(-D.sight*D.sight/(2*sigma[g]*sigma[g]))
      for(j in 1:J.sight[g]){
        if(y.mnoID[g,j]>0){
          marked.inds <- which(mark.states[,g]==1 & tel.z.states[,g]!=0)
          prob <- lamd[marked.inds,j]
          if(sum(prob)<=0) stop(paste("No valid marked candidates for marked-no-ID sightings in year",g,"trap",j))
          prob <- prob/sum(prob)
          add <- as.numeric(rmultinom(1,y.mnoID[g,j],prob=prob))
          y.sight[marked.inds,g,j] <- y.sight[marked.inds,g,j] + add
          y.event[marked.inds,g,j,2] <- y.event[marked.inds,g,j,2] + add
        }
        if(y.um[g,j]>0){
          unmarked.inds <- which(mark.states[,g]==0 & tel.z.states[,g]!=0)
          prob <- lamd[unmarked.inds,j]
          if(sum(prob)<=0) stop(paste("No valid unmarked candidates for unmarked sightings in year",g,"trap",j))
          prob <- prob/sum(prob)
          add <- as.numeric(rmultinom(1,y.um[g,j],prob=prob))
          y.sight[unmarked.inds,g,j] <- y.sight[unmarked.inds,g,j] + add
          y.event[unmarked.inds,g,j,2] <- y.event[unmarked.inds,g,j,2] + add
        }
        if(y.unk[g,j]>0){
          avail.inds <- which(tel.z.states[,g]!=0)
          prob <- lamd[avail.inds,j]
          if(sum(prob)<=0) stop(paste("No valid candidates for unknown-status sightings in year",g,"trap",j))
          prob <- prob/sum(prob)
          add <- as.numeric(rmultinom(1,y.unk[g,j],prob=prob))
          y.sight[avail.inds,g,j] <- y.sight[avail.inds,g,j] + add
          y.event[avail.inds,g,j,3] <- y.event[avail.inds,g,j,3] + add
        }
      }
    }
  }
  
  #final mobile activity centers using allocated latent detections
  s.init <- array(0,dim=c(M,n.primary,2))
  avail.dist.init <- use.dist.init <- array(0,dim=c(M,n.primary-1,n.cells))
  on.inds <- 1:M
  for(i in on.inds){
    obs2D <- rep(0,n.primary)
    for(g in 1:n.primary){
      obs2D[g] <- sum(y.mark[i,g,]) + sum(y.sight[i,g,])
    }
    if(!any(is.na(data$tel.ID))){
      tel.idx <- which(data$tel.ID==i)
      if(length(tel.idx)>0){
        for(tt in 1:data$n.tel.sessions[tel.idx]){
          gg <- data$tel.session[tel.idx,tt]
          if(!is.na(gg) && data$n.locs.ind[tel.idx,tt]>0){
            obs2D[gg] <- obs2D[gg] + data$n.locs.ind[tel.idx,tt]
          }
        }
      }
    }
    dets <- which(obs2D>0)
    if(length(dets)>0){
      first.det <- min(dets)
      last.det <- max(dets)
      for(g in dets){
        locs.g <- matrix(numeric(0),nrow=0,ncol=2)
        trapcaps <- which(y.mark[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.mark[[g]][trapcaps,,drop=FALSE])
        }
        trapcaps <- which(y.sight[i,g,]>0)
        if(length(trapcaps)>0){
          locs.g <- rbind(locs.g,data$X.sight[[g]][trapcaps,,drop=FALSE])
        }
        if(!any(is.na(data$tel.ID))){
          tel.idx <- which(data$tel.ID==i)
          if(length(tel.idx)>0){
            tel.g.idx <- which(data$tel.session[tel.idx,]==g)
            if(length(tel.g.idx)>0){
              nloc <- data$n.locs.ind[tel.idx,tel.g.idx]
              if(nloc>0){
                locs.g <- rbind(locs.g,cbind(data$locs[tel.idx,tel.g.idx,1:nloc,1],data$locs[tel.idx,tel.g.idx,1:nloc,2]))
              }
            }
          }
        }
        mean.loc <- c(mean(locs.g[,1]),mean(locs.g[,2]))
        mean.loc[1] <- min(max(mean.loc[1],xlim[1]+1e-5),xlim[2]-1e-5)
        mean.loc[2] <- min(max(mean.loc[2],ylim[1]+1e-5),ylim[2]-1e-5)
        mean.cell.x <- trunc(mean.loc[1]/data$res) + 1
        mean.cell.y <- trunc(mean.loc[2]/data$res) + 1
        mean.cell <- data$cells[mean.cell.x,mean.cell.y]
        if(data$InSS[mean.cell]==1){
          s.init[i,g,] <- mean.loc
        }else{
          dists <- sqrt((data$dSS[,1] - mean.loc[1])^2 + (data$dSS[,2] - mean.loc[2])^2)
          dists[data$InSS==0] <- Inf
          pick <- which.min(dists)
          s.init[i,g,1] <- runif(1,data$dSS[pick,1] - data$res/2,data$dSS[pick,1] + data$res/2)
          s.init[i,g,2] <- runif(1,data$dSS[pick,2] - data$res/2,data$dSS[pick,2] + data$res/2)
        }
      }
      if(first.det>1){
        for(g in (first.det-1):1){
          avail <- getAvail(s=s.init[i,g+1,],sigma=sigma.move.init,res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.init[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.init[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      if(last.det<n.primary){
        for(g in (last.det+1):n.primary){
          avail <- getAvail(s=s.init[i,g-1,],sigma=sigma.move.init,res=data$res,
                            x.vals=data$x.vals,y.vals=data$y.vals,
                            n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
          use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
          s.cell <- sample(n.cells,1,prob=use)
          s.init[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
          s.init[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
        }
      }
      if(last.det>first.det){
        for(g in first.det:(last.det-1)){
          if(!(g+1)%in%dets){
            avail <- getAvail(s=s.init[i,g,],sigma=sigma.move.init,res=data$res,
                              x.vals=data$x.vals,y.vals=data$y.vals,
                              n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
            use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
            s.cell <- sample(n.cells,1,prob=use)
            s.init[i,g+1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
            s.init[i,g+1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
          }
        }
      }
    }else{
      s.cell <- sample(n.cells,1,prob=pi.cell)
      s.init[i,1,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
      s.init[i,1,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      for(g in 2:n.primary){
        avail <- getAvail(s=s.init[i,g-1,],sigma=sigma.move.init,res=data$res,
                          x.vals=data$x.vals,y.vals=data$y.vals,
                          n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
        use <- getUse(rsf=rsf,avail.dist=avail,z.super=1)
        s.cell <- sample(n.cells,1,prob=use)
        s.init[i,g,1] <- runif(1,data$dSS[s.cell,1] - data$res/2,data$dSS[s.cell,1] + data$res/2)
        s.init[i,g,2] <- runif(1,data$dSS[s.cell,2] - data$res/2,data$dSS[s.cell,2] + data$res/2)
      }
    }
    for(g in 2:n.primary){
      avail.dist.init[i,g-1,] <- getAvail(s=s.init[i,g-1,1:2],sigma=sigma.move.init,res=data$res,
                                          x.vals=data$x.vals,y.vals=data$y.vals,
                                          n.cells.x=n.cells.x,n.cells.y=n.cells.y,z.super=1)
      use.dist.init[i,g-1,] <- getUse(rsf=rsf,avail.dist=avail.dist.init[i,g-1,],z.super=1)
    }
  }
  
  #construct ID/event lists from y.event
  n.samples <- rep(0,n.primary)
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      n.samples[g] <- sum(y.mnoID[g,1:J.sight[g]]) + sum(y.um[g,1:J.sight[g]]) + sum(y.unk[g,1:J.sight[g]])
    }
  }
  n.samples.max <- max(n.samples)
  
  ID <- matrix(1,nrow=n.primary,ncol=n.samples.max)
  this.j <- matrix(0,nrow=n.primary,ncol=n.samples.max)
  event.type <- matrix(0,nrow=n.primary,ncol=n.samples.max)
  match <- array(FALSE,dim=c(n.primary,n.samples.max,M))
  
  for(g in 1:n.primary){
    if(n.samples[g]>0){
      idx.samp <- 1
      #marked no-ID samples: category 2, only currently marked individuals are valid matches
      marked.inds <- which(mark.states[,g]==1)
      for(i in marked.inds){
        for(j in 1:J.sight[g]){
          if(y.event[i,g,j,2]>0){
            for(l in 1:y.event[i,g,j,2]){
              ID[g,idx.samp] <- i
              this.j[g,idx.samp] <- j
              event.type[g,idx.samp] <- 2
              match[g,idx.samp,which(mark.states[,g]==1 & tel.z.states[,g]!=0)] <- TRUE
              idx.samp <- idx.samp + 1
            }
          }
        }
      }
      
      #unmarked samples: category 2, only currently unmarked individuals are valid matches
      unmarked.inds <- which(mark.states[,g]==0)
      for(i in unmarked.inds){
        for(j in 1:J.sight[g]){
          if(y.event[i,g,j,2]>0){
            for(l in 1:y.event[i,g,j,2]){
              ID[g,idx.samp] <- i
              this.j[g,idx.samp] <- j
              event.type[g,idx.samp] <- 2
              match[g,idx.samp,which(mark.states[,g]==0 & tel.z.states[,g]!=0)] <- TRUE
              idx.samp <- idx.samp + 1
            }
          }
        }
      }
      
      #unknown marked-status samples: category 3, any not-known-dead individual is valid
      for(i in 1:M){
        for(j in 1:J.sight[g]){
          if(y.event[i,g,j,3]>0){
            for(l in 1:y.event[i,g,j,3]){
              ID[g,idx.samp] <- i
              this.j[g,idx.samp] <- j
              event.type[g,idx.samp] <- 3
              match[g,idx.samp,which(tel.z.states[,g]!=0)] <- TRUE
              idx.samp <- idx.samp + 1
            }
          }
        }
      }
      
      if(idx.samp != (n.samples[g]+1)) stop(paste("Sample-list construction mismatch in year",g))
    }
  }
  
  #baseline capcounts from known marked-ID category-1 detections
  capcounts.ID <- matrix(0,n.primary,M)
  for(g in 1:n.primary){
    if(J.sight[g]>0){
      for(i in 1:n.marked.all){
        capcounts.ID[g,i] <- sum(y.event[i,g,1:J.sight[g],1])
      }
    }
  }
  
  #initialize z, start with observed/allocated guys
  z.init <- matrix(0,M,n.primary)
  z.start.init <- z.stop.init <- rep(0,M)

  #get y2D constraints for z.start and z.stop update
  y.mark2D <- apply(y.mark,c(1,2),sum)
  y.mID2D <- matrix(0,M,n.primary)
  y.mID2D[1:n.marked.all,] <- apply(y.mID,c(1,2),sum)
  y2D <- y.mark2D + y.mID2D
  #add telemetry states - using these instead of marked states since you can be marked and dead (how you observe telemetry death)
  for(i in 1:n.marked.all){
    idx <- which(tel.z.states[i,]==1)
    if(length(idx)>0){
      y2D[i,idx] <- 1
    }
  }
  #need to use this y2D to initialize z
  y.sight2D <- 1*(apply(y.sight,c(1,2),sum)>0)
  y2D.init <- 1*((y2D + y.sight2D)>0)
  
  z.init <- 1*(y2D.init>0)
  for(i in 1:M){
    det.idx <- which(y2D.init[i,]>0)
    if(length(det.idx)>0){
      z.start.init[i] <- min(det.idx)
      z.stop.init[i] <- max(det.idx)
      z.init[i,z.start.init[i]:z.stop.init[i]] <- 1
    }
  }
  z.super.init <- 1*(z.start.init>0)
  N.super.init <- sum(z.super.init)
  
  if(any(tel.z.states[z.init==1]==0)) stop("At least one z initialized to 1 when tel.z.states=0. Bug in initialization code.")
  
  #initialize N structures from z.init
  N.init <- colSums(z.init[z.super.init==1,,drop=FALSE])
  N.survive.init <- N.recruit.init <- rep(NA,n.primary-1)
  for(g in 2:n.primary){
    N.survive.init[g-1] <- sum(z.init[,g-1]==1 & z.init[,g]==1 & z.super.init==1)
    N.recruit.init[g-1] <- N.init[g] - N.survive.init[g-1]
  }
  
  #basic starting likelihood checks
  logProb <- matrix(0,M,n.primary-1)
  for(i in 1:M){
    if(z.super.init[i]==1){
      for(g in 2:n.primary){
        logProb[i,g-1] <- dHabMove(x=s.init[i,g,1:2],s.prev=s.init[i,g-1,1:2],
                                   use.dist=use.dist.init[i,g-1,1:n.cells],
                                   dSS=data$dSS[1:n.cells,1:2],cells=data$cells[1:n.cells.x,1:n.cells.y],
                                   res=data$res,sigma.move=sigma.move.init,z.super=1,log=TRUE)
      }
    }
  }
  if(!all(is.finite(logProb))){
    stop("Starting logProb for activity centers is not finite, raise sigma.move.init. If that doesnt work, you may need to modify model or initialization algorithm.")
  }
  for(g in 1:n.primary){
    if(J.mark[g]>0){
      D.mark <- e2dist(s.init[,g,],X.mark[g,1:J.mark[g],1:2])
      pd <- p0[g]*exp(-D.mark*D.mark/(2*sigma[g]*sigma[g]))
      logProb <- array(0,dim=c(M,J.mark[g]))
      for(i in 1:M){
        for(j in 1:J.mark[g]){
          logProb[i,j] <- dbinom(y.mark[i,g,j],size=K1D.mark[g,j],prob=pd[i,j],log=TRUE)
        }
      }
      if(!is.finite(sum(logProb))) stop(paste("Starting observation model likelihood not finite. Marking process, year",g))
    }
    if(J.sight[g]>0){
      D.sight <- e2dist(s.init[,g,],X.sight[g,1:J.sight[g],1:2])
      lamd <- lam0[g]*exp(-D.sight*D.sight/(2*sigma[g]*sigma[g]))
      logProb <- matrix(0,M,J.sight[g])
      for(i in 1:M){
        for(j in 1:J.sight[g]){
          if(obsmod=="poisson"){
            logProb[i,j] <- dpois(y.sight[i,g,j],lambda=lamd[i,j]*K1D.sight[g,j],log=TRUE)
          }else if(obsmod=="negbin"){
            logProb[i,j] <- dnbinom(y.sight[i,g,j],mu=lamd[i,j]*K1D.sight[g,j],size=theta.d[g]*K1D.sight[g,j],log=TRUE)
          }
        }
      }
      if(!is.finite(sum(logProb))) stop(paste("Starting observation model likelihood not finite. Sighting process, year",g))
      #check y.event sums back to y.sight
      if(any(apply(y.event[,g,1:J.sight[g],1:3],c(1,2),sum) != y.sight[,g,1:J.sight[g]])){
        stop(paste("y.event does not sum to y.sight in year",g))
      }
    }
  }
  
  
  return(list(s=s.init,z=z.init,N=N.init,N.survive=N.survive.init,N.recruit=N.recruit.init,
              N.super=N.super.init,z.start=z.start.init,z.stop=z.stop.init,z.super=z.super.init,
              K1D.mark=K1D.mark,K1D.sight=K1D.sight,n.marked=n.marked,n.marked.all=n.marked.all,
              ID.marked=ID.marked,X.mark=X.mark,X.sight=X.sight,mark.states=mark.states,
              tel.z.states=tel.z.states,y2D=y2D,
              y.mark=y.mark,y.sight=y.sight,y.event=y.event,
              ID=ID,this.j=this.j,event.type=event.type,match=match,
              n.samples=n.samples,n.samples.max=n.samples.max,capcounts.ID=capcounts.ID,
              y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              xlim=xlim,ylim=ylim))
}
