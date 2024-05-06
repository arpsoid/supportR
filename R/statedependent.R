# This file contains matrix-analytic model and simulation model of the
# multiclass queueing-inventory model where customer classes require
# groups of items from inventory of various size, and each class has a dedicated
# server.


rm(list=ls())

#================================ Matrix-Analytic Model below===================
library(partitions)

tol=1e-15
Sglobal=7
perf=function(s){
  gl=list(S=Sglobal,#3, # number of classes = inventory size
          s=2,#1, # inventory parameter, reorder point
          gamma= 3, # replenishment rate
          psS=0.5, # type of replenishment
          lambda = 0.2, # arrival rate
          p=(1:Sglobal)/sum(1:Sglobal), # customer class arrival probability
          mu = (Sglobal:1)/5, # service rates vector
          eta = 2, # retrial rates vector
          phi=rep(0.1,Sglobal))

  charpart=function(v) sapply(1:gl$S,function(i)sum(v==i))

  numstatesn=1 # for y=0
  npartmatrix=matrix(rep(0,gl$S),nrow=1)
  ynpartmatrix=cbind(0,npartmatrix)
  for(y in 1:gl$S){ # number of possible class combinations in service zone with (y-1) inventory
    numstatesn[y+1]=numstatesn[y]+P(y)
    npartmatrix=rbind(npartmatrix,t(apply(parts(y),2,charpart)))
    ynpartmatrix=rbind(ynpartmatrix,
                       cbind(y,npartmatrix))
  }

  sizeY=sum(numstatesn)
  forwind=function(v) which(apply(ynpartmatrix,1,function(x)all(v==x)))

  piof0=c(1,rep(0,sizeY-1))

  A1=Am1=A0=matrix(0,nrow=sizeY,ncol=sizeY)

  #A1=diag(sum(gl$phi*gl$p*gl$lambda),nrow=sizeY)

  for(y in 1:sizeY){
    was=ynpartmatrix[y,]
    toshare=was[1]-sum(was[-1]*(1:gl$S)) # any free inventory?
    if(toshare>0)
      for(wh in 1:toshare){ # class of successfully retrying customer
        will=was
        will[wh+1]=will[wh+1]+1
        Am1[y,forwind(will)]=gl$p[wh]*gl$eta
      }
  }

  Arr0=matrix(0,nrow=sizeY,ncol=sizeY) # successful arrival
  for(y in 1:sizeY)
    for(wh in 1:gl$S) { # class of arriving customer
      was=ynpartmatrix[y,]
      will=was
      toshare=was[1]-sum(was[-1]*(1:gl$S)) # any free inventory?
      if(toshare>0 & wh <=toshare) { # yes, and arriving customer fits in it
        will[wh+1]=will[wh+1]+1
        Arr0[y,forwind(will)]=gl$p[wh]*gl$lambda
      } else { # no, or the customer does not fit - goes to orbit
        A1[y,y]=A1[y,y]+gl$p[wh]*gl$phi[wh]*gl$lambda
      }
    }

  Gam=matrix(0,nrow=sizeY,ncol=sizeY) # replenishment
  for(y in 1:sizeY)
    if(ynpartmatrix[y,1]<=gl$s)
    {
      was=ynpartmatrix[y,]
      will=was
      will[1]=gl$S
      Gam[y,forwind(will)]=gl$psS*gl$gamma
      will[1]=was[1]+gl$S-gl$s
      Gam[y,forwind(will)]=(1-gl$psS)*gl$gamma
    }

  M=matrix(0,nrow=sizeY,ncol=sizeY) # service
  for(y in 1:sizeY)
  {
    was=ynpartmatrix[y,]
    for(k in 1:gl$S)
      if(was[k+1]>0){
        will=was
        will[k+1]=will[k+1]-1
        will[1]=will[1]-k
        M[y,forwind(will)]=gl$mu[k]
      }
  }

  A0=Gam+M+Arr0-A1-diag(unlist(apply(Am1+Gam+M+Arr0,1,sum)))

  A00=Gam+M+Arr0-A1-diag(unlist(apply(Gam+M+Arr0,1,sum)))

  A=A0+A1+Am1
  alpha=c(rep(0,sizeY-1),1) %*% solve(cbind(A[,-1],rep(1,sizeY)))
  if(sum(alpha %*% A1) >= sum(alpha %*% Am1)) stop("UNSTABLE")

  # now make it time-dependent
  A0=A0-diag(s,nrow=sizeY)
  A00=A00-diag(s,nrow=sizeY)
  A=A0+A1+Am1


  # PERFORMANCE
  Hm=solve(-A0) %*% A1
  Lm=solve(-A0) %*% Am1
  Tm=Hm
  G=Lm
  while(norm(abs(A1 %*% G %*% G + A0 %*% G + Am1),"m")>tol){
    Um=Hm%*%Lm+Lm%*%Hm
    Mm=Hm%*%Hm
    Hm=solve(diag(sizeY)-Um)%*%Mm
    Mm=Lm%*%Lm
    Lm=solve(diag(sizeY)-Um)%*%Mm
    G=G+Tm%*%Lm
    Tm=Tm%*%Hm
  }
  #H=A0
  #for (i in 1:20) {  G=solve(-H) %*% Am1; H=A0+A1 %*% G}

  R=-A1 %*% solve(A0+A1 %*% G)

  if(abs(s)<tol){
    q=c(rep(0,sizeY-1),1) %*% solve(cbind( (A00+R %*% Am1)[,-1], apply( solve(diag(sizeY)-R),1,sum) ))
    qnow=q
    return(sum(q[1,] %*% R %*% solve(diag(sizeY)-R)%*% solve(diag(sizeY)-R)))
    check=1-sum(q)

    while(abs(check)>tol){
      qnow= qnow %*% R # starts from q[1] onwards
      q=rbind(q, qnow) # steady-state probability vector for levels by rows
      check=check-sum(qnow) # calculate the tail probability
    }

  }  else{
    q=-piof0 %*% solve(A00+R %*% Am1)
    return(sum(q %*% R %*% solve(diag(sizeY)-R)%*% solve(diag(sizeY)-R)))
  }

  return(list(
    AIL=sum(apply(q,1,function(v)sum(v*ynpartmatrix[,1]))), # average inventory level
    AOL=sum(q[1,] %*% R %*% solve(diag(sizeY)-R)%*% solve(diag(sizeY)-R)), # average orbit level
    #AOL=sum(apply(q,1,sum)*(1:nrow(q)-1)) # average orbit level
    PIO=sum(q[1,]), # probability of idle orbit
    PES=q[1,1], # probability of idle system
    ERR=sum(q[1,] %*% solve(diag(sizeY)-R) %*% Gam), # effective replenishment rate
    #sum(apply(q,1,function(x) x %*% Gam)) # effective replenishment rate
    ESR=sum(q[1,] %*% solve(diag(sizeY)-R) %*% M), # effective service rate
    EAS=sum(q[1,] %*% solve(diag(sizeY)-R) %*% (Arr0)), # effective arrival into service
    EAO=sum(q[1,] %*% solve(diag(sizeY)-R) %*% (A1)), # effective arrival into orbit
    ETAR=sum(q[1,] %*% solve(diag(sizeY)-R) %*% (A1+Arr0)) # effective total arrival rate
  ))
  #AC1=sum(apply(q,1,function(v)sum(v*ynpartmatrix[,2]))) # average inventory level
}

library(pracma)
Fs=Vectorize(perf)
Li <- invlap(Fs, 0, 80, 500)
plot(Li,type="l",xlab="t",ylab=expression(A[OL](t)),lwd=2,log="x")
abline(h=perf(0),lty=2,lwd=2)

#========================================== Discrete-Event Simultation (GSMP)===

library(simulato)
# This file contains a class for the multiclass inventory-retrial model

#globalcounter=0

# State:
# 1:S -- number of this class customers in service zone
# S+1 -- orbit size
# S+2 -- inventory size
# Clock:
# 1:S -- residual service time of this class
# S+1 -- residual retrial time
# S+2 -- residual interarrival time
# S+3 -- residual replenishment time
queueinventory <- function(gl=list(S=3, # number of classes = inventory size
                                   s=1, # inventory parameter, reorder point
                                   gamma= 3, # replenishment rate
                                   psS=0.5, # type of replenishment
                                   lambda = 0.2, # arrival rate
                                   p=(1:3)/sum(1:3), # customer class arrival probability
                                   mu = (3:1)/5, # service rates vector
                                   eta = 2, # retrial rate
                                   phi=rep(0.1,3))) {# orbit joining probability
  #if(with(gl,sum(lambda*p/mu)+max(lambda*p/(lambda*p+eta))>1)) stop("Unstable")
  m=gsmp()
  class(m) <- append(class(m),"queueinventory")
  m$state <- c(1,rep(0,gl$S-1),0,gl$S) # start by serving a single class-1 customer
  m$gl <- gl
  m$clocks <- rep(Inf,gl$S+3)
  m$clocks[getActiveEvents(m)] <- getNewClocks(m,getActiveEvents(m))
  return(m)
}

isRegeneration.queueinventory <- function(m) {
  return(all(m$state[1:(m$gl$S+1)]==0) & m$state[m$gl$S+2]==m$gl$S)
}
getPerformance.queueinventory <- function(m){# what is the model performance measured at given point
  #m$state[m$gl$N+1]>0 # busy probability
  #m$state[m$gl$N+1]==1:m$gl$N # per class busy probability
  #m$state[m$gl$N+1]==0 & m$state[1:m$gl$N]==0 # per class idle with idle orbit
  #m$state[m$gl$N+1]==0 & m$state[1:m$gl$N]>0 # per class idle with busy orbit
  #return(m$state[m$gl$S+1]) # orbit-queue size
  return(all(m$state==0)) # idle probability
  #return(m$state[m$gl$S+2]) # inventory size #c(m$state[1],m$state[2],m$state[1]^2,m$state[2]^2,m$state[1]*m$state[2]))
}
getRates.queueinventory <- function(m){
  r <- rep(0,length(m$clocks))
  r[getActiveEvents(m)] <- 1
  return(r)
}
getNewGSMP.queueinventory <- function(m,e){
  S <- m$gl$S
  nm <- m
  fits=function(m,clarr) return(sum((1:S)*m$state[1:S])+clarr<=m$state[S+2])
  if(e==S+2){ # arrival
    clarr=sample.int(S,size=1,prob = m$gl$p)
    if(fits(m,clarr)) # can be accomodated in the service zone
      nm$state[clarr] <- m$state[clarr]+1 # this arrival goes to service
    else if(runif(1)<=m$gl$phi[clarr]) # decides to join the orbit
      nm$state[S+1] <- m$state[S+1]+1
  } else if(e<=S){ # departure from service zone
    nm$state[e]=m$state[e]-1
    nm$state[S+2]=m$state[S+2]-e
  } else if(e==S+1){ # retrial
    clarr=sample.int(S,size=1,prob = m$gl$p)
    if(fits(m,clarr)){ # successful
      nm$state[clarr] <- m$state[clarr]+1
      nm$state[e] <- m$state[e]-1
    }
  } else if(e==S+3){ # replenishment
    nm$state[S+2]=ifelse(runif(1)<m$gl$psS,m$gl$S,nm$state[S+2]+m$gl$S-m$gl$s)
    #globalcounter <<- globalcounter+1 # to count replenishments
  }
  return(nm)
}
getActiveEvents.queueinventory <- function(m){
  a=which(m$state[1:(m$gl$S+1)]>0)
  a=c(a,m$gl$S+2)
  if(m$state[m$gl$S+2]<=m$gl$s)
    a=c(a,m$gl$S+3)
  return(a)
}
getNewClocks.queueinventory <- function(m, e)  {# here we need to process the clocks based on the new events (their numbers)
  rates=c(m$gl$mu,m$gl$eta,m$gl$lambda,m$gl$gamma)[e]
  return(rexp(length(e),rate=rates))#ifelse(e==m$gl$S+2,m$gl$lambda, # arrival
  #ifelse(e<=m$gl$S,m$gl$mu[e], # service
  #ifelse(e==m$gl$S+3,m$gl$gamma, # replenishment
  #       m$gl$eta))   ))) # retrial
}
m=queueinventory()
getRegEst(trace(m,10000000))#[,"est"]
