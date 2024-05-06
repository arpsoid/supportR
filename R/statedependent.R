#' Matrix-analytic model of the queueing-inventory system with random-size inventory requirement
#'
#' TransientQIG computes the Laplace-Stieltjes Transform for the
#' performance measures in the queueing-inventory model in transient state
#'
#' @param s The value at which to compute the LST (gives steady-state average if s=0)
#' @return The LST for the various performance measures of the system transient state
#' @export
#' @examples
#' TransientQIG(0) # returns the performance in steady state in a system with lambda=1, mu1=1, mu2=2, p1=0.4
#' invlap(Vectorize(TransientQIG), 0, 80, 500) # returns the transient performance of the model
TransientQIG=function(s){
  tol=1e-15
  Sglobal=7
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
    numstatesn[y+1]=numstatesn[y]+partitions::P(y)
    npartmatrix=rbind(npartmatrix,t(apply(partitions::parts(y),2,charpart)))
    ynpartmatrix=rbind(ynpartmatrix,
                       cbind(y,npartmatrix))
  }

  sizeY=sum(numstatesn)
  forwind=function(v) which(apply(ynpartmatrix,1,function(x)all(v==x)))

  piof0=c(1,rep(0,sizeY-1))

  A1=Am1=A0=matrix(0,nrow=sizeY,ncol=sizeY)

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
    PIO=sum(q[1,]), # probability of idle orbit
    PES=q[1,1], # probability of idle system
    ERR=sum(q[1,] %*% solve(diag(sizeY)-R) %*% Gam), # effective replenishment rate
    ESR=sum(q[1,] %*% solve(diag(sizeY)-R) %*% M), # effective service rate
    EAS=sum(q[1,] %*% solve(diag(sizeY)-R) %*% (Arr0)), # effective arrival into service
    EAO=sum(q[1,] %*% solve(diag(sizeY)-R) %*% (A1)), # effective arrival into orbit
    ETAR=sum(q[1,] %*% solve(diag(sizeY)-R) %*% (A1+Arr0)) # effective total arrival rate
  ))
}
