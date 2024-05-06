#' Matrix-analytic model of the 2-server multiserver job model
#'
#' TransientMJM2 computes the Laplace-Stieltjes Transform for the
#' average number of customers in the system in transient state
#'
#' @param s The value at which to compute the LST (gives steady-state average if s=0)
#' @param lambda The arrival rate
#' @param mu1 The service rate for class-1 customer (requires 1 server)
#' @param mu2 The service rate for class-2 customer (requires 2 servers)
#' @param p1 Probability that an arrival is class-1
#' @return The LST for the mean number of customers in the system transient state
#' @references R. Razumchik et al., A queueing system for performance evaluation of a markovian supercomputer model. 2023 doi: 10.14357/19922264230209.
#' @export
#' @examples
#' TransientMJM2(0) # returns the average number of customers in steady state in a system with lambda=1, mu1=1, mu2=2, p1=0.4
#' invlap(Vectorize(TransientMJM2), 0, 5000000, 1000) # returns the transient performance of the model
TransientMJM2=function(s,lambda=1,mu1=1,mu2=2,p1=0.4){# phases are (1,1), (1,2) and (2,*) - overall 3 phases
  p2=1-p1
  A1=diag(lambda,3)
  A0=diag(c(-lambda-2*mu1,-lambda-mu1,-lambda-mu2)-s)
  Am1=rbind(
    c(2*mu1*p1, 2*mu1*p2, 0),
    c(0, 0, mu1),
    c(mu2*p1*p1, mu2*p1*p2, mu2*p2)
  )
  if(lambda*(p1^2/(2*mu1)+p1*p2/mu1+p2/mu2)>1) stop("UNSTABLE")
  A00=-lambda-s
  A01=matrix(lambda*c(p1,p2), ncol=2 )
  A10=matrix(c(mu1,mu2),nrow=2)
  A11=diag(c(-lambda-mu1,-lambda-mu2)-s)
  A12=lambda*rbind(
    c(p1,p2,0),
    c(0,0,1)
  )
  A21=rbind(
    c(2*mu1,0),
    c(0,mu1),
    c(mu2*p1,mu2*p2)
  )
  A22=A0
  Bsh=rbind( # shifted B matrix after column operations
    c(0,0,lambda,0,0,0),
    c(0,0,0,lambda,0,0),
    c(0,0,0,0,lambda,0),
    c(0,0,lambda,0,0,0),
    c(-lambda-mu1-s,0,0,lambda,0,0),
    c(0,-lambda-mu2-s,0,0,lambda,0)
  )
  Csh=rbind( # shifted C matrix after column operations
    c(2*mu1*p2,  0,      -lambda-2*mu1-s, 0,             0,             0),
    c(0,         mu1,    0,               -lambda-mu1-s, 0,             0),
    c(mu2*p1*p2, mu2*p2, 0,               0,             -lambda-mu2-s, 0),
    c(2*mu1*p2,  0,      0,               0,             0,             -lambda-2*mu1-s),
    c(0,         mu1,    0,               0,             0,             p1*(lambda+mu1+s)/p2),
    c(mu2*p1*p2, mu2*p2, 0,               0,             0,             0)
  )
  # Cinv=rbind(
  #   c(-d3,p*l1/(m2*d2),1/m2+d3,-d4),
  #   c(m1/(beta*d1),0,-m1/(beta*d1),1/beta),
  #   c(-1/d1,0,1/d1,0),
  #   c(0,-1/d2,0,1/d2)
  # )
  W= - Bsh %*% solve(Csh)
  # W=rbind(
  #    c(l2/d1, 0, -l2/d1, 0),
  #    c(0, l2/d2, 0, -l2/d2),
  #    c(l2/d1-d1*d3, p*l1*d1/(m2*d2), d1/m2+d1*d3-l2/d1, -d1*d4),
  #    c(m1*d2/(beta*d1), l2/d2, -m1*d2/(beta*d1), d2/beta-l2/d2)
  # )
  #if(max(abs(W+B%*%solve(C)))>1e-8) print(max(abs(W+B%*%solve(C))))

  B0=rbind(cbind(A00,A01,0*A01,matrix(0,nrow=1,ncol=ncol(A12))),
           cbind(A10,A11,0*A11,A12),
           cbind(A10,A11-diag(nrow(A11))*A11,diag(nrow(A11))*A11,A12),
           cbind(matrix(0,nrow=nrow(A21),ncol=1),0*A21,A21,A22),
           cbind(matrix(0,nrow=nrow(A21),ncol=1),0*A21,A21,A22-diag(nrow(A22))*A22))
  B1=rbind(cbind(matrix(0,nrow=1,ncol=ncol(A12)),matrix(0,nrow=1,ncol=ncol(A1))),
           cbind(matrix(0,nrow=nrow(A10),ncol=ncol(A12)),matrix(0,nrow=nrow(A10),ncol=ncol(A1))),
           cbind(matrix(0,nrow=nrow(A10),ncol=ncol(A12)),matrix(0,nrow=nrow(A10),ncol=ncol(A1))),
           Bsh)
  B0=cbind(B0,
           c(rep(0,8),-lambda-2*mu1-s,p1*(lambda+mu1+s)/p2,0))
  valvec=eigen(W)
  xi=valvec$values
  r=valvec$vectors
  q=solve(r)
  #
  W=W-xi[1]*(r[,1]%o%q[1,])-xi[2]*(r[,2]%o%q[2,])
  Lambda=diag(xi*(abs(xi)<0.9999))
  #
  if(s==0)
    pi0=c(rep(0,nrow(B0)+nrow(Csh)-1),1) %*% solve(
      cbind(
        rbind(cbind(B0,B1),
              cbind(matrix(0,ncol=ncol(B0),nrow=nrow(Csh)),Csh)),
        c(rep(0,nrow(B0)),r[,1]),
        c(rep(1,nrow(B0)), solve(diag(nrow(W))-W) %*% rep(1,nrow(W)))
      ))
  else
    pi0=c(-1,rep(0,nrow(B0)+nrow(Csh)-1)) %*% solve(
      cbind(
        rbind(cbind(B0,B1),
              cbind(matrix(0,ncol=ncol(B0),nrow=nrow(Csh)),Csh)),
        c(rep(0,nrow(B0)),r[,1]),
        c(rep(0,nrow(B0)),r[,2])
      )
    )
  pi0=as.vector(Re(pi0))
  pi1=pi0[2:5]
  pi2=pi0[6:11]
  pi3=pi0[12:17]
  I=diag(nrow(Lambda))
  stat=sum(pi1)+ 2*sum(pi2)+ sum(pi3 %*% r %*% (3*I-2* Lambda) %*% solve(I-Lambda) %*% solve(I-Lambda) %*% q)
  stat=sum(pi1)+ 2*sum(pi2)+ sum(pi3 %*%(3*I-2* W) %*% solve(I-W) %*% solve(I-W))
  return(stat)
}
