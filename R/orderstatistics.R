
#' Moments of order statistics for three-parameter Weibull
#'
#' WeibullOrder is the moments of order statistics from three-parameter Weibull distribution.
#'
#' @param q The order of statistics (q-th out of r)
#' @param r The number of replicates to extract the q-th statistics from (q-th out of r)
#' @param shape The shape of three-parameter Weibull distribution (1 is for exponential)
#' @param scale The scale of three-parameter Weibull distribution
#' @param location The location of three-parameter Weibull distribution (0 default)
#' @param pow The order of the moment to obtain statistics for (default pow=1 stands for expectation)
#' @return The pow-th moment of q-th order statistics of r iid. copies of the
#'    r.v. distribued as the three-parameter Weibull with distribution function
#'    \deqn{F(x)=1-exp\left(-\left( \frac{x-location}{scale}\right)^{shape} \right).}
#' @seealso [supportR::ParetoOrder()] for moments from Pareto distribution.
#' @references C.-D. Lai, Generalized Weibull Distributions. Springer: 2014. doi: 10.1007/978-3-642-39106-4.
#'    J. Lieblein, “On Moments of Order Statistics from the Weibull Distribution,” Ann. Math. Statist., vol. 26, no. 2, pp. 330–333, Jun. 1955, doi: 10.1214/aoms/1177728551.
#' @export
#' @examples
#' WeibullOrder() # should be 1 as the expectation of a standard Exponential
#' WeibullOrder(1,10) # should be 0.1 as the minimum of 10 standard Exponentially distributed r.v.'s
WeibullOrder=function(q=1,r=1,shape=1,scale=1,location=0,pow=1){
  Lieblein=function(pow) q*choose(r,q)*gamma(1+pow/shape)*sum(sapply(0:(q-1),function(j){ ((r-q+1+j)^(-1-pow/shape))*choose(q-1,j)*(-1)^j})) # Order statistics of one-parameter Weibull according to Lieblein 1955 paper
  return(sum(sapply(0:pow,function(i) location^i*scale^(pow-i)*Lieblein(pow-i))))
}

#' Moments of order statistics for two-parameter Pareto
#'
#' ParetoOrder is the moments of order statistics from two-parameter Pareto distribution.
#'
#' @param q The order of statistics (q-th out of r)
#' @param r The number of replicates to extract the q-th statistics from (q-th out of r)
#' @param shape The shape of two-parameter Pareto distribution
#' @param scale The scale of two-parameter Pareto distribution
#' @param pow The order of the moment to obtain statistics for (default pow=1 stands for expectation)
#' @return The pow-th moment of q-th order statistics of r iid. copies of the
#'    r.v. distribued as the two-parameter Pareto with distribution function
#'    \deqn{F(x)=1-\left(\frac{scale}{x}\right)^{shape}, \;\; x\geq scale.}
#' @seealso [supportR::WeibullOrder()] for moments from three-parameter Weibull distribution.
#' @references J. S. Huang, “A Note on Order Statistics from Pareto Distribution,” Scandinavian Actuarial Journal, vol. 1975, no. 3, pp. 187–190, Jul. 1975, doi: 10.1080/03461238.1975.10405095.
#' @export
#' @examples
#' ParetoOrder() # should be 2
ParetoOrder=function(q=1,r=1,scale=1,shape=2,pow=1){
  stopifnot(shape>q/(r-pow+1))
  return((scale^pow)*gamma(r+1)*gamma(r+1-q-pow/shape)/(gamma(r-q+1)*gamma(r+1-pow/shape)))
}

#' Moments of order statistics for Phase-type distribution PH(beta,B)
#'
#' PHOrder calculates the moments of order statistics from Phase-type distribution PH(beta,B).
#'
#' @param q The order of statistics (q-th out of r)
#' @param r The number of replicates to extract the q-th statistics from (q-th out of r)
#' @param B The matrix describing Phase-type distribution
#' @param beta The initial probability vector describing Phase-type distribution
#' @param pow The order of the moment to obtain statistics for (default pow=1) - CURRENTLY UNUSED
#' @return The pow-th moment of q-th order statistics of r iid. copies of the
#'    r.v. distribued as the Phase-type distribution PH\eqn{(\beta,B)} with density
#'    \deqn{F(x)=1-\beta e^{B x}\mathbf{e},}
#'    where \eqn{\mathbf{e}=(1,\dots,1)} is vector of ones.
#' @seealso [supportR::WeibullOrder()] for moments from three-parameter Weibull distribution.
#' @references M. Bladt and B. F. Nielsen, Matrix-Exponential Distributions in Applied Probability, vol. 81. in Probability Theory and Stochastic Modelling, vol. 81. Boston, MA: Springer US, 2017. doi: 10.1007/978-1-4939-7049-0.
#' @export
#' @examples
#' PHOrder() # should be 1 as the expectation of standard Exponential
#' PHOrder(1,10,-1,1) # should be 0.1 as the minimum of 10 standard Exponentially distributed r.v.'s
PHOrder=function(q=1,r=1,B=-1,beta=1,pow=1){
  B=as.matrix(B)
  m=nrow(B)
  b0=-apply(B,1,sum)
  "%p%" = function(A,B){A %x% diag(nrow=nrow(B))+diag(nrow=nrow(A)) %x% B} # Kronecker plus
  B0=function(k) # generate the departure intensity from level from k to k-1
  {
    tmp=sapply(1:k,function(i){diag(m^(i-1)) %x% b0 %x% diag(m^(k-i))},simplify="array")
    if(is.vector(tmp))
      return(matrix(sum(tmp)))
    else
      return(apply(tmp,c(1,2),sum))
  }
  kron.plus.power=function(B,k) # oplus iterated k times
  {
    T=matrix(0)
    for(i in 1:k) T=T %p% B
    T
  }
  kron.mult.power=function(B,k) # oplus iterated k times
  {
    T=matrix(1)
    for(i in 1:k) T=T %x% B
    T
  }
  diagjoin=function(A,B){
    C=matrix(0,ncol=ncol(A)+ncol(B),nrow=nrow(A)+nrow(B))
    C[1:nrow(A),1:ncol(A)]=A
    C[nrow(A)+(1:nrow(B)),ncol(A)+(1:ncol(B))]=B
    return(C)
  }
  Bqr=function(q,r){
    T=kron.plus.power(B,r)
    if(q>1){
      for(i in 2:q)
        T=diagjoin(T,kron.plus.power(B,r-i+1))
      B0r=B0(r)
      if(q>2)
        for(i in 2:(q-1))
          B0r=diagjoin(B0r,B0(r-i+1))
      B0r=cbind(matrix(0,ncol=ncol(kron.plus.power(B,r)),nrow=nrow(B0r)),B0r)
      B0r=rbind(B0r,matrix(0,ncol=ncol(B0r),nrow=nrow(T)-nrow(B0r)))
      T=T+B0r
    }

    return(T)
  }
  betaqr=function(q,r)
    c(kron.mult.power(beta,r),rep(0,ifelse(q>1,sum(m^((r-q+1):(r-1))),0)))
  return(as.vector(-betaqr(q,r) %*% solve(Bqr(q,r)) %*% rep(1,length(betaqr(q,r)))))
}
