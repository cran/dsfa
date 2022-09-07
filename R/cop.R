#' Copula distribution
#'
#' Probablity density function and random number generation for the normal, frank and gumbel bivariate copula.
#'
#' @param U matrix of pseudo observations. Must have two columns.
#' @param Tau matrix of Kendall's tau.
#' @param family integer, defines the copula family:\cr
#' `1` = Gaussian copula \cr
#' @param deriv derivative of order \code{deriv} of the log density. Available are 1,2,3,4.
#' @param disjoint logical; if TRUE, only derivatives with respect to Tau are provided.
#' @param num logical; if TRUE, numerical derivatives are provided.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#'
#' @details For more than 2 dimensions only the gaussian copula is implemented. The functions \code{pcop} and \code{rcop} are wrapper functions for the \code{VineCopula} package.
#'
#' @return \code{dcopula} gives the density, \code{pcop} gives the distribution function for a specified copula and \code{rcop} generates random numbers, with given Tau.
#' If the derivatives are calculated these are provided as the attributes \code{gradient}, \code{hessian}, \code{l3} and \code{l4} of the output of the density.
#'
#' @examples
#' pdf <- dcop(U=matrix(c(0.3,0.7), ncol=2), Tau=matrix(0.5,ncol=1), family=1)
#' cdf <- pcop(U=matrix(c(0.3,0.7), ncol=2), Tau=matrix(0.5,ncol=1), family=1)
#' r <- rcop(n=100, Tau=matrix(0.5,nrow=100), family=1)
#'
#' @references
#' \itemize{
#' \item \insertRef{schepsmeier2014derivatives}{dsfa}
#' \item \insertRef{hofert2018elements}{dsfa}
#' }
#' @export
dcop<-function(U, Tau=0, family=1, deriv = 0, disjoint=TRUE, num=FALSE, log.p = FALSE){
  #density of copula
  if(!is.matrix(U)){
    stop(paste("U must be a matrix", "\n", ""))
  }

  #Correct for numerical inaccuracy
  U[U>1-1e-16]<-1-1e-16
  U[U<0+1e-16]<-0+1e-16

  #Get parameters
  N<-nrow(U)
  D<-ncol(U)

  #Check for correct inputs
  if (!is.numeric(as.numeric(family))){
    stop(paste("family must be an integer", "\n", ""))
  }


  if(!is.matrix(Tau)){
    if(length(Tau)==N){
      Tau<-matrix(Tau, nrow=N, ncol=1)
    } else {
      Tau<-matrix(Tau, nrow=N, ncol=length(Tau), byrow=T)
    }
  }

  if(any(Tau<(-1))|any(Tau>1)){
    stop(paste("Tau must be in [-1,1]", "\n", ""))
  }

  #Check if the correct family is chosen for the dimensions
  if(D>2 & family!=1){
    stop(paste("For more than 3 dimensions only the normal copula is implemented", "\n", ""))
  }

  if(D>2 & disjoint & !num){
    stop(paste("For more than 3 dimensions only numerical disjoint derivatives are implemented", "\n", ""))
  }

  #Initialize out
  out<-matrix(0, ncol=1, nrow=N)

  if(num){
    out<-sapply(1:N, function(n) dcop_copula(U[n,], Tau[n,], family=1, log.p = TRUE))

    if(deriv>0){
      if(disjoint){
        #Gradient
        l1_disjoint<-lapply(1:N, function(n) numDeriv::jacobian(dcop_copula_deriv, x=Tau[n,], U=U[n,], family=family))
        l1<-lapply(1:N, function(n) c(rep(0,D),l1_disjoint[[n]]))

        #Hessian
        l2_disjoint_index<-lower.tri(diag(ncol(Tau)), diag=T)
        l2_disjoint<-lapply(1:N, function(n) numDeriv::hessian(dcop_copula_deriv, x=Tau[n,], U=U[n,], family=family)[l2_disjoint_index])

        DD<-max(mgcv::trind.generator(length(c(U[1, , drop=F],Tau[1, , drop=F])))$i2[,c(1:D)])
        l2<-lapply(1:N, function(n) c(rep(0,DD),l2_disjoint[[n]]))
      } else {
        #Gradient
        l1<-lapply(1:N, function(n) numDeriv::jacobian(dcop_copula_deriv, x=cbind(U[n, , drop=F],Tau[n, , drop=F]), family=family, D=D))

        #Hessian
        l2_index<-lower.tri(diag(length(c(U[1, , drop=F],Tau[1, , drop=F]))), diag=T)
        l2<-lapply(1:N, function(n) numDeriv::hessian(dcop_copula_deriv, x=cbind(U[n, , drop=F],Tau[n, , drop=F]), family=family, D=D)[l2_index])
      }

      #Transfrom l1 and l2 from lists into matrices
      l1<-matrix(unlist(l1), byrow=T, nrow=N)
      l2<-matrix(unlist(l2), byrow=T, nrow=N)
    }
  } else {
    if(disjoint){
      u<-U[,1]
      v<-U[,2]

      out<-NULL

      if(family==1){
        #Gaussian copula

        #Compute Tau
        p<-sin(Tau*pi/2)
        out<-(p * (p * erfcinv(2 *u)^2 - 2 * erfcinv(2 * u) *erfcinv(2 *v) +
                     p * erfcinv(2 * v)^2))/(-1 + p^2) - 1/2 *log(1 - p^2)

        if(deriv>0){
          d1taudp1<-1/2*pi*cos((pi* Tau)/2)

          f1<-(p - p^3 + 2 * (p *erfcinv(2 * u) - erfcinv(2 *v)) * (-erfcinv(2 *u) + p *erfcinv(2 * v)))/(-1 + p^2)^2

          d2taudp2<--(1/4) * pi^2 * sin((pi * Tau)/2)

          f2<-(-1 + p^4 + (2 + 6 * p^2) * erfcinv(2 * u)^2 -
                 4 * p * (3 + p^2) * erfcinv(2 *u) * erfcinv(
                   2 * v) + (2 + 6 * p^2)*erfcinv(2 * v)^2)/(-1 + p^2)^3
          if(deriv>2){
            d3taudp3<--(1/8) * pi^3 * cos((pi * Tau)/2)

            f3<-(-2 * p * (-3 + 2 *p^2 + p^4) +
                   12 * ((1 + p^2) * erfcinv(2 * u) -
                           2 * p * erfcinv(2 *v)) * (-2 * p * erfcinv(
                             2 * u) + (1 + p^2) * erfcinv(2 * v)))/(-1 + p^2)^4

            d4taudp4<-1/16 * pi^4 * sin((pi * Tau)/2)

            f4<-(6 * (-1 - 5 * p^2 + 5 * p^4 + p^6 +
                        4 * (1 + 5 * p^2 * (2 + p^2)) * erfcinv(2 * u)^2 -
                        8 * p * (5 + 10 * p^2 + p^4) * erfcinv(2 * u) * erfcinv(2 * v) +
                        4 * (1 + 5 * p^2 * (2 + p^2)) * erfcinv(2 * v)^2))/(-1 + p^2)^5
          } else {
            g3<-NULL
            g4<-NULL
            d3taudp3<-NULL
            d4taudp4<-NULL
          }
        }
      }

      if(family==3){
        #Clayton copula

        #Compute Tau
        p<-2*(Tau/(1 - Tau))
        out<-log((1 + p) * (u * v)^(-1 - p) * (-1 + u^-p + v^-p)^(-2 - 1/p))

        if(deriv>0){
          d1taudp1<-2/(-1 + Tau)^2

          f1<--((p * (1 + p) * (1 + 2 * p) * v^p * log(u) + p * (1 + p) * (1 + 2 * p) * u^ p * log(v) + (-u^p + (-1 + u^p) * v^p) * (p^2 * (-1 + (1 + p) * log(u * v)) - (1 + p) * log(-1 + u^-p +
                                                                                                                                                                                         v^-p)))/(p^2 * (1 + p) * (-v^p + u^p * (-1 + v^p))))

          d2taudp2<--(4/(-1 + Tau)^3)

          f2<-((1 + 2 * p) * (p + p^2)^2 * u^p * v^p * (-1 + v^p) * log(u)^2 + 2 * p * (1 + p)^2 * u^p * (-v^p + u^p * (-1 + v^p)) * log(v) + (1 + 2 * p) * (p + p^2)^2 * u^
                 p * (-1 + u^p) * v^p * log(v)^2 + 2 * p * (1 + p)^2 * v^p * log(u)* (-v^p + u^p * (-1 + v^p + p * (1 + 2 * p) * log(v))) - (v^p - u^p * (-1 + v^p))^2 * (p^3 +
                                                                                                                                                                            2 * (1 + p)^2 * log(-1 + u^-p + v^-p)))/(p^3 * (1 + p)^2 * (v^p - u^p * (-1 + v^p))^2)
          if(deriv>2){
            d3taudp3<-12/(-1 + Tau)^4

            f3<--((2 * p^4 * (v^p - u^p * (-1 + v^p))^3 + (1 +p)^3 * (p^3 * (1 + 2 * p) * u^p * v^p * (-1 + v^p) * (v^p + u^p * (-1 + v^p)) * log(u)^3 +
                                                                        6 * p * u^p * (v^p - u^p * (-1 + v^p))^2 * log(v) + 3 * p^2 * u^p * (-1 + u^p) * v^p * (-v^p + u^p * (-1 + v^p)) * log(v)^2 +
                                                                        p^3 * (1 + 2 * p) * u^p * (-1 + u^p) * v^p * (u^p + (-1 + u^p) * v^p) * log(v)^3 + 3 * p^2 * u^p * v^p * log(u)^2 * ((-1 + v^p) * (-v^p + u^p * (-1 + v^p)) +
                                                                                                                                                                                               p * (1 + 2 * p) * (v^p + u^p * (-1 + v^p)) * log(v)) +  3 * p * v^p * log(u) * (2 * (v^p - u^p * (-1 + v^p))^2 +
                                                                                                                                                                                                                                                                                 p * u^p * log(v) * (u^p * (-2 + p * (1 + 2 * p) * log(v)) + (-1 + u^p) * v^p * (2 + p * (1 + 2 * p) * log(v)))) +
                                                                        6 * (v^p - u^p * (-1 + v^p))^3 * log(-1 + u^-p + v^-p)))/(p^4 * (1 + p)^3 * (-v^p + u^p * (-1 + v^p))^3))

            d4taudp4<--(48/(-1 + Tau)^5)

            f4<--(6/(1 + p)^4) + (p^4 * (1 + 2 *p) * u^p * v^p * (-1 + v^p) * (v^(2 * p) + 4 * u^p * v^p * (-1 + v^p) +
                 u^(2 * p) * (-1 + v^p)^2) * log(u)^4 - 24 * p * u^p * (u^p + v^p - u^p * v^p)^3 * log(v) +
                 12 * p^2 * u^p * (-1 + u^p) * v^p * (u^p + v^p - u^p * v^p)^2 * log(v)^2 + 4 * p^3 * u^p * (-1 + u^p) * v^p * (-u^(2 * p) + (-1 + u^p)^2 * v^(2 * p)) * log(v)^3 + p^4 * (1 + 2 * p) * u^p * (-1 + u^p) * v^
                  p * (u^(2 * p) + 4 * u^p * (-1 + u^p) * v^p + (-1 + u^p)^2 * v^(2 * p)) * log(v)^4 + 4 * p^3 * u^p * v^p * log(u)^3 * (-v^(2 * p) * (-1 + v^p) + u^(2 * p) * (-1 + v^p)^3 +  p * (1 + 2*  p) * (v^(2 * p) + 4 * u^p * v^p * (-1 + v^p) +
                  u^(2 * p) * (-1 + v^p)^2) * log(v)) + 6 * p^2 * u^p * v^p * log(u)^2 * (2 * (-1 + v^p) * (v^p - u^p * (-1 + v^p))^2 +
                  p * log(v)* (-2 * v^(2 * p) + 2 * u^(2 * p) * (-1 + v^p)^2 + p * (1 + 2 * p) * (4 * u^p * v^p - v^(2 * p) +
                                                                                                                                                                                                                                                                                                                                                                                                  u^(2 * p) * (-1 + v^(2 * p))) * log(v))) + 4 * p * v^p * log(u) * (-6 * (v^p - u^p * (-1 + v^p))^3 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       p * u^p * log(v) * (4 * u^p * (-1 + u^p) * v^p * (-3 + p^2 * (1 + 2 * p) * log(v)^2) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             u^(2 * p) * (6 + p * log(v) * (-3 + p * (1 + 2 * p) * log(v))) + (-1 + u^p)^2 * v^(2 * p) * (6 + p * log(v)* (3 + p * (1 + 2 * p) * log(v))))) -
                 24 * (u^p + v^p - u^p * v^p)^4 * log(-1 + u^-p + v^-p))/( p^5 * (v^p - u^p * (-1 + v^p))^4)

          } else {
            g3<-NULL
            g4<-NULL
            d3taudp3<-NULL
            d4taudp4<-NULL
          }
        }
      }

      if(deriv>0){
        l<-chainrule(f1=f1, f2=f2, f3=f3, f4=f4,
                     g1=as.matrix(d1taudp1), g2=as.matrix(d2taudp2), g3=as.matrix(d3taudp3), g4=as.matrix(d4taudp4),
                     deriv=deriv, xg=NULL)

        l1<-matrix(cbind(0,0,l$h1), nrow=N)
        l2<-matrix(cbind(0,0,0,0,0,l$h2), nrow=N)

        if(deriv>2){
          attr(out,"l3")<-l$h3
          attr(out,"l4")<-l$h4
        }
      }
    } else {
      if(family==1){
        #Quantile transformation
        X<-stats::qnorm(U)

        #Get parameters
        npar_tau<-ncol(Tau)
        npar_all<-D+npar_tau
        ID<-diag(D)

        xg_all<-mgcv::trind.generator(npar_all)
        xg_x<-mgcv::trind.generator(D)

        # xg_P<-trind.generator(npar_tau)
        # diag(xg_P$i2)<-0

        xN<-matrix(0,D,1)
        PN<-matrix(0,D,D)

        xg_P<-list()
        xg_P$i1<-PN
        xg_P$i1[upper.tri(ID)]<-1:npar_tau
        xg_P$i1<-xg_P$i1+t(xg_P$i1)-diag(diag(xg_P$i1))

        #Initialize out, l1 and l2
        out<-rep(0,N)
        l1<-matrix(0, nrow=N, ncol=D+npar_tau)
        l2<-matrix(0, nrow=N, ncol=max(xg_all$i2))

        for(n in 1:N){
          # print(n)

          #Get pseudo observations n
          x<-t(X[n, , drop=FALSE])
          Z<-x%*%t(x)

          #Create correlation matrix P
          P<-ID
          P[upper.tri(P)]<-sin(Tau[n,]*pi/2)#Tau[n,]
          P<-P+t(P)-diag(diag(P))
          P<-pos_eigen(P)$I
          P_inv<-solve(P)

          #Numerical convience
          h0<-(P_inv-ID)
          h4<-h0%*%P_inv

          #Calculate loglike
          out[n]<--1/2*log(det(P))-1/2*sum(diag(h0%*%Z))#-1/2*log(det(P))-1/2*sum(diag((P_inv-ID)%*%Z))

          if(deriv>0){
              #Derivatives of qnorm wrt U
              d1qnormdu1<-d1qnormdx1(U[n, , drop=FALSE])
              d2qnormdu2<-diag(c(d2qnormdx2(U[n, , drop=FALSE])))

              #Derivatives of Tau wrt p
              d1taudp1<-1/2*pi*cos((pi* Tau[n,])/2)
              d2taudp2<--(1/4) * pi^2 * sin((pi * Tau[n,, drop=T])/2)

              if(length(d2taudp2)==1){
                d2taudp2<-as.matrix(d2taudp2, nrow=1)
              } else {
                d2taudp2<-diag(d2taudp2)
              }

              h9<-P_inv%*%Z
              #Deriv of lcop wrt U
              for(i in 1:npar_all){
                if(i<=D){
                  #First deriv wrt U
                  d1xdui1<-xN
                  d1xdui1[i]<-1
                  h1<-(d1xdui1*d1qnormdu1[,i])%*%t(x)
                  l1[n,i]<--1/2*sum(diag(h0%*%(h1+t(h1))))

                  #Second deriv wrt U and Tau
                  for(j in i:npar_all){
                    if(j<=D){
                      d1xdxi1<-xN
                      d1xdxi1[i]<-1
                      d1xidui1<-d1qnormdu1[,i]

                      d1xduj1<-xN
                      d1xduj1[j]<-1
                      d1xjduj1<-d1qnormdu1[,j]

                      h6<-d1xdui1%*%t(d1xduj1)
                      l2[n,xg_all$i2[i,j]]<--1/2*sum(diag(h0%*%(h6+t(h6))))#-1/2*sum(diag((P_inv-ID)%*%((d1xdui1%*%t(d1xduj1))+(d1xduj1%*%t(d1xdui1)))))

                      d1xdui1%*%t(x)+t(d1xdui1%*%t(x))
                      h7<-d1xdui1%*%t(x)
                      h7+t(h7)
                      l2[n,xg_all$i2[i,j]]<-l2[n,xg_all$i2[i,j]]*d1qnormdu1[,j]*d1qnormdu1[,i]-1/2*sum(diag(h0%*%(h7+t(h7))))*d2qnormdu2[i,j]
                      # l2[n,xg_all$i2[i,j]]<-l2[n,xg_all$i2[i,j]]*d1qnormdu1[,j]*d1qnormdu1[,i]-1/2*sum(diag(h0%*%(d1xdui1%*%t(x)+t(d1xdui1%*%t(x)))))*d2qnormdu2[i,j]
                    } else {
                      d1Pdp1<-PN
                      d1Pdp1[xg_P$i1==(j-D)]<-1
                      l2[n,xg_all$i2[i,j]]<-1/2*sum(diag(P_inv%*%d1Pdp1%*%P_inv%*%(h1+t(h1))))

                      l2[n,xg_all$i2[i,j]]<-l2[n,xg_all$i2[i,j]]*d1taudp1[j-D]
                      # l2[n,xg_all$i2[i,j]]*d1qnormdu1[n,i]
                    }
                  }
                } else {
                  #First deriv wrt Tau
                  d1Pdp1<-PN
                  d1Pdp1[xg_P$i1==(i-D)]<-1
                  h5<-P_inv%*%d1Pdp1

                  l1[n,i]<--1/2*sum(diag(h5%*%(ID-h9)))#-1/2*sum(diag(P_inv%*%d1Pdp1%*%(ID-P_inv%*%Z)))
                  l1[n,i]<-l1[n,i]*d1taudp1[i-D]

                  #Second deriv wrt Tau
                  for(j in i:npar_all){
                    d1Pdq1<-PN
                    d1Pdq1[xg_P$i1==(j-D)]<-1
                    h8<-P_inv%*%d1Pdq1
                    l2[n,xg_all$i2[i,j]]<--1/2*sum(diag(-h8%*%h5%*%(diag(D)-h9)+
                                                          h5%*%h8%*%h9))
                    #-1/2*sum(diag(-P_inv%*%d1Pdq1%*%P_inv%*%d1Pdp1%*%(diag(D)-P_inv%*%Z)+
                    #              P_inv%*%d1Pdp1%*%P_inv%*%d1Pdq1%*%P_inv%*%Z))

                    l2[n,xg_all$i2[i,j]]<-l2[n,xg_all$i2[i,j]]*d1taudp1[j-D]*d1taudp1[i-D]+(-1/2*sum(diag(h5%*%(ID-h9))))*d2taudp2[i-D, j-D]#l2[n,xg_all$i2[i,j]]*d1taudp1[j-D]*d1taudp1[i-D]+(-1/2*sum(diag(h5%*%(ID-P_inv%*%Z))))*d2taudp2[i-D, j-D]
                  }
                }
            }
          }
        }
      }
    }
  }
  if(deriv>0){
    #Set gradient and hessian as attributes
    attr(out,"gradient")<-l1
    attr(out,"hessian")<-l2
  }


  if(!log.p){
   out<-exp(out)
  }

  return(out)
}

#' @describeIn dcop distribution function for the joint distribution.
#' @export
pcop<-function(U, Tau=0, family=1, log.p = FALSE){
  #Density of bivariate copula

  if(!is.matrix(U)){
    stop(paste("U must be a matrix", "\n", ""))
  }

  #Correct for numerical inaccuracy
  U[U>1-1e-16]<-1-1e-16
  U[U<0+1e-16]<-0+1e-16
  N<-nrow(U)

  #Check for correct inputs
  if (!is.numeric(as.numeric(family))){
    stop(paste("family must be an integer", "\n", ""))
  }

  if(!is.matrix(Tau)){
    if(length(Tau)==N){
      Tau<-matrix(Tau, nrow=N, ncol=1)
    } else {
      Tau<-matrix(Tau, nrow=N, ncol=length(Tau), byrow=T)
    }
  }

  if(any(Tau<(-1))|any(Tau>1)){
    stop(paste("Tau must be in [-1,1]", "\n", ""))
  }

  out<-sapply(1:N, function(n) pcop_copula(U[n,], Tau[n,], family=1, log.p=log.p))

  #Return out
  return(out)
}

#' @describeIn dcop random number generation for the joint distribution.
#' @param n number of observations.
#' @export
rcop<-function(n, Tau=0, family=1){
  #Random number generation function of copula

  #Check for correct inputs
  if (!is.numeric(as.numeric(family))){
    stop(paste("family must be an integer", "\n", ""))
  }

  if(!is.matrix(Tau)){
    if(length(Tau)==n){
      Tau<-matrix(Tau, nrow=n, ncol=1)
    } else {
      Tau<-matrix(Tau, nrow=n, ncol=length(Tau), byrow=T)
    }
  }

  if(any(Tau<(-1))|any(Tau>1)){
    stop(paste("Tau must be in [-1,1]", "\n", ""))
  }

  N<-n
  out<-lapply(1:N, function(i) rcop_copula(Tau[i,], family=family))
  out<-matrix(unlist(out), byrow=T, nrow=N)

  return(out)
}
