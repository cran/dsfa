#' Copula function
#'
#' Probablitiy density function, distribution and random number generation for copulas.
#'
#' @param W matrix of pseudo observations. Must have at least two columns.
#' @param delta matrix of copula parameter. Must have at least one column.
#' @param family_cop string, defines the copula family:\cr
#' `independent` = Independence copula \cr
#' `normal` = Gaussian copula \cr
#' `clayton` = Clayton copula \cr
#' `gumbel` = Gumbel copula \cr
#' `frank` = Frank copula \cr
#' `joe` = Joe copula \cr
#' @param deriv derivative of order \code{deriv} of the log density. Available are \code{0},\code{2}.
#' @param log.p logical; if \code{TRUE}, probabilities p are given as \code{log(p)}.
#'
#' @details For more than 2 dimensions only the gaussian copula is implemented. The functions \code{pcop} and \code{rcop} are wrapper functions for 'copula' package.
#' The functions \code{pcop} and \code{rcop} are wrapper functions for the \code{\link[copula:pCopula]{pCopula()}} and \code{\link[copula:rCopula]{rCopula()}}.
#' Although the parameter space is larger in theory for some copulas, numeric under- and overflow limit the parameter space. The intervalls for the parameter \code{delta} are given as follows:
#' \enumerate{
#'  \item  `independent`, min=0 and max=1
#' \item  `normal`, min=-1 and max=1
#' \item  `clayton`, min=1e-16 and max=28
#' \item  `gumbel`, min=1 and max=17
#' \item  `frank`, min=-35 and max=35
#' \item  `joe`, min=1e-16 and max=30
#' }
#'
#' @return \code{dcop} gives the density, \code{pcop} gives the distribution function for a specified copula and \code{rcop} generates random numbers, with given delta.
#' If the derivatives are calculated these are provided as the attributes \code{gradient}, \code{hessian} of the output of the density.
#'
#' @examples
#' u=0.3; v=0.7; p=0.5
#' pdf <- dcop(W=matrix(c(u,v), ncol=2), delta=matrix(p,ncol=1), family_cop="normal")
#' cdf <- pcop(W=matrix(c(u,v), ncol=2), delta=matrix(p,ncol=1), family_cop="normal")
#' r <- rcop(n=100, delta=matrix(p,nrow=100), family_cop="normal")
#'
#' @references
#' \itemize{
#' \item \insertRef{schepsmeier2014derivatives}{dsfa}
#' \item \insertRef{hofert2018elements}{dsfa}
#' }
#' @export
dcop<-function(W, delta, family_cop="normal", log.p=FALSE, deriv=0){
  #density of copula

  #Check for correct inputs
  if (!family_cop%in%c("independent","normal","clayton","gumbel","frank","joe")){
    stop(paste("Incorrect family_cop", "\n", ""))
  }

  if(!is.matrix(W)){
    stop(paste("W must be a matrix", "\n", ""))
  }

  #Correct for numerical inaccuracy
  W[W>1-1e-16]<-1-1e-16
  W[W<0+1e-16]<-0+1e-16

  if(!is.matrix(delta)){
    stop(paste("delta must be a matrix", "\n", ""))
  }
  #Get parameters
  N<-nrow(W)
  D<-ncol(W)
  p<-delta

  #Check if the dimension is 2
  if(D>2 & family_cop!="normal"){
    stop(paste("For more than two dimensions only the normal copula is implemented", "\n", ""))
  } else {
    u<-W[,1]
    v<-W[,2]
  }

  #Initialize value
  value<-NULL

  if(family_cop=="independent"){
    value <- rep(0,N)

    if(deriv>0){
      gradient <- matrix(0, nrow=N, ncol=3)
      hessian <- matrix(0, nrow=N, ncol=6)
    }
  }


  if(family_cop == "normal") {
    #Check for correct inputs
    if(any(-1>=p|p>=1)){ #p in (-1,1)
      stop(paste("p must be in (-1,1)", "\n", ""))
    }

    #Quantile transformation
    X<-stats::qnorm(W)

      #Get parameters
      npar_delta<-ncol(delta)
      npar_all<-D+npar_delta
      ID<-diag(D)

      tri_all<-mgcv::trind.generator(npar_all)
      tri_x<-mgcv::trind.generator(D)

      # tri_P<-trind.generator(npar_delta)
      # diag(tri_P$i2)<-0

      xN<-matrix(0,D,1)
      PN<-matrix(0,D,D)

      tri_P<-list()
      tri_P$i1<-PN
      tri_P$i1[upper.tri(ID)]<-1:npar_delta
      tri_P$i1<-tri_P$i1+t(tri_P$i1)-diag(diag(tri_P$i1))

      #Initialize out, l1 and l2
      out<-rep(0,N)
      l1<-matrix(0, nrow=N, ncol=D+npar_delta)
      l2<-matrix(0, nrow=N, ncol=max(tri_all$i2))

      for(n in 1:N){
        # print(n)

        #Get pseudo observations n
        x<-t(X[n, , drop=FALSE])
        Z<-x%*%t(x)

        #Create correlation matrix P
        P<-ID
        P[upper.tri(P)]<-delta[n,]
        P<-P+t(P)-diag(diag(P))
        P<-pos_eigen(P)$I
        P_inv<-solve(P)

        #Numerical convience
        h0<-(P_inv-ID)
        h4<-h0%*%P_inv

        #Calculate loglike
        out[n]<--1/2*log(det(P))-1/2*sum(diag(h0%*%Z))#-1/2*log(det(P))-1/2*sum(diag((P_inv-ID)%*%Z))

        if(deriv>0){
          #Derivatives of qnorm wrt W
          d1qnormdu1<-d1qnormdx1(W[n, , drop=FALSE])
          d2qnormdu2<-diag(c(d2qnormdx2(W[n, , drop=FALSE])))

          h9<-P_inv%*%Z
          #Deriv of lcop wrt W
          for(i in 1:npar_all){
            if(i<=D){
              #First deriv wrt W
              d1xdui1<-xN
              d1xdui1[i]<-1
              h1<-(d1xdui1*d1qnormdu1[,i])%*%t(x)
              l1[n,i]<--1/2*sum(diag(h0%*%(h1+t(h1))))

              #Second deriv wrt W and delta
              for(j in i:npar_all){
                if(j<=D){
                  d1xdxi1<-xN
                  d1xdxi1[i]<-1
                  d1xidui1<-d1qnormdu1[,i]

                  d1xduj1<-xN
                  d1xduj1[j]<-1
                  d1xjduj1<-d1qnormdu1[,j]

                  h6<-d1xdui1%*%t(d1xduj1)
                  l2[n,tri_all$i2[i,j]]<--1/2*sum(diag(h0%*%(h6+t(h6))))
                  h7<-d1xdui1%*%t(x)
                  l2[n,tri_all$i2[i,j]]<-l2[n,tri_all$i2[i,j]]*d1qnormdu1[,j]*d1qnormdu1[,i]-1/2*sum(diag(h0%*%(h7+t(h7))))*d2qnormdu2[i,j]
                } else {
                  d1Pdp1<-PN
                  d1Pdp1[tri_P$i1==(j-D)]<-1
                  l2[n,tri_all$i2[i,j]]<-1/2*sum(diag(P_inv%*%d1Pdp1%*%P_inv%*%(h1+t(h1))))
                }
              }
            } else {
              #First deriv wrt delta
              d1Pdp1<-PN
              d1Pdp1[tri_P$i1==(i-D)]<-1
              h5<-P_inv%*%d1Pdp1
              l1[n,i]<--1/2*sum(diag(h5%*%(ID-h9)))

              #Second deriv wrt delta
              for(j in i:npar_all){
                d1Pdq1<-PN
                d1Pdq1[tri_P$i1==(j-D)]<-1
                h8<-P_inv%*%d1Pdq1
                l2[n,tri_all$i2[i,j]]<--1/2*sum(diag(-h8%*%h5%*%(diag(D)-h9)+
                                                       h5%*%h8%*%h9))
              }
            }
          }
        }
      }
    # }
    # value <- -1/2*log(1-p^2)+p/(1-p^2)*qnorm(u)*qnorm(v)-p^2/(2*(1-p^2))*(qnorm(u)^2 +qnorm(v)^2)
    #
    # if(deriv>0) {
    #   d1qnormdu1 <- d1qnormdx1(u)
    #   d1qnormdv1 <- d1qnormdx1(v)
    #
    #   d2qnormdu2 <- d2qnormdx2(u)
    #   d2qnormdv2 <- d2qnormdx2(v)
    #
    #   gradient <- matrix(0, nrow = N, ncol = 3)
    #   hessian <- matrix(0, nrow = N, ncol = 6)
    #   gradient[, 1] <- p/(1-p^2)*d1qnormdu1*qnorm(v)-p^2/(1-p)*qnorm(u)*d1qnormdu1
    #   gradient[, 2] <- p/(1-p^2)*qnorm(u)*d1qnormdv1-p^2/(1-p)*qnorm(v)*d1qnormdv1
    #   gradient[, 3] <- 1/(1-p^2)+(p^2+1)/(1-p^2)^2*qnorm(u)*qnorm(v)-p/(1-p^2)^2*(qnorm(u)^2+qnorm(v)^2)
    #   hessian[, 1] <- p/(1-p^2)*d2qnormdu2*qnorm(v)-p^2/(1-p^2)*(d1qnormdu1^2+qnorm(u)*d2qnormdu2)
    #   hessian[, 2] <- p/(1-p^2)*d1qnormdu1*d1qnormdv1
    #   hessian[, 3] <- (p^2+1)/(1-p^2)^2*d1qnormdu1*qnorm(v)-2*p/(1-p^2)^2*qnorm(u)*d1qnormdu1
    #   hessian[, 4] <- p/(1-p^2)*qnorm(u)*d2qnormdv2-p^2/(1-p^2)*(d1qnormdv1^2+qnorm(v)*d2qnormdv2)
    #   hessian[, 5] <- (p^2+1)/(1-p^2)^2*qnorm(u)*d1qnormdv1-2*p/(1-p^2)^2*qnorm(v)*d1qnormdv1
    #   hessian[, 6] <- (p^2+1)/(1-p^2)^2+2*p*(p^2+3)/(1-p^2)^3*qnorm(u)*qnorm(v)-(-3*p-1)/(1-p^2)^3*(qnorm(u)^2+qnorm(v)^2)
    # }
    value<-out
    gradient<-l1
    hessian<-l2
  }

  if(family_cop=="clayton"){
    #Check for correct inputs
    if(any(0>=p|p>28)){ #p in (0,28]
      stop(paste("p must be in (0,28]", "\n", ""))
    }

    expr1 <- 1 + p
    expr3 <- -1
    expr4 <- expr3 - p
    expr5 <- u * v
    expr6 <- log(expr5)
    expr10 <- expr3/p - 2
    expr11 <- -p
    expr12 <- u^expr11
    expr13 <- v^expr11
    expr15 <- expr12 + expr13 - 1
    expr16 <- log(expr15)
    expr19 <- v/expr5
    expr21 <- p + 1
    expr22 <- -expr21
    expr23 <- u^expr22
    expr24 <- expr23 * p
    expr25 <- expr24/expr15
    expr29 <- expr5^2
    expr33 <- -(expr21 + 1)
    expr39 <- expr15^2
    expr50 <- v^expr22
    expr51 <- expr50 * p
    expr56 <- p^2
    expr57 <- 1/expr56
    expr59 <- log(u)
    expr64 <- log(v)
    expr65 <- expr13 * expr64
    expr66 <- expr12 * expr59
    expr67 <- expr65 + expr66
    expr75 <- u/expr5
    expr77 <- expr51/expr15
    expr108 <- expr67/expr15
    expr112 <- expr57 * expr108
    value <- log(expr1) + expr4 * expr6 + expr10 * expr16

    if(deriv>0){
      gradient <- matrix(0, nrow=N, ncol=3)
      hessian <- matrix(0, nrow=N, ncol=6)
      gradient[, 1] <- expr4 * expr19 - expr10 * expr25
      hessian[, 1] <- -(expr4 * (v * v/expr29) - expr10 * (u^expr33 * expr21 * p/expr15 - expr24 *  expr24/expr39))
      hessian[, 2] <-  expr4 * (1/expr5 - v * u/expr29) - expr10 * (expr24 * expr51/expr39)
      hessian[, 3] <- -(expr19 + (expr57 * expr25 + expr10 * ((expr23 - expr23 * expr59 * p)/expr15 + expr24 * expr67/expr39)))
      gradient[, 2] <- expr4 * expr75 - expr10 * expr77
      hessian[, 4] <- -(expr4 * (u * u/expr29) - expr10 * (v^expr33 * expr21 * p/expr15 - expr51 * expr51/expr39))
      hessian[, 5] <- -(expr75 + (expr57 * expr77 + expr10 * ((expr50 - expr50 * expr64 * p)/expr15 + expr51 * expr67/expr39)))
      gradient[, 3] <- 1/expr1 - expr6 + (expr57 * expr16 - expr10 * expr108)
      hessian[, 6] <- -(expr112 + 2 * p/expr56^2 * expr16 + (expr112 - expr10 * ((expr66 * expr59 + expr65 * expr64)/expr15 - expr67 * expr67/expr39)) + 1/expr1^2)
    }
  }

  if(family_cop=="gumbel"){
    #Check for correct inputs
    if(any(1>p|p>17)){ #p in [1,17]
      stop(paste("p must be in [1,17]", "\n", ""))
    }

    expr1 <- log(u)
    expr2 <- -expr1
    expr3 <- expr2^p
    expr4 <- log(v)
    expr5 <- -expr4
    expr6 <- expr5^p
    expr7 <- expr3 + expr6
    expr8 <- 1/p
    expr9 <- expr7^expr8
    expr11 <- u * v
    expr16 <- -2 + 2/p
    expr17 <- log(expr7)
    expr20 <- p - 1
    expr21 <- expr1 * expr4
    expr22 <- log(expr21)
    expr26 <- -1/p
    expr27 <- expr7^expr26
    expr29 <- 1 + expr20 * expr27
    expr32 <- expr8 - 1
    expr33 <- expr7^expr32
    expr34 <- expr2^expr20
    expr35 <- 1/u
    expr36 <- p * expr35
    expr37 <- expr34 * expr36
    expr38 <- expr8 * expr37
    expr42 <- expr37/expr7
    expr45 <- expr35 * expr4
    expr46 <- expr45/expr21
    expr49 <- expr26 - 1
    expr50 <- expr7^expr49
    expr51 <- expr26 * expr37
    expr52 <- expr50 * expr51
    expr53 <- expr20 * expr52
    expr57 <- 1/u^2
    expr61 <- expr21^2
    expr67 <- expr20 - 1
    expr72 <- expr34 * (p * expr57) + expr2^expr67 * (expr20 * expr35) * expr36
    expr76 <- expr7^(expr32 - 1)
    expr82 <- expr11^2
    expr87 <- expr7^2
    expr96 <- expr7^(expr49 - 1)
    expr104 <- expr29^2
    expr109 <- 1/v
    expr112 <- expr1 * expr109
    expr117 <- expr5^expr20
    expr118 <- p * expr109
    expr119 <- expr117 * expr118
    expr121 <- expr76 * (expr32 * expr119)
    expr134 <- expr96 * (expr49 * expr119)
    expr138 <- expr26 * expr119
    expr139 <- expr50 * expr138
    expr140 <- expr20 * expr139
    expr145 <- log(expr2)
    expr146 <- expr3 * expr145
    expr147 <- log(expr5)
    expr148 <- expr6 * expr147
    expr149 <- expr146 + expr148
    expr152 <- p^2
    expr153 <- 1/expr152
    expr154 <- expr17 * expr153
    expr156 <- expr76 * (expr32 * expr149) - expr33 * expr154
    expr161 <- expr34 * expr145 * expr36 + expr34 * expr35
    expr163 <- expr153 * expr37
    expr172 <- 2/expr152
    expr180 <- expr96 * (expr49 * expr149) + expr50 * expr154
    expr189 <- expr26 * expr149
    expr192 <- expr50 * expr189 + expr27 * expr154
    expr194 <- expr27 + expr20 * expr192
    expr199 <- expr8 * expr119
    expr203 <- expr119/expr7
    expr206 <- expr112/expr21
    expr212 <- 1/v^2
    expr225 <- expr117 * (p * expr212) + expr5^expr67 * (expr20 * expr109) * expr118
    expr255 <- expr117 * expr147 * expr118 + expr117 * expr109
    expr257 <- expr153 * expr119
    expr282 <- expr149/expr7
    expr286 <- expr8 * expr149
    expr289 <- expr33 * expr286 - expr9 * expr154
    expr296 <- expr146 * expr145 + expr148 * expr147
    expr302 <- expr172 * expr282
    expr304 <- 2 * p
    expr306 <- expr152^2
    expr313 <- expr153 * expr149
    expr321 <- expr282 * expr153 - expr17 * (expr304/expr306)
    value <- -expr9 - log(expr11) + expr16 * expr17 + expr20 * expr22 + log(expr29)

    if(deriv>0){
      gradient <- matrix(0, nrow=N, ncol=3)
      hessian <- matrix(0, nrow=N, ncol=6)
      gradient[, 1] <- expr33 * expr38 - v/expr11 - expr16 * expr42 + expr20 * expr46 - expr53/expr29
      hessian[, 1] <- -(expr20 * (expr57 * expr4/expr21 + expr45 * expr45/expr61) + (expr33 * (expr8 * expr72) + expr76 * (expr32 * expr37) * expr38 - v * v/expr82 - expr16 * (expr72/expr7 - expr37 * expr37/expr87)) - (expr20 * (expr50 * (expr26 *
                                                                                                                                                                                                                                                                         expr72) + expr96 * (expr49 * expr37) * expr51)/expr29 - expr53 * expr53/expr104))
      hessian[, 2]  <- expr20 * (expr35 * expr109/expr21 - expr45 * expr112/expr61) - (expr121 * expr38 + (1/expr11 -  v * u/expr82) + expr16 * (expr37 * expr119/expr87)) + (expr20 * (expr134 * expr51)/expr29 - expr53 * expr140/expr104)
      hessian[, 3]  <- expr156 * expr38 + expr33 * (expr8 * expr161 - expr163) - (expr16 * (expr161/expr7 - expr37 * expr149/expr87) - expr172 * expr42) + expr46 - ((expr52 + expr20 * (expr180 * expr51 + expr50 * (expr163 + expr26 * expr161)))/expr29 - expr53 *  expr194/expr104)
      gradient[, 2] <- expr33 * expr199 - u/expr11 - expr16 *   expr203 + expr20 * expr206 - expr140/expr29
      hessian[, 4] <- -(expr20 * (expr1 *  expr212/expr21 + expr112 * expr112/expr61) + (expr33 * (expr8 * expr225) + expr121 * expr199 - u * u/expr82 -
                                                                                                   expr16 * (expr225/expr7 - expr119 * expr119/expr87)) - (expr20 * (expr50 * (expr26 * expr225) + expr134 * expr138)/expr29 - expr140 * expr140/expr104))
      hessian[, 5]  <- expr156 * expr199 + expr33 * (expr8 * expr255 - expr257) - (expr16 * (expr255/expr7 - expr119 * expr149/expr87) - expr172 * expr203) + expr206 - ((expr139 + expr20 * (expr180 * expr138 + expr50 * (expr257 + expr26 * expr255)))/expr29 - expr140 * expr194/expr104)
      gradient[, 3] <- expr16 * expr282 - expr172 * expr17 - expr289 + expr22 + expr194/expr29
      hessian[, 6] <- expr16 * (expr296/expr7 - expr149 * expr149/expr87) - expr302 - (expr302 -  2 * expr304/expr306 * expr17) - (expr156 * expr286 + expr33 * (expr8 * expr296 - expr313) - (expr289 * expr154 + expr9 * expr321)) + ((expr192 + (expr192 +  expr20 * (expr180 * expr189 + expr50 * (expr313 +  expr26 * expr296) + (expr192 * expr154 + expr27 *  expr321))))/expr29 - expr194 * expr194/expr104)
    }
  }

  if(family_cop=="frank"){
    #Check for correct inputs
    if(any(-35>p|p>35)){ #p in [-35,35]
      stop(paste("p must be in [-35,35]", "\n", ""))
    }

    expr2 <- -p
    expr3 <- exp(expr2)
    expr4 <- 1 - expr3
    expr7 <- u + v
    expr11 <- exp(expr2 * u)
    expr13 <- exp(expr2 * v)
    expr14 <- 1 - expr13
    expr17 <- expr4 - (1 - expr11 * expr14)
    expr18 <- expr17^2
    expr21 <- expr11 * p
    expr22 <- expr21 * expr14
    expr24 <- 2 * (expr22 * expr17)
    expr36 <- expr18^2
    expr40 <- expr13 * p
    expr43 <- expr11 * expr40
    expr49 <- 2 * (expr43 * expr17)
    expr53 <- expr11 * u
    expr57 <- expr13 * v
    expr64 <- expr3 + (expr11 * expr57 - expr53 * expr14)
    expr70 <- 2 * (expr64 * expr17)
    expr106 <- expr3/expr4
    expr121 <- expr53 * expr57
    value <- log(p) + log(expr4) - p * expr7 - log(expr18)

    if(deriv>0){
      gradient <- matrix(0, nrow=N, ncol=3)
      hessian <- matrix(0, nrow=N, ncol=6)

      gradient[, 1] <- -(p - expr24/expr18)
      hessian[, 1] <- -(2 * (expr22 * expr22 +  expr21 * p * expr14 * expr17)/expr18 - expr24 * expr24/expr36)
      hessian[, 2] <-  2 * (expr21 * expr40 * expr17 + expr22 * expr43)/expr18 - expr24 * expr49/expr36
      hessian[, 3] <-  -(1 - (2 * (((expr11 - expr53 * p) * expr14 + expr21 * expr57) * expr17 + expr22 * expr64)/expr18 - expr24 * expr70/expr36))
      gradient[, 2] <- -(p + expr49/expr18)
      hessian[, 4] <- -(2 * (expr43 * expr43 -  expr11 * (expr40 * p) * expr17)/expr18 - expr49 *  expr49/expr36)
      hessian[, 5] <- -(1 + (2 * ((expr11 * (expr13 - expr57 * p) - expr53 * expr40) * expr17 + expr43 * expr64)/expr18 - expr49 * expr70/expr36))
      gradient[, 3] <- 1/p + expr106 - expr7 - expr70/expr18
      hessian[,6] <- -(expr106 + expr3 * expr3/expr4^2 + 1/p^2 + (2 * (expr64 * expr64 - (expr11 * (expr57 * v) + expr121 + (expr121 - expr53 * u * expr14) + expr3) * expr17)/expr18 - expr70 * expr70/expr36))
    }
  }

  if(family_cop=="joe"){
    #Check for correct inputs
    if(any(1>=p|p>30)){ #p in (1,30]
      stop(paste("p must be in (1,30]", "\n", ""))
    }

    expr2 <- 1/p - 2
    expr3 <- 1 - u
    expr4 <- expr3^p
    expr6 <- -log(v)
    expr7 <- expr6^p
    expr9 <- expr4 * expr7
    expr10 <- expr4 + expr7 - expr9
    expr11 <- log(expr10)
    expr13 <- p - 1
    expr16 <- expr13 + expr4 + expr7 - expr9
    expr19 <- log(expr3)
    expr22 <- 1 - v
    expr23 <- log(expr22)
    expr26 <- 1/expr3
    expr28 <- expr3^expr13
    expr29 <- expr28 * p
    expr31 <- expr29 - expr29 * expr7
    expr33 <- expr31/expr10
    expr41 <- expr13 - 1
    expr44 <- expr3^expr41 * expr13 * p
    expr46 <- expr44 - expr44 * expr7
    expr48 <- expr31 * expr31
    expr49 <- expr10^2
    expr54 <- expr16^2
    expr60 <- expr6^expr13
    expr61 <- 1/v
    expr62 <- p * expr61
    expr63 <- expr60 * expr62
    expr64 <- expr29 * expr63
    expr67 <- expr63 - expr4 * expr63
    expr68 <- expr31 * expr67
    expr79 <- expr28 * expr19 * p + expr28
    expr81 <- log(expr6)
    expr82 <- expr7 * expr81
    expr85 <- expr79 - (expr79 * expr7 + expr29 * expr82)
    expr87 <- expr4 * expr19
    expr92 <- expr87 * expr7 + expr4 * expr82
    expr93 <- 1 + expr87 + expr82 - expr92
    expr99 <- expr87 + expr82 - expr92
    expr104 <- p^2
    expr105 <- 1/expr104
    expr111 <- 1/expr22
    expr114 <- expr67/expr10
    expr130 <- expr60 * (p * (1/v^2)) + expr6^expr41 * (expr13 * expr61) * expr62
    expr132 <- expr130 - expr4 * expr130
    expr134 <- expr67 * expr67
    expr147 <- expr60 * expr81 * expr62 + expr60 * expr61
    expr151 <- expr147 - (expr87 * expr63 + expr4 * expr147)
    expr166 <- expr99/expr10
    expr174 <- expr87 * expr19
    expr175 <- expr82 * expr81
    expr178 <- expr87 * expr82
    expr183 <- expr174 + expr175 - (expr174 * expr7 + expr178 + (expr178 + expr4 * expr175))
    expr189 <- expr105 * expr166
    value <- expr2 * expr11 + log(expr16) + expr13 * expr19 + expr13 * expr23

    if(deriv>0){
      gradient <- matrix(0, nrow=N, ncol=3)
      hessian <- matrix(0, nrow=N, ncol=6)
      gradient[, 1] <- -(expr13 * expr26 + (expr31/expr16 +  expr2 * expr33))
      hessian[, 1] <- -(expr13 * (1/expr3^2) - (expr2 * (expr46/expr10 - expr48/expr49) + (expr46/expr16 - expr48/expr54)))
      hessian[, 2]  <- -(expr64/expr16 +   expr68/expr54 + expr2 * (expr64/expr10 + expr68/expr49))
      hessian[, 3]  <- -(expr26 +   (expr85/expr16 - expr31 * expr93/expr54 + (expr2 * (expr85/expr10 - expr31 * expr99/expr49) - expr105 * expr33)))
      gradient[, 2] <- -(expr13 * expr111 + (expr67/expr16 + expr2 * expr114))
      hessian[, 4] <- -(expr13 * (1/expr22^2) - (expr2 * (expr132/expr10 - expr134/expr49) + (expr132/expr16 - expr134/expr54)))
      hessian[, 5]  <- -(expr111 + (expr151/expr16 - expr67 * expr93/expr54 + (expr2 * (expr151/expr10 - expr67 * expr99/expr49) - expr105 * expr114)))
      gradient[, 3] <- expr2 * expr166 - expr105 * expr11 + expr93/expr16 + expr19 + expr23
      hessian[, 6] <- expr2 * (expr183/expr10 - expr99 *  expr99/expr49) - expr189 - (expr189 - 2 * p/expr104^2 * expr11) + (expr183/expr16 - expr93 * expr93/expr54)
    }
  }

  if(!log.p){
    value<-exp(value)
  }

  if(deriv>0){
    attr(value, "gradient") <- gradient
    attr(value, "hessian") <- hessian
  }

  #return value
  return(value)
}

#' @describeIn dcop distribution function for copula.
#' @export
pcop<-function(W, delta=0, family_cop="normal", log.p = FALSE){
  #Distribution function of bivariate copula
  if(!is.matrix(delta)){
    if(length(delta)==N){
      delta<-matrix(delta, nrow=N, ncol=1)
    } else {
      delta<-matrix(delta, nrow=N, ncol=length(delta), byrow=TRUE)
    }
  }

  if(!is.matrix(W)){
    stop(paste("W must be a matrix", "\n", ""))
  }

  #Correct for numerical inaccuracy
  W[W>1-1e-16]<-1-1e-16
  W[W<0+1e-16]<-0+1e-16
  N<-nrow(W)

  #Check for correct inputs
  if (!family_cop%in%c("independent","normal","clayton","gumbel","frank","joe")){
    stop(paste("Incorrect family_cop", "\n", ""))
  }


  out<-sapply(1:N, function(n) pcop_copula(W=W[n,], delta=delta[n,], family_cop=family_cop, log.p=log.p))

  #Return out
  return(out)
}

#' @describeIn dcop random number generation for copula.
#' @param n number of observations.
#' @export
rcop<-function(n, delta=0, family_cop="normal"){
  #Random number generation function of copula

  #Check for correct inputs
  if (!family_cop%in%c("independent","normal","clayton","gumbel","frank","joe")){
    stop(paste("Incorrect family_cop", "\n", ""))
  }


  if(!is.matrix(delta)){
    if(length(delta)==n){
      delta<-matrix(delta, nrow=n, ncol=1)
    } else {
      delta<-matrix(delta, nrow=n, ncol=length(delta), byrow=TRUE)
    }
  }

  N<-n
  out<-lapply(1:N, function(i) rcop_copula(delta=delta[i,], family_cop=family_cop))
  out<-matrix(unlist(out), byrow=TRUE, nrow=N)

  return(out)
}

pcop_copula<-function(W, delta=0, family_cop="normal", log.p = FALSE){
  #dcop wrapper function for distribution function for copula with scalar inputs

  if(family_cop=="independent"){
    cop_object<-copula::indepCopula(param=rep(0,ncol(delta)), dim = 2)
  }

  if(family_cop=="normal"){
    cop_dim<-nrow(copula::p2P(c(delta)))
    cop_object<-copula::normalCopula(param=delta, dim = cop_dim)
  }

  if(family_cop=="clayton"){
    cop_object<-copula::claytonCopula(param=delta, dim = 2)
  }

  if(family_cop=="gumbel"){
    cop_object<-copula::gumbelCopula(param=delta, dim = 2)
  }

  if(family_cop=="frank"){
    cop_object<-copula::frankCopula(param=delta, dim = 2)
  }

  if(family_cop=="joe"){
    cop_object<-copula::joeCopula(param=delta, dim = 2)
  }

  #Wrapper for Copula package function
  out<-copula::pCopula(u=W, copula=cop_object)

  if(log.p){
    out<-log(out)
  }

  return(out)
}

rcop_copula<-function(delta=0, family_cop="normal"){
  #wrapper function for random number generation for copula with scalar inputs

  if(family_cop=="independent"){
    cop_object<-copula::indepCopula(param=rep(0,ncol(delta)), dim = 2)
  }

  if(family_cop=="normal"){
    cop_dim<-nrow(copula::p2P(c(delta)))
    cop_object<-copula::normalCopula(param=delta, dim = cop_dim)
  }

  if(family_cop=="clayton"){
    cop_object<-copula::claytonCopula(param=delta, dim = 2)
  }

  if(family_cop=="gumbel"){
    cop_object<-copula::gumbelCopula(param=delta, dim = 2)
  }

  if(family_cop=="frank"){
    cop_object<-copula::frankCopula(param=delta, dim = 2)
  }

  if(family_cop=="joe"){
    cop_object<-copula::joeCopula(param=delta, dim = 2)
  }

  #Wrapper for Copula package function
  out<-copula::rCopula(n=1, copula=cop_object)

  return(out)
}
