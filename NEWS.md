*   denotes small changes
**  denotes quite substantial/important changes
*** denotes really big changes 
inp "x" -> inp "y" denotes that input argument x has changed to input argument y

## version 2.0.2
* cop_distr -> distr to be more consistent for copulas

* Changed output of elasticity and efficiency function such that confidence intervals are now included.

* comper predict function bug fixed, defined npar.

* Added BurkinaFarms data set, which is used as an example in comper_mv.

* Bounded the inputs for W matrix to be in [0,1]^2 via clamp for more stability

* Improved starting values for the copula parameters via mgcv_cop

* Added "normalized" residuals, which are probability integral transformed residuals which are then transformed via qnorm()

* Introduced a wrapper function for mgcv::gam called "dsfa" in order to be in line with the naming convention of R packages

* Changed link function in comper and comper_mv from log to logshift. This is done for more numerical stability.

## version 2.0.1
* Added rotations 90, 180 and 270 for archimedian copulas and check parameter bounds of delta for dcop, pcop and rcop

* reworked elasticity() output to make it more consistent with efficiency()

* Added cost function example for dsfa and comper

## version 2.0.0
*** Major overhaul of the package.

* dnormhnorm, pnormhnorm, dnormexp, pnormexp, dcop are now wrappers for functions written via RcppArmadillo.
As these are called in the optimization the computation speed could be significantly increased. 

* Reducing dependencies to other packages

* Streamlining naming, comperr -> comper to be in line with erf.

* Adding examples from 'RiceFarms' data set

* Added data set 'manuf'

* Reparametrize is now mom2par and par2mom

## version 1.0.1
* Added rsp function (inverse link function)

* Added pnorm and qnorm functions which handle like the base R function but allow for the derivatives 0 and 2.

* OwenT function has now the option for deriv = 2.

* Added function "check_arguments" to check input arguments for pdf, cdf, quantile function and random number generation of normhnorm and normexp.

** Rewrote "chainrule" function such that the output is a vector with attributes "gradient", "hessian", "l3" and "l4" to be more consistent with the other functions. Further input
   g can be a list now.

* normexp and normhnorm can now take values of the parameters (mu, sigma_v, sigma_u, lambda) which differ in length from x,p,q,n. They are recycled to the appropriate length.

* Added functions "outlier_correct_column" and "outlier_correct" in "utility_functions" which are used to make the calculation of the gradient and hessian of copula more stable.

* Correct loglike of "joe" copula

* Simplified inputs for reparametrize.

* Added postproc for normexp_mgcv, normhnorm_mgcv and comperr_mv_mgcv to add attribute "s" to predicted values, such that qq.gam can now access information on s

* Bivariate normal copula compuation is now faster

* normexp_mgcv, normhnorm_mgcv and comperr_mv_mgcv have now simplified starting values


## version 1.0.0
*** inp"xg" -> inp "tri" in all functions, to be consistent with 'mgcv'

*** dcop, pcop and rcop have been rewritten and allow for independent, clayton, gumbel, frank and joe copula
  inp "U" -> inp "W"
  inp "Tau" -> inp "delta", delta is the copula parameter not kendalls tau anymore
  inp "family" -> inp "family_cop"
  Removed option for disjoint and numerical derivatives.
  Adjusted documentation

* removed cop_utility.R and moved functions into cop.R

*** OwenT function added for less dependencies from "sn" package

*** zeta function added for less dependencies from "sn" package

* erf.R added which calculates error function, complementary error function, inverse error function and inverse complementary error function.
  Further, d1qnormdx1 and d2qnormdx2 were also moved the same way. Formerly these were in utility_functions.R

* Rewrote the help page for 'dsfa'. Added an example.

## version 0.0.3
* initial release version