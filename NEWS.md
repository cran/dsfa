*   denotes small changes
**  denotes quite substantial/important changes
*** denotes really big changes 
inp "x" -> inp "y" denotes that input argument x has changed to input argument y

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