## Weak instrument-robust inference in multivariable Mendelian randomization with two-sample summary data ##
#### by Ashish Patel, James Lane, & Stephen Burgess #### 

##### Installing the R package
* devtools::install_github("ash-res/mvmr.weakiv")

##### Inputs (genetic association summary data)
 * bx: J x K matrix of genetic variant associations with K exposures (number of variants J cannot be less than number of exposures K)
 * by: J vector of genetic variant associations with outcome
 * sx: J x K matrix of corresponding standard errors from variant-exposure regressions
 * sy: J vector of corresponding standard errors from variant-outcome regressions
 * nx: sample size for genetic variant-exposure regressions; the exposure associations are assumed to be measured from the same sample
 * ny: sample size for genetic variant-outcome regressions
 * ld (optional): J x J genetic variant correlation matrix; if not supplied, the default is to assume uncorrelated variants
 * cor.x (optional): K x K matrix exposure correlation matrix; if not supplied, the default is to assume uncorrelated exposures
 * alpha (optional): nominal significance level; 1-alpha is the nominal coverage probability of confidence sets; if not supplied, the default is 0.05
 * gamma.min (optional): minimum allowable coverage distortion for non-robust inference; if not supplied, the default is 0.05
 * max.search (optional): c(-max.search, max.search) are the limits for the search if a positive scalar is provided, or else a vector with two numbers may be provided so that c(max.search1,max.search2) are the limits for each exposure. If not supplied, default is max.search = 2
 * len.search (optional): number of potential parameter values for each exposure to be tested to form confidence sets. The runtime substantially increases with larger len.search, and exponentially so with the number of exposures. If not supplied, default is len.search=100
 * robust (optional): logical TRUE/FALSE. If TRUE, the Kleibergen-OH confidence set is computed, but runtime may take a lot longer. The default is FALSE
 * kappa.limit (optional): c(-kappa.limit, kappa.limit) are the limits for the overdispersion parameter estimate. Only relevant if robust is TRUE and for the Kleibergen-OH confidence set. If not supplied, default is kappa.limit=250
 * seed (optional): the seed for random sampling involved in computing critical values for the Andrews robust confidence set (default value is 100)

