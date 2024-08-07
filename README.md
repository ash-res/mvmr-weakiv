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

##### Output is a list containing...
 * condFstat: K vector of conditional F statistics (one for each exposure)
 * gmm_est: K vector of causal effect estimates using GMM estimation
 * cs_n: non-robust confidence set based on inverting Wald test statistics; each column is a K vector of exposure effect values not rejected the test
 * cs_ar: robust confidence set based on inverting Anderson-Rubin test statistics; each column is a K vector of exposure effect values not rejected by the test
 * cs_k: robust confidence set based on inverting Kleibergen test statistics; each column is a K vector of exposure effect values not rejected by the test
 * cs_koh: robust confidence set based on inverting Kleibergen-OH test statistics that account for overdispersion heterogeneity; each column is a K vector of exposure effect values not rejected by the test. Reported only if the input robust is TRUE
 * cs_r: robust confidence set based on inverting Andrews linear combination test statistics; each column is a K vector of exposure effect values not rejected by the test
 * cutoff: maximum coverage distortion for non-robust Wald-based confidence set cs_n; minimum value is gamma.min
 * kappa: the average estimated overdispersion heterogeneity parameter over non-rejected values of exposure effects. Reported only if the input robust is TRUE

##### Example
 * res <- mvmr.weakiv(bx=cbind(dat$beta.map,dat$beta.pp),by=dat$beta.stroke,sx=cbind(dat$se.map,dat$se.pp),sy=dat$se.stroke,nx=dat$nx,ny=dat$ny,cor.x=dat$cor.x,max.search=0.3,len.search=80)
  (this may take a few minutes to run)
 * plot(res$cs_n[1,],res$cs_n[2,],pch=16,main="Non-robust (Wald) confidence set", xlab = "MAP effect", ylab = "PP effect", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
   grid();abline(v=0,lty=2);abline(h=0,lty=2)
 * plot(res$cs_r[1,],res$cs_r[2,],pch=16,main="Robust (Andrews) confidence set", xlab = "MAP effect", ylab = "PP effect", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
   grid();abline(v=0,lty=2);abline(h=0,lty=2)

   

##### Paper
 * for further details, please see the preprint: "Weak instruments in multivariable Mendelian randomization: methods and practice"
