#' Weak instrument-robust inference for multivariable Mendelian randomization with two-sample summary data
#'
#' @description Compute weak instrument-robust confidence sets in multivariable Mendelian randomization with two-sample summary data.
#'
#' @param bx J x K matrix of genetic variant associations with K exposures (number of variants J cannot be less than number of exposures K)
#' @param sx J x K matrix of corresponding standard errors from variant-exposure regressions
#' @param by J vector of genetic variant associations with outcome
#' @param sy J vector of corresponding standard errors from variant-outcome regressions
#' @param nx sample size for genetic variant-exposure regressions; the exposure associations are assumed to be measured from the same sample
#' @param ny sample size for genetic variant-outcome regressions
#' @param ld (optional) J x J genetic variant correlation matrix; if not supplied, the default is to assume uncorrelated variants
#' @param cor.x (optional) K x K matrix exposure correlation matrix; if not supplied, the default is to assume uncorrelated exposures
#' @param alpha (optional) nominal significance level; 1-alpha is the nominal coverage probability of confidence sets
#' @param gamma.min (optional) minimum allowable coverage distortion for non-robust inference; if not supplied, the default is 0.05
#' @param max.search (optional) c(-max.search, max.search) are the limits for the search if a positive scalar is provided, or else a vector with two numbers may be provided so that c(max.search1,max.search2) are the limits for each exposure. If not supplied, default is max.search = 2
#' @param len.search (optional) number of potential parameter values for each exposure to be tested to form confidence sets. The runtime substantially increases with larger len.search, and exponentially so with the number of exposures. If not supplied, default is len.search=100
#' @param kappa.limit (optional) c(-kappa.limit, kappa.limit) are the limits for the overdispersion parameter estimate. Only relevant if \code{robust} is \code{TRUE} and for the Kleibergen-OH confidence set. If not supplied, default is kappa.limit=250
#' @param robust (optional) logical \code{TRUE}/\code{FALSE}. If \code{TRUE}, the Kleibergen-OH confidence set is computed, but runtime may take a lot longer. The default is \code{FALSE}
#' @param seed (optional) the seed for random sampling involved in computing critical values for the Andrews robust confidence set (default value is 100)
#'
#' @details Weak instrument-robust confidence sets for multivariable Mendelian randomization with two-sample summary data. \cr
#' \cr
#' The calculated distortion cutoff is the estimated maximum coverage distortion relating to only the non-robust confidence set (based on inverting Wald test statistics). \cr
#' \cr
#' The runtime for computing confidence sets substantially increases with larger values of len.search, especially with large numbers of instruments, and the runtime is exponentially longer with every additional exposure included in the model.
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{condFstat} \cr
#'  K vector of conditional F statistics (one for each exposure)
#'  \item \code{gmm_est} \cr
#'  K vector of causal effect estimates using GMM estimation
#'  \item \code{cs_n} \cr
#'  non-robust confidence set based on inverting Wald test statistics; each column is a K vector of exposure effect values not rejected the test
#'  \item \code{cs_ar} \cr
#'  robust confidence set based on inverting Anderson-Rubin test statistics; each column is a K vector of exposure effect values not rejected by the test
#'  \item \code{cs_k} \cr
#'  robust confidence set based on inverting Kleibergen test statistics; each column is a K vector of exposure effect values not rejected by the test
#'  \item \code{cs_koh} \cr
#'  robust confidence set based on inverting Kleibergen-OH test statistics that account for overdispersion heterogeneity; each column is a K vector of exposure effect values not rejected by the test. Reported only if the input \code{robust} is \code{TRUE}
#'  \item \code{cs_r} \cr
#'  robust confidence set based on inverting Andrews linear combination test statistics; each column is a K vector of exposure effect values not rejected by the test
#'  \item \code{cutoff} \cr
#'  maximum coverage distortion for non-robust Wald-based confidence set \code{cs_n}; minimum value is gamma.min
#'  \item \code{kappa} \cr
#'  the average estimated overdispersion heterogeneity parameter over non-rejected values of exposure effects. Reported only if the input \code{robust} is \code{TRUE}
#' }
#'
#' @examples res <- mvmr.weakiv(bx=cbind(dat$beta.map,dat$beta.pp),by=dat$beta.stroke,sx=cbind(dat$se.map,dat$se.pp),sy=dat$se.stroke,nx=dat$nx,ny=dat$ny,cor.x=dat$cor.x,max.search=0.3,len.search=80)
#' # this may take a few minutes to run
#'
#' # note that the 95% Anderson-Rubin confidence set res$cs_ar is empty. The Anderson-Rubin test doubles as an overidentification test, so the confidence set being empty may suggest excessive heterogeneity in genetic variant-outcome associations in the model
#'
#' # we can plot the 95% non-robust (Wald) and robust (Andrews) confidence sets
#' plot(res$cs_n[1,],res$cs_n[2,],pch=16,main="Non-robust (Wald) confidence set", xlab = "MAP effect", ylab = "PP effect", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#' grid();abline(v=0,lty=2);abline(h=0,lty=2)
#' plot(res$cs_r[1,],res$cs_r[2,],pch=16,main="Robust (Andrews) confidence set", xlab = "MAP effect", ylab = "PP effect", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))
#' grid();abline(v=0,lty=2);abline(h=0,lty=2)
#'
#' # the distortion cutoff res$cutoff is around 15%, which means we may value the evidence from the non-robust set if we are willing to accept that the non-robust set may have coverage of at least 80%
#'
#' @references Weak instruments in multivariable Mendelian randomization: methods and practice. Preprint.
#'
#' @author Ashish Patel and Stephen Burgess
#'
#' @export

mvmr.weakiv <- function(bx,by,sx,sy,nx,ny,ld=NULL,cor.x=NULL,alpha=NULL,gamma.min=NULL,max.search=NULL,len.search=NULL,robust=FALSE,kappa.limit=NULL,seed=100){
J = nrow(bx); K = ncol(bx)
if(missing(ld)){ld=diag(J)}
if(missing(cor.x)){cor.x=diag(K)}
if(missing(alpha)){alpha=0.05}
if(missing(gamma.min)){gamma.min=0.05}
if(missing(max.search)){max.search=2}
if(missing(len.search)){len.search=100}
if(missing(kappa.limit)){kappa.limit=250}

if( exists(".Random.seed") ) {
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
}
if (!is.na(seed)) { set.seed(seed) }

if(J<K){stop("number of instruments is less than the number of exposures.")}
if(robust==TRUE & J<=K){stop("The Kleibergen-OH method requires the number of instruments to be greater than the number of exposures")}

# exposures quantities of interest
ax <- matrix(NA,nrow=J,ncol=K); Ax <- list()
for (k in 1:K){ax[,k] <- 1/((nx*sx[,k]^2)+bx[,k]^2)}
for (k in 1:K){Ax[[k]] <- (sqrt(ax[,k])%*%t(sqrt(ax[,k])))*ld}
Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)
sqrt.Ax <- function(k){
  evec <- eigen(Ax[[k]])$vectors; eval <- eigen(Ax[[k]])$values
  return((evec%*%diag(sqrt(eval))%*%t(evec)))
}
sqrt.Ax <- lapply(1:K,sqrt.Ax)
SigX <- function(k,l){
  solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))*(cor.x[k,l]-(as.numeric(t(Bx[,k])%*%solve(sqrt.Ax[[k]]%*%t(sqrt.Ax[[l]]))%*%Bx[,l])))
}
SigX2 <- function(m){
  SigX2a <- list()
  for (m1 in 1:K){SigX2a[[m1]] <- SigX(m1,m)}
  SigX2a <- do.call(rbind, SigX2a)
  return(SigX2a)
}
SigX3 <- list()
for (m1 in 1:K){SigX3[[m1]] <- SigX2(m1)}
SigX3 <- do.call(cbind, SigX3)
# check: SigX3[(2*J+1):(3*J),(1*J+1):(2*J)] == SigX(3,2)
SigX <- SigX3; rm(SigX2, SigX3)
gamX_est <- function(k){as.vector(solve(Ax[[k]])%*%Bx[,k])}
gamX_est <- sapply(1:K,gamX_est)

# outcome quantities of interest
ay <- 1/((ny*sy^2)+by^2)
Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
By <- ay*by
SigY <- solve(Ay)*(1-as.numeric(t(By)%*%solve(Ay)%*%By))*(nx/ny)
gamY_est <- as.vector(solve(Ay)%*%By)


## compute conditional F-statistics
if(K>2){
  condF <- function(j){
    SigXX <- SigX[c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))])),c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))]))]
    g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]%*%tet))}
    Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
    tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
    Om.nr <- function(tet){as.matrix(cbind(diag(J),kronecker(t(-tet),diag(J)))%*%SigXX%*%t(cbind(diag(J),kronecker(t(-tet),diag(J)))))}
    Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
    G <- -gamX_est[,-j]
    DQ.nr <- function(tet){2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet))}
    condF <- nx*(nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(J-K+1))
    return(condF)
  }
  condF <- sapply(1:K,condF)
}

if(K==2){
  condF <- function(j){
    SigXX <- SigX[c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))])),c(((((j-1)*J)+1):(j*J)),((1:(J*K))[-((((j-1)*J)+1):(j*J))]))]
    g <- function(tet){as.vector(gamX_est[,j] - (gamX_est[,-j]*tet))}
    Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
    tet.gg <- nlminb(rep(0,(K-1)),objective=Q.gg)$par
    Om.nr <- function(tet){as.matrix(cbind(diag(J),kronecker(t(-tet),diag(J)))%*%SigXX%*%t(cbind(diag(J),kronecker(t(-tet),diag(J)))))}
    Q.nr <- function(tet){as.numeric(t(g(tet))%*%solve(Om.nr(tet))%*%g(tet))}
    G <- -gamX_est[,-j]
    DQ.nr <- function(tet){as.numeric(2*as.matrix(t(G)%*%solve(Om.nr(tet))%*%g(tet)))}
    condF <- nx*(nlminb(tet.gg,objective=Q.nr,gradient=DQ.nr)$objective/(J-K+1))
    return(condF)
  }
  condF <- sapply(1:K,condF)
}


## continuously-updating GMM inference

# preliminary estimate for initial value
g <- function(tet){as.vector(gamY_est - (gamX_est%*%tet))}
Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
tet.gg <- nlminb(rep(0,K),objective=Q.gg)$par

# GMM estimate
phi <- function(tet){kronecker(t(tet),diag(J))}
Om <- function(tet){as.matrix(SigY+(phi(tet)%*%SigX%*%t(phi(tet))))}
Om.r <- function(tet,kappa){as.matrix(SigY+(diag(J)*kappa*(nx/ny))+(phi(tet)%*%SigX%*%t(phi(tet))))}
Q <- function(tet){as.numeric(t(g(tet))%*%solve(Om(tet))%*%g(tet))}
G <- -gamX_est
DQ <- function(tet){2*as.matrix(t(G)%*%solve(Om(tet))%*%g(tet))}
gmm <- nlminb(tet.gg,objective=Q,gradient=DQ)$par
Sig.tet <- as.matrix(solve(t(G)%*%solve(Om(gmm))%*%G))
tet.bar <- function(tet){tet - (solve(t(G)%*%solve(Om(tet))%*%G)%*%(t(G)%*%solve(Om(tet))%*%g(tet)))}
W <- function(tet){nx*as.numeric(t(tet.bar(tet)-tet)%*%solve(Sig.tet)%*%(tet.bar(tet)-tet))}

## identification-robust inference

# create grid of potential effect values to evaluate when constructing confidence sets
if(length(max.search)==1)
{tet.seq <- DescTools::CombSet(seq(-max.search,max.search,max.search/len.search),K,repl=TRUE,ord=TRUE)}
if(length(max.search)==2)
{tet.seq <- DescTools::CombSet(seq(max.search[1],max.search[2],(max.search[2]-max.search[1])/len.search),K,repl=TRUE,ord=TRUE)}
tet.seq <- t(tet.seq)

if(robust==FALSE){
  
  # evaluate various test statistics over the grid of potential effect values
  test.stats <- function(a){
    theta0 = tet.seq[,a]
    Om0 <- Om(theta0)
    g0 <- g(theta0)
    del <- SigX%*%t(phi(theta0))
    ind <- matrix(1:(J*K),ncol=J,byrow=TRUE)
    D <- list()
    for (k in 1:K){
      D[[k]] <- G[,k]-(del[ind[k,],]%*%solve(Om0)%*%g0)
    }
    D0 <- do.call(cbind, D); rm(D)
    theta.star <- theta0 - (solve(t(D0)%*%solve(Om0)%*%D0)%*%t(D0)%*%solve(Om0)%*%g0)
    M0 <- solve(Om0)%*%D0%*%solve(t(D0)%*%solve(Om0)%*%D0)
    K0 <- nx*as.numeric(t(theta.star-theta0)%*%solve(t(M0)%*%Om0%*%M0)%*%(theta.star-theta0))
    K01 <- nx*as.numeric((t(g0)%*%solve(Om0)%*%D0)%*%solve(t(D0)%*%solve(Om0)%*%D0)%*%(t(D0)%*%solve(Om0)%*%g0))
    W0 <- W(theta0)
    S0 <- nx*Q(theta0)
    return(c(W0,S0,K0,K01))
  }
  
  # collect Wald, AR, and K statistic results
  test.stats <- sapply(1:ncol(tet.seq),test.stats)
  W0 <- test.stats[1,]; S0 <- test.stats[2,]; K0 <- test.stats[3,]; K01 <- test.stats[4,]; rm(test.stats)
  
  # for a gamma.min cutoff, calculate a_hat for the linear combination of K0 and S0
  chi1 <- rchisq(1e6,df=K); chi2 <- rchisq(1e6,df=(J-K))
  
  a_hat <- function(a){mean((((1+a)*chi1)+(a*chi2)) <= qchisq(1-alpha,K))-(1-alpha-gamma.min)}
  sol <- try(uniroot(a_hat,lower=0,upper=100))
  a_hat <- if (inherits(sol, "try-error")) 0.01 else sol$root
  
  # calculate robust and non-robust sets of Andrews (2018)
  cv.r <- quantile(((1+a_hat)*chi1)+(a_hat*chi2),1-alpha)[[1]]
  cs_robust <- which((K0+(a_hat*S0))<=cv.r)
  if(length(cs_robust)>0){cs_robust <- tet.seq[,cs_robust]} else {cs_robust <- NA}
  cs_nonrobust <- which(W0<=qchisq(1-alpha,K))
  if(length(cs_nonrobust)>0){cs_nonrobust <- tet.seq[,cs_nonrobust]} else {cs_nonrobust <- NA}
  
  # calculate distortion cutoff of Andrews (2018)
  a_tilde <- max(((qchisq(1-alpha,K)-K0)/S0)*ifelse(W0>qchisq(1-alpha,K),1,0))
  gamma_tilde <- 1-alpha-mean((((1+a_tilde)*chi1)+(a_tilde*chi2)) <= qchisq(1-alpha,K))
  gamma_hat <- max(gamma_tilde,gamma.min)
  
  # calculate Anderson-Rubin (1949) confidence sets
  ar <- which(S0<=qchisq(1-alpha,J))
  if(length(ar)>0){ar <- tet.seq[,ar]} else {ar <- NA}
  
  # calculate Kleibergen (2005) confidence sets
  kleibergen <- which(K01<=qchisq(1-alpha,K))
  if(length(kleibergen)>0){kleibergen <- tet.seq[,kleibergen]} else {kleibergen <- NA}
  
  
  # one-dimensional confidence intervals through the projection method
  tet.diff <- tet.seq[1,2]-tet.seq[1,1]
  
  # wald intervals
  if(sum(is.na(cs_nonrobust))>0){cs_nonrobust.ci <- NA}
  if(sum(is.na(cs_nonrobust))==0){
    cs_nonrobust1 <- function(s){unique(sort(cs_nonrobust[s,]))}; cs_nonrobust1 <- sapply(1:K,cs_nonrobust1)
    cs_nonrobust.ci <- function(s){
      cs_nonrobust.diff <- c(NA, diff(cs_nonrobust1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((cs_nonrobust.diff[-1]<(tet.diff+0.5*tet.diff))&(cs_nonrobust.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){cs_nonrobust.ci <- c(min(cs_nonrobust1[[s]]),max(cs_nonrobust1[[s]]))}
      if(length(jump)>0){
        cs_nonrobust.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        cs_nonrobust.ci[1,] <- cs_nonrobust1[[s]][c(1,(jump[1]-1))]
        cs_nonrobust.ci[(length(jump)+1),] <- cs_nonrobust1[[s]][c(jump[length(jump)],length(cs_nonrobust1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            cs_nonrobust.ci[(i+1),] <- cs_nonrobust1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(cs_nonrobust.ci)
    }
    cs_nonrobust.ci <- lapply(1:K,cs_nonrobust.ci)
  }
  
  # anderson-rubin intervals
  if(sum(is.na(ar))==1){ar.ci <- NA}
  if(sum(is.na(ar))==0){
    ar1 <- function(s){unique(sort(ar[s,]))}; ar1 <- sapply(1:K,ar1)
    ar.ci <- function(s){
      ar.diff <- c(NA, diff(ar1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((ar.diff[-1]<(tet.diff+0.5*tet.diff))&(ar.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){ar.ci <- c(min(ar1[[s]]),max(ar1[[s]]))}
      if(length(jump)>0){
        ar.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        ar.ci[1,] <- ar1[[s]][c(1,(jump[1]-1))]
        ar.ci[(length(jump)+1),] <- ar1[[s]][c(jump[length(jump)],length(ar1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            ar.ci[(i+1),] <- ar1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(ar.ci)
    }
    ar.ci <- lapply(1:K,ar.ci)
  }
  
  # kleibergen intervals
  if(sum(is.na(kleibergen))==1){kleibergen.ci <- NA}
  if(sum(is.na(kleibergen))==0){
    kleibergen1 <- function(s){unique(sort(kleibergen[s,]))}; kleibergen1 <- sapply(1:K,kleibergen1)
    kleibergen.ci <- function(s){
      kleibergen.diff <- c(NA, diff(kleibergen1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((kleibergen.diff[-1]<(tet.diff+0.5*tet.diff))&(kleibergen.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){kleibergen.ci <- c(min(kleibergen1[[s]]),max(kleibergen1[[s]]))}
      if(length(jump)>0){
        kleibergen.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        kleibergen.ci[1,] <- kleibergen1[[s]][c(1,(jump[1]-1))]
        kleibergen.ci[(length(jump)+1),] <- kleibergen1[[s]][c(jump[length(jump)],length(kleibergen1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            kleibergen.ci[(i+1),] <- kleibergen1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(kleibergen.ci)
    }
    kleibergen.ci <- lapply(1:K,kleibergen.ci)
  }
  
  # andrews intervals
  if(sum(is.na(cs_robust))==1){cs_robust.ci <- NA}
  if(sum(is.na(cs_robust))==0){
    cs_robust1 <- function(s){unique(sort(cs_robust[s,]))}; cs_robust1 <- sapply(1:K,cs_robust1)
    cs_robust.ci <- function(s){
      cs_robust.diff <- c(NA, diff(cs_robust1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((cs_robust.diff[-1]<(tet.diff+0.5*tet.diff))&(cs_robust.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){cs_robust.ci <- c(min(cs_robust1[[s]]),max(cs_robust1[[s]]))}
      if(length(jump)>0){
        cs_robust.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        cs_robust.ci[1,] <- cs_robust1[[s]][c(1,(jump[1]-1))]
        cs_robust.ci[(length(jump)+1),] <- cs_robust1[[s]][c(jump[length(jump)],length(cs_robust1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            cs_robust.ci[(i+1),] <- cs_robust1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(cs_robust.ci)
    }
    cs_robust.ci <- lapply(1:K,cs_robust.ci)
  }
  
  decimals <- function(number, places) format(round(number, places), nsmall = places)
  lim.reached <- matrix(0,nrow=K,ncol=4)
  for (k in 1:K){
    cat("\n\n-----------------------------------------------------------------------------------------\n")
    cat((1-alpha)*100,"% confidence interval (via the projection method) for exposure ",k," effect")
    cat("\n-----------------------------------------------------------------------------------------")
    cat("\nWald (non-robust):       ")
    if(sum(is.na(cs_nonrobust.ci))>0) {cat("empty")} else {
      if(is.vector(cs_nonrobust.ci[[k]])){cs_nonrobust.ci[[k]] <- t(as.matrix(cs_nonrobust.ci[[k]]))}
      for(j in 1:nrow(cs_nonrobust.ci[[k]])) cat("[", ifelse(cs_nonrobust.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(cs_nonrobust.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(cs_nonrobust.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(cs_nonrobust.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,1] <- sum(cs_nonrobust.ci[[k]]==min(tet.seq[k,]))+sum(cs_nonrobust.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nAnderson-Rubin:          ")
    if(sum(is.na(ar.ci))>0) {cat("empty")} else {
      if(is.vector(ar.ci[[k]])){ar.ci[[k]] <- t(as.matrix(ar.ci[[k]]))}
      for(j in 1:nrow(ar.ci[[k]])) cat("[", ifelse(ar.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(ar.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(ar.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(ar.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,2] <- sum(ar.ci[[k]]==min(tet.seq[k,]))+sum(ar.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nKleibergen:              ") 
    if(sum(is.na(kleibergen.ci))>0) {cat("empty")} else {
      if(is.vector(kleibergen.ci[[k]])){kleibergen.ci[[k]] <- t(as.matrix(kleibergen.ci[[k]]))}
      for(j in 1:nrow(kleibergen.ci[[k]])) cat("[", ifelse(kleibergen.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(kleibergen.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(kleibergen.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(kleibergen.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,3] <- sum(kleibergen.ci[[k]]==min(tet.seq[k,]))+sum(kleibergen.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nAndrews:                 ")
    if(sum(is.na(cs_robust.ci))>0) {cat("empty")} else {
      if(is.vector(cs_robust.ci[[k]])){cs_robust.ci[[k]] <- t(as.matrix(cs_robust.ci[[k]]))}
      for(j in 1:nrow(cs_robust.ci[[k]])) cat("[", ifelse(cs_robust.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(cs_robust.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(cs_robust.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(cs_robust.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,4] <- sum(cs_robust.ci[[k]]==min(tet.seq[k,]))+sum(cs_robust.ci[[k]]==max(tet.seq[k,]))
    }
  }
  cat("\n ")
  if(sum(is.na(ar.ci))>0){cat("\nThe Anderson-Rubin confidence set is empty, which indicates that no parameter values of exposure effects support the assumed model given the data.\n")}
  if(sum(is.na(cs_nonrobust.ci))>0){cat("\nThe Wald (non-robust) confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(is.na(kleibergen.ci))>0){cat("\nThe Kleibergen confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(is.na(cs_robust.ci))>0){cat("\nThe Andrews confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(lim.reached)>0){cat("\nConfidence intervals indicated with a '<' or '>' sign are constrained by the range of values considered, and so do not reflect the full interval. For a full interval, please increase the max.search parameter and try again.\n")}
  cat("\n ")
  res.list <- list("condFstat"=condF, "gmm_est"=gmm, "cs_ar"=ar, "cs_k"=kleibergen, "cs_n"=cs_nonrobust, "cs_r"=cs_robust, "cutoff"=gamma_hat)
  return(res.list)
}

if(robust==TRUE){
  
  # evaluate various test statistics over the grid of potential effect values
  test.stats <- function(a){
    theta0 = tet.seq[,a]
    Om0 <- Om(theta0)
    Om.r0 <- function(kappa){Om.r(theta0,kappa)}
    g0 <- g(theta0)
    del <- SigX%*%t(phi(theta0))
    ind <- matrix(1:(J*K),ncol=J,byrow=TRUE)
    D <- list()
    for (k in 1:K){
      D[[k]] <- G[,k]-(del[ind[k,],]%*%solve(Om0)%*%g0)
    }
    D0 <- do.call(cbind, D); rm(D)
    theta.star <- theta0 - (solve(t(D0)%*%solve(Om0)%*%D0)%*%t(D0)%*%solve(Om0)%*%g0)
    M0 <- solve(Om0)%*%D0%*%solve(t(D0)%*%solve(Om0)%*%D0)
    K0 <- nx*as.numeric(t(theta.star-theta0)%*%solve(t(M0)%*%Om0%*%M0)%*%(theta.star-theta0))
    K01 <- nx*as.numeric((t(g0)%*%solve(Om0)%*%D0)%*%solve(t(D0)%*%solve(Om0)%*%D0)%*%(t(D0)%*%solve(Om0)%*%g0))
    W0 <- W(theta0)
    S0 <- nx*Q(theta0)
    Q.kap <- function(kappa){(as.numeric(nx*(t(g0)%*%solve(Om.r0(kappa))%*%g0))-J)^2}
    kappa.est <- optim(0,Q.kap,method="Brent",lower=-kappa.limit,upper=kappa.limit)$par
    Om.r0 <- Om.r0(max(kappa.est,0))
    D.r <- list()
    for (k in 1:K){
      D.r[[k]] <- G[,k]-(del[ind[k,],]%*%solve(Om.r0)%*%g0)
    }
    D0.r <- do.call(cbind, D.r); rm(D.r)
    K0.r <- nx*as.numeric((t(g0)%*%solve(Om.r0)%*%D0.r)%*%solve(t(D0.r)%*%solve(Om.r0)%*%D0.r)%*%(t(D0.r)%*%solve(Om.r0)%*%g0))
    return(c(W0,S0,K0,K01,K0.r,kappa.est))
  }
  
  
  # collect Wald, AR, and K statistic results
  test.stats <- sapply(1:ncol(tet.seq),test.stats)
  W0 <- test.stats[1,]; S0 <- test.stats[2,]; K0 <- test.stats[3,]; K01 <- test.stats[4,]; K0.r <- test.stats[5,]; kappa.est <- test.stats[6,]; rm(test.stats)
  
  # for a gamma.min cutoff, calculate a_hat for the linear combination of K0 and S0
  chi1 <- rchisq(1e6,df=K); chi2 <- rchisq(1e6,df=(J-K))
  
  a_hat <- function(a){mean((((1+a)*chi1)+(a*chi2)) <= qchisq(1-alpha,K))-(1-alpha-gamma.min)}
  sol <- try(uniroot(a_hat,lower=0,upper=100))
  a_hat <- if (inherits(sol, "try-error")) 0.01 else sol$root
  
  # calculate robust and non-robust sets of Andrews (2018)
  cv.r <- quantile(((1+a_hat)*chi1)+(a_hat*chi2),1-alpha)[[1]]
  cs_robust <- which((K0+(a_hat*S0))<=cv.r)
  if(length(cs_robust)>0){cs_robust <- tet.seq[,cs_robust]} else {cs_robust <- NA}
  cs_nonrobust <- which(W0<=qchisq(1-alpha,K))
  if(length(cs_nonrobust)>0){cs_nonrobust <- tet.seq[,cs_nonrobust]} else {cs_nonrobust <- NA}
  
  # calculate distortion cutoff of Andrews (2018)
  a_tilde <- max(((qchisq(1-alpha,K)-K0)/S0)*ifelse(W0>qchisq(1-alpha,K),1,0))
  gamma_tilde <- 1-alpha-mean((((1+a_tilde)*chi1)+(a_tilde*chi2)) <= qchisq(1-alpha,K))
  gamma_hat <- max(gamma_tilde,gamma.min)
  
  # calculate Anderson-Rubin (1949) confidence sets
  ar <- which(S0<=qchisq(1-alpha,J))
  if(length(ar)>0){ar <- tet.seq[,ar]} else {ar <- NA}
  
  # calculate Kleibergen (2005) confidence sets
  kleibergen <- which(K01<=qchisq(1-alpha,K))
  if(length(kleibergen)>0){kleibergen <- tet.seq[,kleibergen]} else {kleibergen <- NA}
  
  # calculate Kleibergen-OH confidence sets
  kleibergen.r <- which(K0.r<=qchisq(1-alpha,K))
  if(length(kleibergen.r)>0){kappa.est <- kappa.est[kleibergen.r]; kleibergen.r <- tet.seq[,kleibergen.r]} else {kleibergen.r <- NA; kappa.est <- NA}
  
  # one-dimensional confidence intervals through the projection method
  tet.diff <- tet.seq[1,2]-tet.seq[1,1]
  
  # wald intervals
  if(sum(is.na(cs_nonrobust))>0){cs_nonrobust.ci <- NA}
  if(sum(is.na(cs_nonrobust))==0){
    cs_nonrobust1 <- function(s){unique(sort(cs_nonrobust[s,]))}; cs_nonrobust1 <- sapply(1:K,cs_nonrobust1)
    cs_nonrobust.ci <- function(s){
      cs_nonrobust.diff <- c(NA, diff(cs_nonrobust1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((cs_nonrobust.diff[-1]<(tet.diff+0.5*tet.diff))&(cs_nonrobust.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){cs_nonrobust.ci <- c(min(cs_nonrobust1[[s]]),max(cs_nonrobust1[[s]]))}
      if(length(jump)>0){
        cs_nonrobust.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        cs_nonrobust.ci[1,] <- cs_nonrobust1[[s]][c(1,(jump[1]-1))]
        cs_nonrobust.ci[(length(jump)+1),] <- cs_nonrobust1[[s]][c(jump[length(jump)],length(cs_nonrobust1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            cs_nonrobust.ci[(i+1),] <- cs_nonrobust1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(cs_nonrobust.ci)
    }
    cs_nonrobust.ci <- lapply(1:K,cs_nonrobust.ci)
  }
  
  # anderson-rubin intervals
  if(sum(is.na(ar))==1){ar.ci <- NA}
  if(sum(is.na(ar))==0){
    ar1 <- function(s){unique(sort(ar[s,]))}; ar1 <- sapply(1:K,ar1)
    ar.ci <- function(s){
      ar.diff <- c(NA, diff(ar1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((ar.diff[-1]<(tet.diff+0.5*tet.diff))&(ar.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){ar.ci <- c(min(ar1[[s]]),max(ar1[[s]]))}
      if(length(jump)>0){
        ar.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        ar.ci[1,] <- ar1[[s]][c(1,(jump[1]-1))]
        ar.ci[(length(jump)+1),] <- ar1[[s]][c(jump[length(jump)],length(ar1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            ar.ci[(i+1),] <- ar1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(ar.ci)
    }
    ar.ci <- lapply(1:K,ar.ci)
  }
  
  # kleibergen intervals
  if(sum(is.na(kleibergen))==1){kleibergen.ci <- NA}
  if(sum(is.na(kleibergen))==0){
    kleibergen1 <- function(s){unique(sort(kleibergen[s,]))}; kleibergen1 <- sapply(1:K,kleibergen1)
    kleibergen.ci <- function(s){
      kleibergen.diff <- c(NA, diff(kleibergen1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((kleibergen.diff[-1]<(tet.diff+0.5*tet.diff))&(kleibergen.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){kleibergen.ci <- c(min(kleibergen1[[s]]),max(kleibergen1[[s]]))}
      if(length(jump)>0){
        kleibergen.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        kleibergen.ci[1,] <- kleibergen1[[s]][c(1,(jump[1]-1))]
        kleibergen.ci[(length(jump)+1),] <- kleibergen1[[s]][c(jump[length(jump)],length(kleibergen1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            kleibergen.ci[(i+1),] <- kleibergen1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(kleibergen.ci)
    }
    kleibergen.ci <- lapply(1:K,kleibergen.ci)
  }
  
  # kleibergen-OH intervals
  if(sum(is.na(kleibergen.r))==1){kleibergen.r.ci <- NA}
  if(sum(is.na(kleibergen.r))==0){
    kleibergen.r1 <- function(s){unique(sort(kleibergen.r[s,]))}; kleibergen.r1 <- sapply(1:K,kleibergen.r1)
    kleibergen.r.ci <- function(s){
      kleibergen.r.diff <- c(NA, diff(kleibergen.r1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((kleibergen.r.diff[-1]<(tet.diff+0.5*tet.diff))&(kleibergen.r.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){kleibergen.r.ci <- c(min(kleibergen.r1[[s]]),max(kleibergen.r1[[s]]))}
      if(length(jump)>0){
        kleibergen.r.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        kleibergen.r.ci[1,] <- kleibergen.r1[[s]][c(1,(jump[1]-1))]
        kleibergen.r.ci[(length(jump)+1),] <- kleibergen.r1[[s]][c(jump[length(jump)],length(kleibergen.r1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            kleibergen.r.ci[(i+1),] <- kleibergen.r1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(kleibergen.r.ci)
    }
    kleibergen.r.ci <- lapply(1:K,kleibergen.r.ci)
  }
  
  # andrews intervals
  if(sum(is.na(cs_robust))==1){cs_robust.ci <- NA}
  if(sum(is.na(cs_robust))==0){
    cs_robust1 <- function(s){unique(sort(cs_robust[s,]))}; cs_robust1 <- sapply(1:K,cs_robust1)
    cs_robust.ci <- function(s){
      cs_robust.diff <- c(NA, diff(cs_robust1[[s]], lag = 1, differences = 1))
      jump <- (c(TRUE,((cs_robust.diff[-1]<(tet.diff+0.5*tet.diff))&(cs_robust.diff[-1]>(tet.diff-0.5*tet.diff)))))
      jump <- which(jump==FALSE)
      if(length(jump)==0){cs_robust.ci <- c(min(cs_robust1[[s]]),max(cs_robust1[[s]]))}
      if(length(jump)>0){
        cs_robust.ci <- matrix(,ncol=2,nrow=(length(jump)+1))
        cs_robust.ci[1,] <- cs_robust1[[s]][c(1,(jump[1]-1))]
        cs_robust.ci[(length(jump)+1),] <- cs_robust1[[s]][c(jump[length(jump)],length(cs_robust1[[s]]))]
        if(length(jump)>1){
          for (i in 1:(length(jump)-1)){
            cs_robust.ci[(i+1),] <- cs_robust1[[s]][c(jump[i],(jump[(i+1)]-1))]
          }
        }
      }
      return(cs_robust.ci)
    }
    cs_robust.ci <- lapply(1:K,cs_robust.ci)
  }
  
  decimals <- function(number, places) format(round(number, places), nsmall = places)
  lim.reached <- matrix(0,nrow=K,ncol=5)
  for (k in 1:K){
    cat("\n\n-----------------------------------------------------------------------------------------\n")
    cat((1-alpha)*100,"% confidence interval (via the projection method) for exposure ",k," effect")
    cat("\n-----------------------------------------------------------------------------------------")
    cat("\nWald (non-robust):       ")
    if(sum(is.na(cs_nonrobust.ci))>0) {cat("empty")} else {
      if(is.vector(cs_nonrobust.ci[[k]])){cs_nonrobust.ci[[k]] <- t(as.matrix(cs_nonrobust.ci[[k]]))}
      for(j in 1:nrow(cs_nonrobust.ci[[k]])) cat("[", ifelse(cs_nonrobust.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(cs_nonrobust.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(cs_nonrobust.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(cs_nonrobust.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,1] <- sum(cs_nonrobust.ci[[k]]==min(tet.seq[k,]))+sum(cs_nonrobust.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nAnderson-Rubin:          ")
    if(sum(is.na(ar.ci))>0) {cat("empty")} else {
      if(is.vector(ar.ci[[k]])){ar.ci[[k]] <- t(as.matrix(ar.ci[[k]]))}
      for(j in 1:nrow(ar.ci[[k]])) cat("[", ifelse(ar.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(ar.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(ar.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(ar.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,2] <- sum(ar.ci[[k]]==min(tet.seq[k,]))+sum(ar.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nKleibergen:              ") 
    if(sum(is.na(kleibergen.ci))>0) {cat("empty")} else {
      if(is.vector(kleibergen.ci[[k]])){kleibergen.ci[[k]] <- t(as.matrix(kleibergen.ci[[k]]))}
      for(j in 1:nrow(kleibergen.ci[[k]])) cat("[", ifelse(kleibergen.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(kleibergen.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(kleibergen.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(kleibergen.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,3] <- sum(kleibergen.ci[[k]]==min(tet.seq[k,]))+sum(kleibergen.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nAndrews:                 ")
    if(sum(is.na(cs_robust.ci))>0) {cat("empty")} else {
      if(is.vector(cs_robust.ci[[k]])){cs_robust.ci[[k]] <- t(as.matrix(cs_robust.ci[[k]]))}
      for(j in 1:nrow(cs_robust.ci[[k]])) cat("[", ifelse(cs_robust.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(cs_robust.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(cs_robust.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(cs_robust.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,4] <- sum(cs_robust.ci[[k]]==min(tet.seq[k,]))+sum(cs_robust.ci[[k]]==max(tet.seq[k,]))
    }
    cat("\nKleibergen-OH:           ") 
    if(sum(is.na(kleibergen.r.ci))>0) {cat("empty")} else {
      if(is.vector(kleibergen.r.ci[[k]])){kleibergen.r.ci[[k]] <- t(as.matrix(kleibergen.r.ci[[k]]))}
      for(j in 1:nrow(kleibergen.r.ci[[k]])) cat("[", ifelse(kleibergen.r.ci[[k]][j,1]==min(tet.seq[k,]), "<", ""), decimals(kleibergen.r.ci[[k]][j,1], max(ceiling(-log10(tet.diff)), 1)),",", ifelse(kleibergen.r.ci[[k]][j,2]==max(tet.seq[k,]), ">", ""), decimals(kleibergen.r.ci[[k]][j,2],max(ceiling(-log10(tet.diff)))),"]  ")
      lim.reached[k,5] <- sum(kleibergen.r.ci[[k]]==min(tet.seq[k,]))+sum(kleibergen.r.ci[[k]]==max(tet.seq[k,]))
    }
  }
  cat("\n ")
  if(sum(is.na(ar.ci))>0){cat("\nThe Anderson-Rubin confidence set is empty, which indicates that no parameter values of exposure effects support the assumed model given the data.\n")}
  if(sum(is.na(cs_nonrobust.ci))>0){cat("\nThe Wald (non-robust) confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(is.na(kleibergen.ci))>0){cat("\nThe Kleibergen confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(is.na(cs_robust.ci))>0){cat("\nThe Andrews confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(is.na(kleibergen.r.ci))>0){cat("\nThe Kleibergen-OH confidence set is empty. Consider increasing len.search and/or max.search parameters.\n")}
  if(sum(lim.reached)>0){cat("\nConfidence intervals indicated with a '<' or '>' sign are constrained by the range of values considered, and so do not reflect the full interval. For a full interval, please increase the max.search parameter and try again.\n")}
  cat("\n ")
  res.list <- list("condFstat"=condF, "gmm_est"=gmm, "cs_ar"=ar, "cs_k"=kleibergen, "cs_koh"=kleibergen.r, "cs_n"=cs_nonrobust, "cs_r"=cs_robust, "cutoff"=gamma_hat, "kappa"=mean(kappa.est))
  return(res.list)
}
}
