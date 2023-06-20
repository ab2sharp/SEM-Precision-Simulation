#Robust Regression and SEM
library(Rfast)
library(R.utils)

# Initialization ----------------------------------------------------------

rpar <- function(p = NULL, signu=NULL)
{
  ### generate random parameters
  
  par = vector("list")
  par$beta = rnorm(p + 1) 
  if(is.null(signu))
  {
    par$sigma = rchisq(1, df = 1)
    par$nu = runif(1,1,5)
  }else
  {
    par$sigma = signu$sigma
    par$nu = signu$nu
  }
  return(par)
}

rbpar = function( p=NULL, sigma=NULL, nu=NULL, var=1 )
{
  par=list()
  
  par$beta = rnorm(p+1, sd=var)
  
  if( is.null(sigma) ) sigma = rchisq(1, df = 1)
  if( is.null(nu) ) nu = 1
  
  par$sigma = sigma
  par$nu = nu
  
  return( par )
  
}


# Likelihood --------------------------------------------------------------

robust.loglik <- function(par=NULL, x=NULL, y=NULL)
{
  # density for robust regression.  x is a matrix of covariates returns the
  # log-lik
  
  t = as.numeric(y - x %*% par$beta)/par$sigma
  val = sum(dt(t, df = par$nu, log = TRUE) - log(par$sigma))
  return(val)
}

rb.ll.betas <- function( 
  beta1=NULL, beta2=NULL, xd=NULL, yd=NULL, sigma=NULL, nu=NULL
) {
  # density for robust regression.  x is a matrix of covariates returns the
  # log-lik
  betas = c(beta1,beta2)
  t     = as.numeric(yd - xd %*% betas)/sigma
  val   = sum(dt(t, df = nu, log = TRUE) - log(sigma))
  return(val)
}

rb.ll.betav <- function(beta=NULL, xd=NULL, yd=NULL, signu=NULL)
{
  # density for robust regression.  x is a matrix of covariates returns the
  # log-lik
  # Dimension of yd, xd, and beta must match
  
  #Rename for readability
  sigma = signu$sigma
  nu    = signu$nu
  
  
  #loglik1
  t = as.numeric(yd - xd %*% beta)/sigma
  val1 = sum(dt(t, df = nu, log = TRUE) - log(sigma))
  
  return( val1 )
}

loglik.nu <- function(nu0=NULL, x=NULL, y=NULL, par0=NULL)
{
  # nu0 is a vector of nu0; MLE is the MLE of all the parameters; returns val
  # which is the profile log-likelihood at each nu
  val = numeric(length(nu0))
  for (i in 1:length(nu0)) {
    par0$nu = nu0[i]
    temp = EMn( par0=par0, x=x, y=y, n=1000, updatenu=FALSE )
    
    par0 = temp$par
    val[i] = temp$maxloglik
  }
  return(val)
}

q2.fn <- function(nu=NULL, u.hat=NULL, nu.old=NULL, wt=NULL)
{
  ## wt are weights
  if (is.null(wt)) 
    wt = rep(1, length(u.hat))
  
  n = sum(wt)
  sum.log.u = sum((log(u.hat) - u.hat) * wt)
  c1 = digamma((nu.old + 1)/2) - log((nu.old + 1)/2)
  
  val = -2 * lgamma(nu/2) + nu * log(nu/2) + nu * (c1 + sum.log.u/n)
  
  return(val)
}

q2.fn1 <- function(nu=NULL, u.hat=NULL, nu.old=NULL, wt=NULL)
{
  # The derivative of the q2.fn
  if (is.null(wt)) 
    wt = rep(1, length(u.hat))
  
  sum.log.u = sum((log(u.hat) - u.hat) * wt)/sum(wt)
  
  val = 1 - digamma(nu/2) + log(nu/2) + sum.log.u + digamma((nu.old + 1)/2) - 
    log((nu.old + 1)/2)
  return(val)
  
}



# Algorithm ---------------------------------------------------------------

EM.step <- function( par=NULL, x=NULL, y=NULL, wt=NULL, updatenu=TRUE )
{
  ## perform one iteration of EM; wt are weights
  
  if (is.null(wt)){ wt = rep(1, length(y)) }
  
  t2 = as.numeric(y - x %*% par$beta)^2/par$sigma^2
  ## E-step
  u.hat = (par$nu + 1)/(par$nu + t2)
  
  Xtw = t(sweep(x, 1, u.hat * wt, FUN = "*"))
  par$beta = solve(Xtw %*% x) %*% (Xtw %*% y)
  
  n = sum(wt)
  par$sigma = sqrt(sum(wt * u.hat * (y - x %*% par$beta)^2)/n)
  
  # should we update nu?
  if (updatenu) {
    par$nu = optimize(
      q2.fn, interval = c(1, 1000), u.hat = u.hat, nu.old = par$nu, 
      wt = wt, maximum = TRUE
    )$maximum
    
  }
  
  return(par)
}

EMn <- function( par0=NULL, x=NULL, y=NULL, n=NULL, updatenu=TRUE )
{
  ## perform n iterations of EM using the par0 starting value
  
  loglik = numeric(n)
  par = par0
  for (i in 1:length(loglik)) {
    par = EM.step(par = par, x = x, y = y, updatenu = updatenu)
    loglik[i] = robust.loglik(par = par, x = x, y = y)
  }
  
  val = list(loglik = loglik, par = par)
  return(val)
}

EMnv = function( 
  par0=NULL, x=NULL, y=NULL, n=NULL, updatenu=TRUE, EM.step1=NULL 
) {
  ## perform n iterations of EM using the par0 starting value
  
  loglik = numeric(n)
  par = par0
  for (i in 1:length(loglik)) {
    par = EM.step1(par = par, x = x, y = y, updatenu = updatenu)
    loglik[i] = robust.loglik(par = par, x = x, y = y)
  }
  
  val = list(loglik = loglik, par = par)
  return(val)
}

EMc <- function(
  par0 = NULL, x = NULL, y = NULL, atol = NULL, max.iter = 1000, updatenu = TRUE
) {
  ## perform n iterations of EM using the par0 starting value
  par = par0
  loglik = numeric(max.iter)
  for (i in 1:3) {
    par = EM.step(par = par, x = x, y = y, updatenu = updatenu)
    loglik[i] = robust.loglik(par = par, x = x, y = y)
  }
  
  i = 3
  while (aitken.acceleration(loglik[(i - 2):i]) > atol & i < max.iter) {
    i = i + 1
    par = EM.step(par = par, x = x, y = y, updatenu = updatenu)
    loglik[i] = robust.loglik(par = par, x = x, y = y)
  }
  
  val = list(loglik = loglik[1:i], atol = atol, max.iter = max.iter, 
             aa = aitken.acceleration(loglik[1:i]), 
             iter = i, par = par, maxloglik = loglik[i])
  return(val)
}

#SEM
mce.step <- function( par=NULL, x=NULL, y=NULL, m=1 )
{
  ## perform S-step
  t2 = as.numeric(y - x %*% par$beta)^2/par$sigma^2
  
  u.hat = sapply(seq_along(y), FUN=function(i)
  {
    gdata = rgamma( m, shape=(par$nu +1)/2, rate=(par$nu + t2[i])/2 )
    mean( gdata )
  })
  
  return( u.hat )
}

m.step <- function( par=NULL, x=NULL, y=NULL, ua.hat=NULL )
{
  # update beta and sigma
  
  Xtw = t( sweep(x, 1, ua.hat, FUN = "*") )
  par$beta = solve(Xtw %*% x) %*% (Xtw %*% y)
  
  return( par )
}

MCEMn <- function( par0=NULL, x=NULL, y=NULL, n=NULL, m=1 )
{
  ## starting at par0 perform n iterations of EM; m is number of expectations to
  ## update.
  
  # initialize the expectations
  par    = par0
  maxp   = par0
  parbar = par0
  maxll  = -Inf
  ibest  = 0
  loglik = numeric(n)
  
  for ( i in 1:length(loglik) ) 
  {
    # E-step
    u.hat = mce.step(par=par, x=x, y=y, m=m)
    
    # M-step
    par = m.step(par=par, x=x, y=y, ua.hat=u.hat)
    loglik[i] = robust.loglik(par=par, x=x, y=y)
    
    if( loglik[i] > maxll ){ maxp=par; maxll=loglik[i]; ibest=i; }
    
    parbar = (unlist(parbar)*i + unlist(par))/(i+1)
  }
  
  val = list(loglik=loglik, max=maxp, parbar=parbar, ibest=ibest)
  
  return(val)
}

MCEMchain <- function( par0=NULL, x=NULL, y=NULL, ndraws=NULL, m=1 )
{
  ## starting at par0 perform n iterations of EM m is number of expectations to
  ## update.
  
  # initialize the expectations
  par    = par0
  loglik = numeric(ndraws)
  chain  = vector("list", length = ndraws)
  
  
  for ( i in seq_along(loglik) ) 
  {
    # E-step
    u.hat = mce.step(par=par, x=x, y=y, m=m)
    
    # M-step
    par = m.step(par=par, x=x, y=y, ua.hat=u.hat)
    
    loglik[i] = robust.loglik(par=par, x=x, y=y)
    chain[[i]] = par$beta
    
  }
  
  chain |> 
    unlist() |> 
    matrix(ncol=ndraws)->
    cdat
  
  return( cdat )
}

#Profile
EM.profile.beta <- function(
  par=NULL, x=NULL, y=NULL, k=NULL, wt=NULL, updatenu=FALSE 
) {
  ## k which beta variables do we hold fixed. 
  ## perform one iteration of EM
  ## wt are weights
  if (is.null(wt)) wt =rep(1, length(y))
  
  t2   = as.numeric( y - x %*% par$beta )^2/par$sigma^2
  ## E-step
  u.hat = (par$nu + 1)/(par$nu + t2 )
  
  Xtw =  t(sweep(x, 1, u.hat*wt, FUN="*"))
  
  par$beta[-k] = solve( Xtw[-k,] %*% x[,-k] ) %*% 
    ( Xtw[-k,] %*% (y - par$beta[k]*x[,k] ) )
  
  # should we update nu?
  if (updatenu) {
    par$nu = optimize(
      q2.fn,  interval=c(1, 1000), u.hat =u.hat, nu.old=par$nu, wt=wt, 
      maximum=TRUE
    )$maximum
    
  }
  
  return(par)    
}

EMProfn <- function( 
  parb=NULL, parns=NULL, k=NULL, x=NULL, y=NULL, n=NULL, updatenu=FALSE 
) {
  ## perform n iterations of EM using the profile likelihood for the kth beta
  
  loglik = numeric(n)
  par = c( list( "beta"=as.matrix(parb) ), parns)
  for ( i in 1:length(loglik) ) {
    par = EM.profile.beta(par=par, x=x, y=y, k=k, updatenu=updatenu)
    loglik[i] = robust.loglik(par=par, x=x, y=y)
  }
  
  val = loglik[ length(loglik) ]

  return(val)
}

#convergence
aitken.acceleration <- function(loglik = NULL)
{
  if (length(loglik) >= 3) {
    m = length(loglik)
    lk_2 = loglik[-c(m - 1, m)]
    lk_1 = loglik[-c(1, m)]
    lk = loglik[-c(1, 2)]
    
    ak_1 = (lk - lk_1)/(lk_1 - lk_2)
    
    linf.k = lk_1 + (lk - lk_1)/(1 - ak_1)
    
    val = linf.k - lk
  } else {
    val = 1
  }
  
  return(val)
}


# Estimation Functions ----------------------------------------------------

#Get best according to all permutations
all.perm = function(betas=NULL, xdat=NULL, ydat=NULL, sig=NULL, tnu=NULL)
{
  
  val = sapply( 
    seq_along(betas[1,]), FUN=function(ib1)
    {
      sapply(
        seq_along(betas[2,]), FUN=function(ib2)
        {
            rb.ll.betas(
              beta1=betas[1,ib1], 
              beta2=betas[2,ib2], 
              xd=xdat, yd=ydat, sigma=sig, nu=tnu)
        }
      )

    }            
  )
  
}

all.perms.p = function(bchain=NULL, mle=NULL)
{
  sweep( bchain, 1, as.numeric(mle) ) |> 
    (\(u){  u^2  })() |>                  #squared differences
    rowMins() |> 
    (\(u){ cbind(seq_along(u), u) })() |> #grab correct vals from each row
    (\(u){ bchain[u]  })() ->
    val
  
  val
}

# Average of top n terms according to likelihood
topavg = function(betas=NULL, loglik=NULL, nterms=NULL)
{
  ibest = sort(loglik, decreasing=TRUE, index=TRUE)$ix[1:nterms]
  if(nterms==1)
  {
    val = betas[,ibest]
  }
  else
  {
    val = rowmeans(betas[,ibest])
  }
  
  val
}

#Best according to maximum profile likelihood in each component
profll.beta = function(betas=NULL, xdat=NULL, ydat=NULL, par.ns=NULL)
{
  llvals = numeric( length(betas) )
  for( k in seq_along(betas) )
  {
    llvals[k]=EMProfn(parb=betas, parns=par.ns, k=k, x=xdat, y=ydat, n=30)
  }
  llvals
}

profll.chain = function(chain=NULL, x=NULL, y=NULL, nusig=NULL)
{
  chainplls = apply( chain, 2, profll.beta, xdat=x, ydat=y, par.ns=nusig )
  ibest     = rowMaxs( chainplls )
  ibest     = cbind( seq_along(ibest), ibest )
  val       = chain[ ibest ]
}

get.profll.chain = function(chain=NULL, x=NULL, y=NULL, nusig=NULL)
{
  chainplls = apply( chain, 2, profll.beta, xdat=x, ydat=y, par.ns=nusig )
}

#mc integrated
int.ll = function(beta=NULL, bchain=NULL, idx=NULL, llfun=NULL)
{
  bchain[idx,] = beta
  apply(bchain, 2, llfun) |> mean() -> val
  val
}

max.int.ll = function(row=NULL, chain=NULL, index=NULL, lfun=NULL)
{
  sapply( 
    row, FUN=int.ll, bchain=chain, idx=index, llfun=lfun 
  ) |> 
    which.max() -> 
    maxdex
    
  val = c(index, maxdex)
}

int.ll.est = function(bchain=NULL, llfun=NULL)
{
  nrows = length( bchain[,1] )
  ivals = vector("list", nrows)
  
  for( i in seq_len(nrows) )
  {
    ivals[[i]]=max.int.ll( bchain[i,],  chain=bchain, index=i, lfun=llfun )
  }
  ivals |> unlist() |> matrix(nrow=(nrows), byrow=TRUE)->ibestll
  
  bchain[ ibestll ]
}


# Data Generation ---------------------------------------------------------


genX = function(nobs=NULL, p=NULL)
{
  X = mvtnorm::rmvnorm( nobs, rep(0,p), sigma=diag(p) )
  X = cbind(rep(1,nobs), X)
}

genB = function(p=NULL)
{
  beta= runif(p+1, -2, 2)
}

genY = function(xmat=NULL, b=NULL, nusig=NULL)
{
  xmat %*%b + nusig$sigma * rt( dim(xmat)[1], df=nusig$nu  )
}

genRBdata = function(nobs=NULL, p=NULL, nusig=NULL, seed=90053)
{
  if( is.null(nusig) ) nusig = list("sigma"=1, "nu"=2.5)
  
  set.seed(seed)
  
  X = genX(nobs, p)
  b = genB(p)
  y = genY(X, b, nusig)
  
  val = list(
    "x" = X,
    "y" = y,
    "par"= list("beta"=b, "sigma"=nusig$sigma, "nu"=nusig$nu)
  )
  
}



#