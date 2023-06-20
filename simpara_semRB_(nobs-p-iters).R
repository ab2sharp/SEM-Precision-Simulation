# SEM with Robust Regression in parallel

#load
source("semRB_functions.R")
library(parallel)
library(doParallel)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)  #args for nob, p, iters from cl

#Simulation parameters
set.seed( 90053 )
nobs  = as.numeric( args[1] )
p     = as.numeric( args[2] )  # number of covariates (not including intercept)
iters = as.numeric( args[3] )  # length of SEM chain
reps  = 1000                 # simulation repetitions
nusig = list( "sigma" = 2, 
              "nu"    = 2.5
) #specified by user



# Data generation ---------------------------------------------------------

all.dat = genRBdata( nobs, p, nusig )


# EM on simulated data ----------------------------------------------------

# get mle
ipar = rpar( p, nusig )
em100 = EMn( par0=ipar, x=all.dat$x, y=all.dat$y, n=100 )
mlpars = em100$par
mle    = mlpars[[1]]

# relabel data
xdat = all.dat$x
ydat = all.dat$y


# Run simulation ----------------------------------------------------------

cls = makeCluster(120, type="FORK") #change arguments as appropriate
registerDoParallel(cls)
clusterSetRNGStream(cls)

allests = foreach( i=seq_len(reps) ) %dopar%
  {
    schain = MCEMchain( par0=mlpars, x=xdat, y=ydat, ndraws=(iters+20), m=1 )
    schain = schain[ ,-c(1:20) ]
    signu  = mlpars[ c(2,3) ]
    
    #readable log-likelihood function
    rbll = function(ebeta=NULL)
    {
      rb.ll.betav( beta=ebeta, x=xdat, y=ydat, signu=signu )
    }
    
    llbest = rbll( mle ) #llval of mle
    
    
    #All permutations
    apbeta = all.perms.p( schain, mle )
    
    # Likelihood chain and weights
    llchain = apply( schain, 2, FUN=function(acol)
    {
      rbll( acol )
    })
    
    llchain |> 
      (\(u){  exp(u-mean(u)) / ( length(u) * mean(exp(u-mean(u))) )  })() -> 
      llws
    
    ## min llr and average of top p by llr
    semmax = topavg( schain, llchain, nterms=1 )
    top10  = topavg( schain, llchain, nterms=10 )
    
    # Profile
    profbeta = profll.chain( schain, x=xdat, y=ydat, nusig=signu )
    
    #integrated
    intbeta = int.ll.est(schain, rbll)
    
    #Average
    bm = rowmeans(schain)
    
    #Weighted likelihood
    schain |> t() |> (\(u){ u*llws })() |> t() |> rowsums() -> wavebeta
    
    #Negative loglikelihood ratio values
    lls = llbest - c( rbll(semmax),   rbll(bm),     rbll(top10),
                      rbll(wavebeta), rbll(apbeta),
                      rbll(profbeta), rbll(intbeta)
    )
    
    list(
      "betas"=cbind( semmax,   bm,     top10,
                     wavebeta, apbeta, profbeta,
                     intbeta,  mle
      ),
      "lls"=lls
    )
    
    
  }

stopCluster(cls)


filename = paste( c(paste(args, collapse = ""), "v3.RData"), collapse = "")

save( allests, file=filename )

#