#####################################################################################################################################################
## MCMC Model Comparisons for results section Model Fitting. 
## 
## In this repo folder there is multiple mcmc model fitting procedures for the models tested in the paper. Follow the code laid out here 
## and ensure you save the RData panel to save the MCMC outputs for futther processing the in follow code folders. 
## 
## Code contained here runs MCMC model fitting for 10 seperate sporozoite infection model 
## using different combinations of immune markers and dose-response curves 
##
## The model selected following fitting was Model_M5 whereby antibody titre dose relationship was modelled using a Hill function and avidty 
## an exponential dose response curve
##
## List of model names and their descriptions: ab = titre, av = avidity, Expon = Exponential dose-response curve, Hill = Hill fucntion dose-response
## model_M2 : ab-Expon, av:Expon 
## model_M3 : ab-Expon, av-Hill 
## model_M4 : ab-Hill , av-Hill 
## model_M5 : ab-Hill , av-Expon
## model_m6 : Interaction dose-response 
## model_m7 : ab-Expon, 
## model_m8 : ab-Hill ,
## model_m9 : 	      ,	av-Expon  
## model_m10:	        ,	av-Hill 
## model_m11: base model - no antibody dose-response curves. 
##
## Each model is run with a Robbin-Munro algorithm for adaptive tuning of the proposal distribution 
## 
## Many thanks to Michael White, Institut Pasteur, for guidence and developmental help in designing the MCMC algorithm 
## that was adapted and used here. 
##
##################################################################################################################################################### 

set.seed(123)


##------------------------------Robbins-munro step scaler for setting acceptance rate-------------------------------------------
##  returns adjustment to step_scale
##  step_scale = Scaler for step size
##  mc = iteration number 
##  N_adapt = number of iterations when adjustment halved from value at i=0
##  dd = outcome at iteration mc (the min of the acceptance ratio and 1)
##  0.23 is the desired acceptance probability for dd
rm_scale <- function(step_scale, mc, log_prob){
  
  dd = exp(log_prob)
  if( dd < -30 ){ dd = 0 }
  dd = min( dd, 1 )
  
  rm_temp = ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
  
  out = step_scale*exp(rm_temp)
  
  out = max( out, 0.02 ) # smaller jumps scale by 0.02 
  out = min( out, 2)     # larger jumps scale by 2 
  out
}

##-------------------------------------model_M2 : ab-Expon, av:Expon------------------------------------------------------------ 
## LIKELIHOOD
model_M2 <- function(par_M2){
  n          <- par_M2[1]
  sig_n      <- par_M2[2]
  sig_mu     <- par_M2[3]
  b_nanp     <- par_M2[4]
  b_av       <- par_M2[5]
  
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  ##########################
  ## pfT 
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- exp( -log(2)*(nanp/b_nanp) )*exp( -log(2)*(avidity/b_av) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M2 = cmpfun(model_M2, options=list(optimize=3))

loglike_M2_cs = function(par_M2){
  -model_M2(par_M2)
}
loglike_M2_cs = cmpfun(loglike_M2_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M2 = function( par_M2 ){
  n          <- par_M2[1]
  sig_n      <- par_M2[2]
  sig_mu     <- par_M2[3]
  b_nanp     <- par_M2[4]
  b_av       <- par_M2[5]
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_av 
  
  if( b_av>0.001 && b_av<100 )
  {
    prior_b_av = log(0.001/100)#log(dgamma(b_av, shape=6, rate=0.5))
  }else{
    prior_b_av = -LARGE
  }
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp + prior_b_av 
  
  prior
}
prior_M2 = cmpfun(prior_M2, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par           = matrix(NA, nrow=N_mcmc, ncol=7)
colnames(MCMC_par) = c("n", "sig_n", "sig_mu","b_nanp", "b_av", "loglike", "prior")

## Implement MCMC iterations
par_MC = c(150, 200, 67200, 15867, 10 )       ## (n, sig_n, sig_mu, b_ab, bav)
Sigma_MC = diag( (0.8*par_MC)^2 )            ## Initial guess of covariance of MVN proposal dist
prior_MC   = prior_M2( par_MC )                 
loglike_MC = loglike_M2_cs( par_MC ) + prior_MC

for(mc in 1:N_mcmc){
  par_MCp1 = mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)
  
  prior_MCp1 = prior_M2( par_MCp1 )
  
  if( prior_MCp1 > -0.5*LARGE )
  {
    loglike_MCp1 = loglike_M2_cs( par_MCp1 ) + prior_MCp1 
    
    log_prob = min( loglike_MCp1-loglike_MC, 0 )           
    
    if( log(runif(1)) < log_prob ) 
    {
      par_MC = par_MCp1
      
      loglike_MC = loglike_MCp1
      prior_MC   = prior_MCp1
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par[1:(mc-1),1:5] )
    }
  }
  
  MCMC_par[mc,1:5] = par_MC
  MCMC_par[mc,6]   = loglike_MC
  MCMC_par[mc,7]   = prior_MC
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine posterior distributions  
## remove burn in 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]

par(mfrow=c(2,5))

hist(MCMC_burn[,1], breaks=50, main="n")
hist(MCMC_burn[,2], breaks=50, main="sig_n")
hist(MCMC_burn[,3], breaks=50, main="sig_mu")
hist(MCMC_burn[,4], breaks=50, main="b_nanp")
hist(MCMC_burn[,5], breaks=50, main="b_av")

plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,4], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,5], pch=19, cex=0.25, col="grey", main="b_av")

## Median and 95% Credible intervals 
quantile( MCMC_burn[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

#---------------------------------------- model_M3 : ab-Expon, av-Hill----------------------------------------------------------
## LIKELIHOOD 
model_M3 <- function(par_M3){
  n          <- par_M3[1]
  sig_n      <- par_M3[2]
  sig_mu     <- par_M3[3]
  b_nanp     <- par_M3[4]
  b_av       <- par_M3[5]
  a_av       <- par_M3[6]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  ##########################
  ## pfT 
  
  
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- exp( -log(2)*(nanp/b_nanp) )*( 1/(1+(avidity/b_av)^a_av) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M3 = cmpfun(model_M3, options=list(optimize=3))

loglike_M3_cs = function(par_M3){
  -model_M3(par_M3)
}
loglike_M3_cs = cmpfun(loglike_M3_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M3 = function( par_M3 ){
  n          <- par_M3[1]
  sig_n      <- par_M3[2]
  sig_mu     <- par_M3[3]
  b_nanp     <- par_M3[4]
  b_av       <- par_M3[5]
  a_av       <- par_M3[6]
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_av 
  
  if( b_av>0.001 && b_av<100 )
  {
    prior_b_av = log(0.001/100)
  }else{
    prior_b_av = -LARGE
  }
  
  ######################################
  ##  Uniform prior on a_av ~ U(0,30)
  
  if( a_av>0 && a_av<30 )
  {
    prior_a_av = log(0.001/30)
  }else{
    prior_a_av = -LARGE
  }
  
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp + prior_b_av + prior_a_av
  
  prior
}
prior_M3 = cmpfun(prior_M3, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M3           = matrix(NA, nrow=N_mcmc, ncol=8)
colnames(MCMC_par_M3) = c("n", "sig_n", "sig_mu","b_nanp", "b_av", "a_av", "loglike", "prior")

## Implement MCMC iterations
par_MC3 = c(150, 194, 67200, 15000, 10 , 1 )       ## (n, sig_n, sig_mu, b_nanp, b_av, a_av)
Sigma_MC3 = diag( (0.45*par_MC3)^2 )               ## Initial guess of covariance of MVN proposal dist
prior_MC3   = prior_M3( par_MC3 )                 
loglike_MC3 = loglike_M3_cs( par_MC3 ) + prior_MC3

for(mc in 1:N_mcmc){
  par_MCp3 = mvrnorm(n=1, mu=par_MC3, Sigma=step_scale*Sigma_MC3)
  
  prior_MCp3 = prior_M3( par_MCp3 )
  
  if( prior_MCp3 > -0.5*LARGE )
  {
    loglike_MCp3 = loglike_M3_cs( par_MCp3 ) + prior_MCp3 
    
    
    log_prob3 = min( loglike_MCp3-loglike_MC3, 0 )           
    
    if( log(runif(1)) < log_prob3 ) 
    {
      par_MC3 = par_MCp3
      
      loglike_MC3 = loglike_MCp3
      prior_MC3   = prior_MCp3
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob3)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M3[1:(mc-1),1:6] )
    }
  }
  
  
  MCMC_par_M3[mc,1:6] = par_MC3
  MCMC_par_M3[mc,7]   = loglike_MC3
  MCMC_par_M3[mc,8]   = prior_MC3
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and positerior distributions 
par(mfrow=c(2,6))
MCMC_burn_M3 <- MCMC_par_M3[floor(0.2*nrow(MCMC_par_M3)):(nrow(MCMC_par_M3)-1),]

hist(MCMC_burn_M3[,1], breaks=50, main="n")
hist(MCMC_burn_M3[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M3[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M3[,4], breaks=50, main="b_nanp")
hist(MCMC_burn_M3[,5], breaks=50, main="b_av")
hist(MCMC_burn_M3[,6], breaks=50, main="a_av")

plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,4], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,5], pch=19, cex=0.25, col="grey", main="b_av")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,6], pch=19, cex=0.25, col="grey", main="a_av")

## median parameter values and 95% Credible Intervals 
quantile( MCMC_burn_M3[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,6], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

#------------------------------------model_M4 : ab-Hill, av-Hill-----------------------------------------------------------------
# LIKELIHOOD
model_M4 <- function(par_M4){
  n          <- par_M4[1]
  sig_n      <- par_M4[2]
  sig_mu     <- par_M4[3]
  b_nanp     <- par_M4[4]
  a_nanp     <- par_M4[5]
  b_av       <- par_M4[6]
  a_av       <- par_M4[7]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )*( 1/(1+(avidity/b_av)^a_av) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M4 = cmpfun(model_M4, options=list(optimize=3))

loglike_M4_cs = function(par_M4){
  -model_M4(par_M4)
}
loglike_M4_cs = cmpfun(loglike_M4_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M4 = function( par_M4 ){
  n          <- par_M4[1]
  sig_n      <- par_M4[2]
  sig_mu     <- par_M4[3]
  b_nanp     <- par_M4[4]
  a_nanp     <- par_M4[5]
  b_av       <- par_M4[6]
  a_av       <- par_M4[7]
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Uniform prior on a_nanp ~ U(0,20)
  
  if( a_nanp>0.001 && a_nanp<30 )
  {
    prior_a_nanp = log(0.001/30)
  }else{
    prior_a_nanp = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_av 
  
  if( b_av>0.001 && b_av<100 )
  {
    prior_b_av = log(0.001/100)#log(dgamma(b_av, shape=6, rate=0.5))
  }else{
    prior_b_av = -LARGE
  }
  ######################################
  ## Uniform prior on a_av ~ U(0,30)
  
  if( a_av>0 && a_av<30 )
  {
    prior_a_av = log(0.001/30)
  }else{
    prior_a_av = -LARGE
  }
  
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp + prior_a_nanp + prior_b_av + prior_a_av
  
  prior
}
prior_M4 = cmpfun(prior_M4, options=list(optimize=3))

##  MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 3           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M4           = matrix(NA, nrow=N_mcmc, ncol=9)
colnames(MCMC_par_M4) = c("n", "sig_n", "sig_mu","b_nanp", "a_nanp", "b_av", "a_av", "loglike", "prior")

## Implement MCMC iterations
par_MC4 = c(150, 194, 67200, 15000, 1, 10 , 1 )       ## (n, sig_n, sig_mu, b_nanp, a_nanp, b_av, a_av)
Sigma_MC4 = diag( (0.6*par_MC4)^2 )                  ## Initial guess of covariance of MVN proposal dist
prior_MC4   = prior_M4( par_MC4 )                 
loglike_MC4 = loglike_M4_cs( par_MC4 ) + prior_MC4

for(mc in 1:N_mcmc){
  par_MCp4 = mvrnorm(n=1, mu=par_MC4, Sigma=step_scale*Sigma_MC4)
  
  prior_MCp4 = prior_M4( par_MCp4 )
  
  if( prior_MCp4 > -0.5*LARGE )
  {
    loglike_MCp4 = loglike_M4_cs( par_MCp4 ) + prior_MCp4 
    
    
    log_prob4 = min( loglike_MCp4-loglike_MC4, 0 )           
    
    if( log(runif(1)) < log_prob4 ) 
    {
      par_MC4 = par_MCp4
      
      loglike_MC4 = loglike_MCp4
      prior_MC4   = prior_MCp4
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob4)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M4[1:(mc-1),1:7] )
    }
  }
  
  
  MCMC_par_M4[mc,1:7] = par_MC4
  MCMC_par_M4[mc,8]   = loglike_MC4
  MCMC_par_M4[mc,9]   = prior_MC4
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distributions 
par(mfrow=c(2,7))

MCMC_burn_M4 <- MCMC_par_M4[floor(0.2*nrow(MCMC_par_M4)):(nrow(MCMC_par_M4)-1),]

hist(MCMC_burn_M4[,1], breaks=50, main="n")
hist(MCMC_burn_M4[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M4[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M4[,4], breaks=50, main="b_nanp")
hist(MCMC_burn_M4[,5], breaks=50, main="a_nanp")
hist(MCMC_burn_M4[,6], breaks=50, main="b_av")
hist(MCMC_burn_M4[,7], breaks=50, main="a_av")

plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,4], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,5], pch=19, cex=0.25, col="grey", main="a_nanp")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,6], pch=19, cex=0.25, col="grey", main="b_av")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,7], pch=19, cex=0.25, col="grey", main="a_av")

## median parameter values and 95% Credible Intervals
quantile( MCMC_burn_M4[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,6], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,7], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

##-----------------------------model_M5 ab-Hill , av-Expon-------------------------------------------------------------------
# LIKELIHOOD
model_M5 <- function(par_M5){
  n          <- par_M5[1]
  sig_n      <- par_M5[2]
  sig_mu     <- par_M5[3]
  b_nanp     <- par_M5[4]
  a_nanp     <- par_M5[5]
  b_av       <- par_M5[6]
  
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )* exp( -log(2)*(avidity/b_av) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M5 = cmpfun(model_M5, options=list(optimize=3))

loglike_M5_cs = function(par_M5){
  -model_M5(par_M5)
}
loglike_M5_cs = cmpfun(loglike_M5_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M5 = function( par_M5 ){
  n          <- par_M5[1]
  sig_n      <- par_M5[2]
  sig_mu     <- par_M5[3]
  b_nanp     <- par_M5[4]
  a_nanp     <- par_M5[5]
  b_av       <- par_M5[6]
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Uniform prior on a_nanp ~ U(0,20)
  
  if( a_nanp>0.001 && a_nanp<30 )
  {
    prior_a_nanp = log(0.001/30)
  }else{
    prior_a_nanp = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_av 
  
  if( b_av>0.001 && b_av<100 )
  {
    prior_b_av = log(0.001/100)#log(dgamma(b_av, shape=6, rate=0.5))
  }else{
    prior_b_av = -LARGE
  }
  
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp + prior_a_nanp + prior_b_av 
  
  prior
}
prior_M5 = cmpfun(prior_M5, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 
step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M5           = matrix(NA, nrow=N_mcmc, ncol=8)
colnames(MCMC_par_M5) = c("n", "sig_n", "sig_mu","b_nanp", "a_nanp", "b_av",  "loglike", "prior")

## Implement MCMC iterations
par_MC5 = c(150, 194, 67200, 15000, 1, 10  )       ## (n, sig_n, sig_mu, b_nanp, a_nanp, b_av )
Sigma_MC5 = diag( (0.25*par_MC5)^2 )               ## Initial guess of covariance of MVN proposal dist
prior_MC5   = prior_M5( par_MC5 )                 
loglike_MC5 = loglike_M5_cs( par_MC5 ) + prior_MC5

for(mc in 1:N_mcmc){
  par_MCp5 = mvrnorm(n=1, mu=par_MC5, Sigma=step_scale*Sigma_MC5)
  
  prior_MCp5 = prior_M5( par_MCp5 )
  
  if( prior_MCp5 > -0.5*LARGE )
  {
    loglike_MCp5 = loglike_M5_cs( par_MCp5 ) + prior_MCp5 
    
    
    log_prob5 = min( loglike_MCp5-loglike_MC5, 0 )           
    
    if( log(runif(1)) < log_prob5 ) 
    {
      par_MC5 = par_MCp5
      
      loglike_MC5 = loglike_MCp5
      prior_MC5  = prior_MCp5
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob5)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M5[1:(mc-1),1:6] )
    }
  }
  
  MCMC_par_M5[mc,1:6] = par_MC5
  MCMC_par_M5[mc,7]   = loglike_MC5
  MCMC_par_M5[mc,8]   = prior_MC5
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distributions 
MCMC_burn_M5 <- MCMC_par_M5[floor(0.2*nrow(MCMC_par_M5)):(nrow(MCMC_par_M5)-1),]

#png("Fig_3.png", width = 600, height = 400, units = "px")#width = 80, height = 40, units = "mm", res=600)#, width=800, height=400, res=600)
par(oma=c(1,1,1,1)) 
par(mar=c(4,1,1,1) + 0.1)
par(mfrow=c(2,6))
options(scipen=10)

hist(MCMC_burn_M5[,1], breaks=50, main=" ", xlab=expression("n"), ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M5[,2], breaks=50, xlab=expression(sigma["n"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M5[,3], breaks=50, xlab=expression(sigma[mu]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M5[,4], breaks=50, xlab=expression(beta["t"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M5[,5], breaks=50, xlab=expression(alpha["t"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M5[,6], breaks=50, xlab=expression(beta["ai"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)

plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,1], pch=19, cex=0.25, col="grey", main=" ", ylab = expression("n"), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,2], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(sigma["n"]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,3], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(sigma[mu]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,4], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(beta["t"]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,5], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(alpha["t"]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,6], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(beta["ai"]), xlab="Iteration")

#dev.off()

### Median parameter values and 95% Credible Interval 
quantile( MCMC_burn_M5[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,6], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

##------------------------------------------model_m6 : Interaction dose-response--------------------------------------------------
## LIKELIHOOD
model_M6 <- function(par_M6){
  n         <- par_M6[1]
  sig_n     <- par_M6[2]
  sig_mu    <- par_M6[3]
  b_nanp    <- par_M6[4]
  b_av      <- par_M6[5]
  gamma     <- par_M6[6]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- exp( -0.5*log(2)*( (nanp/b_nanp+avidity/b_av) +
                             sqrt( (nanp/b_nanp+avidity/b_av)^2 + 4*gamma*(nanp/b_nanp)*(avidity/b_av) ) ) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M6 = cmpfun(model_M6, options=list(optimize=3))

loglike_M6_cs = function(par_M6){
  -model_M6(par_M6)
}
loglike_M6_cs = cmpfun(loglike_M6_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M6 = function( par_M6 ){
  n         <- par_M6[1]
  sig_n     <- par_M6[2]
  sig_mu    <- par_M6[3]
  b_nanp    <- par_M6[4]
  b_av      <- par_M6[5]
  gamma     <- par_M6[6]
 
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_av 
  
  if( b_av>0.001 && b_av<100 )
  {
    prior_b_av = log(0.001/100)#log(dgamma(b_av, shape=6, rate=0.5))
  }else{
    prior_b_av = -LARGE
  }
  
  ######################################
  ## Normal prior on gamma ~ limit at -1 : N(1,100)
  
  if( gamma>-1 && gamma<100 )
  {
    prior_gamma = log(dnorm(x = gamma, mean = 0, sd=10))
  }else{
    prior_gamma = -LARGE
  }
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp +  prior_b_av + prior_gamma
  
  prior
}
prior_M6 = cmpfun(prior_M6, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 3           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M6           = matrix(NA, nrow=N_mcmc, ncol=8)
colnames(MCMC_par_M6) = c("n", "sig_n", "sig_mu","b_nanp", "b_av", "gamma", "loglike", "prior")

## Implement MCMC iterations
par_MC6 = c(150, 194, 67200, 15000,  10 , 10  )       ## (n, sig_n, sig_mu, b_nanp, b_av, gamma)
Sigma_MC6 = diag( (0.8*par_MC6)^2 )      ## Initial guess of covariance of MVN proposal dist
prior_MC6   = prior_M6( par_MC6 )                 
loglike_MC6 = loglike_M6_cs( par_MC6 ) + prior_MC6

for(mc in 1:N_mcmc){
  par_MCp6 = mvrnorm(n=1, mu=par_MC6, Sigma=step_scale*Sigma_MC6)
  
  prior_MCp6 = prior_M6( par_MCp6 )
  
  if( prior_MCp6 > -0.5*LARGE )
  {
    loglike_MCp6 = loglike_M6_cs( par_MCp6 ) + prior_MCp6 
    
    
    log_prob6 = min( loglike_MCp6-loglike_MC6, 0 )           
    
    if( log(runif(1)) < log_prob6 ) 
    {
      par_MC6 = par_MCp6
      
      loglike_MC6 = loglike_MCp6
      prior_MC6  = prior_MCp6
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob6)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M6[1:(mc-1),1:6] )
    }
  }
  
  
  MCMC_par_M6[mc,1:6] = par_MC6
  MCMC_par_M6[mc,7]   = loglike_MC6
  MCMC_par_M6[mc,8]   = prior_MC6
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distrbutions  
MCMC_burn_M6 <- MCMC_par_M6[floor(0.2*nrow(MCMC_par_M6)):(nrow(MCMC_par_M6)-1),]

par(mfrow=c(2,6))

hist(MCMC_burn_M6[,1], breaks=50, main="n")
hist(MCMC_burn_M6[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M6[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M6[,4], breaks=50, main="b_nanp")
hist(MCMC_burn_M6[,5], breaks=50, main="b_av")
hist(MCMC_burn_M6[,6], breaks=50, main="gamma")

plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,4], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,5], pch=19, cex=0.25, col="grey", main="b_av")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,6], pch=19, cex=0.25, col="grey", main="gamma")

## Median paramter values and 95% Credible Intervals 
quantile( MCMC_burn_M6[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,6], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )


##------------------------------------------model_m7 : ab-Expon-----------------------------------------------------------------
## LIKELIHOOD
model_M7 <- function(par_M7){
  n          <- par_M7[1]
  sig_n      <- par_M7[2]
  sig_mu     <- par_M7[3]
  b_nanp     <- par_M7[4]
  
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- exp( -log(2)*(nanp/b_nanp) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M7 = cmpfun(model_M7, options=list(optimize=3))

loglike_M7_cs = function(par_M7){
  -model_M7(par_M7)
}
loglike_M7_cs = cmpfun(loglike_M7_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M7 = function( par_M7 ){
  n          <- par_M7[1]
  sig_n      <- par_M7[2]
  sig_mu     <- par_M7[3]
  b_nanp     <- par_M7[4]
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp 
  
  prior
}
prior_M7 = cmpfun(prior_M7, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M7           = matrix(NA, nrow=N_mcmc, ncol=6)
colnames(MCMC_par_M7) = c("n", "sig_n", "sig_mu","b_nanp",  "loglike", "prior")

## Implement MCMC iterations
par_MC7 = c(150, 194, 67200, 15000 )       ## (n, sig_n, sig_mu, b_nanp)
Sigma_MC7 = diag( (0.25*par_MC7)^2 )       ## Initial guess of covariance of MVN proposal dist
prior_MC7   = prior_M7( par_MC7 )                 
loglike_MC7 = loglike_M7_cs( par_MC7 ) + prior_MC7

for(mc in 1:N_mcmc){
  par_MCp7 = mvrnorm(n=1, mu=par_MC7, Sigma=step_scale*Sigma_MC7)
  
  prior_MCp7 = prior_M7( par_MCp7 )
  
  if( prior_MCp7 > -0.5*LARGE )
  {
    loglike_MCp7 = loglike_M7_cs( par_MCp7 ) + prior_MCp7 
    
    
    log_prob7 = min( loglike_MCp7-loglike_MC7, 0 )           
    
    if( log(runif(1)) < log_prob7 ) 
    {
      par_MC7 = par_MCp7
      
      loglike_MC7 = loglike_MCp7
      prior_MC7   = prior_MCp7
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob7)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M7[1:(mc-1),1:4] )
    }
  }
  
  
  MCMC_par_M7[mc,1:4] = par_MC7
  MCMC_par_M7[mc,5]   = loglike_MC7
  MCMC_par_M7[mc,6]   = prior_MC7
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distributions  
MCMC_burn_M7 <- MCMC_par_M7[floor(0.2*nrow(MCMC_par_M7)):(nrow(MCMC_par_M7)-1),]

par(mfrow=c(2,4))

hist(MCMC_burn_M7[,1], breaks=50, main="n")
hist(MCMC_burn_M7[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M7[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M7[,4], breaks=50, main="b_nanp")

plot(x=1:nrow(MCMC_burn_M7), y=MCMC_burn_M7[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M7), y=MCMC_burn_M7[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M7), y=MCMC_burn_M7[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M7), y=MCMC_burn_M7[,4], pch=19, cex=0.25, col="grey", main="b_nanp")

## Median parameter values and 95% Credible Intervals 
quantile( MCMC_burn_M7[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M7[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M7[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M7[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

##---------------------------------------------model_M8 : ab-Hill--------------------------------------------------------------
## LIKELIHOOD
model_M8 <- function(par_M8){
  n          <- par_M8[1]
  sig_n      <- par_M8[2]
  sig_mu     <- par_M8[3]
  b_nanp     <- par_M8[4]
  a_nanp     <- par_M8[5]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M8 = cmpfun(model_M8, options=list(optimize=3))

loglike_M8_cs = function(par_M8){
  -model_M8(par_M8)
}
loglike_M8_cs = cmpfun(loglike_M8_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M8 = function( par_M8 ){
  n          <- par_M8[1]
  sig_n      <- par_M8[2]
  sig_mu     <- par_M8[3]
  b_nanp     <- par_M8[4]
  a_nanp     <- par_M8[5]
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ## Uniform prior on a_nanp ~ U(0,30)
  
  if( a_nanp>0 && a_nanp<30 )
  {
    prior_a_nanp = log(0.001/30)
  }else{
    prior_a_nanp = -LARGE
  }
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_nanp + prior_a_nanp
  
  prior
}
prior_M8 = cmpfun(prior_M8, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M8           = matrix(NA, nrow=N_mcmc, ncol=7)
colnames(MCMC_par_M8) = c("n", "sig_n", "sig_mu","b_nanp", "a_nanp", "loglike", "prior")

## Implement MCMC iterations
par_MC8 = c(150, 194, 67200, 15000, 1)       ## (n, sig_n, sig_mu, b_nanp, a_nanp)
Sigma_MC8 = diag( (0.5*par_MC8)^2 )      ## Initial guess of covariance of MVN proposal dist
prior_MC8   = prior_M8( par_MC8 )                 
loglike_MC8 = loglike_M8_cs( par_MC8 ) + prior_MC8

for(mc in 1:N_mcmc){
  par_MCp8 = mvrnorm(n=1, mu=par_MC8, Sigma=step_scale*Sigma_MC8)
  
  prior_MCp8 = prior_M8( par_MCp8 )
  
  if( prior_MCp8 > -0.5*LARGE )
  {
    loglike_MCp8 = loglike_M8_cs( par_MCp8 ) + prior_MCp8 
    
    
    log_prob8 = min( loglike_MCp8-loglike_MC8, 0 )           
    
    if( log(runif(1)) < log_prob8 ) 
    {
      par_MC8 = par_MCp8
      
      loglike_MC8 = loglike_MCp8
      prior_MC8  = prior_MCp8
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob8)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M8[1:(mc-1),1:5] )
    }
  }
  
  
  MCMC_par_M8[mc,1:5] = par_MC8
  MCMC_par_M8[mc,6]   = loglike_MC8
  MCMC_par_M8[mc,7]   = prior_MC8
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distributions  
MCMC_burn_M8 <- MCMC_par_M8[floor(0.2*nrow(MCMC_par_M8)):(nrow(MCMC_par_M8)-1),]

par(mfrow=c(2,5))

hist(MCMC_burn_M8[,1], breaks=50, main="n")
hist(MCMC_burn_M8[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M8[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M8[,4], breaks=50, main="b_nanp")
hist(MCMC_burn_M8[,5], breaks=50, main="a_nanp")

plot(x=1:nrow(MCMC_burn_M8), y=MCMC_burn_M8[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M8), y=MCMC_burn_M8[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M8), y=MCMC_burn_M8[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M8), y=MCMC_burn_M8[,4], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M8), y=MCMC_burn_M8[,5], pch=19, cex=0.25, col="grey", main="a_nanp")

## Median parameter values and 95% Credible Interval 
quantile( MCMC_burn_M8[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M8[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M8[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M8[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M8[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

##--------------------------------------model_m9 : ,	av-Expon----------------------------------------------------------------- 
## LIKELIHOOD
model_M9 <- function(par_M9){
  n          <- par_M9[1]
  sig_n      <- par_M9[2]
  sig_mu     <- par_M9[3]
  b_av       <- par_M9[4]
  
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- exp( -log(2)*(avidity/b_av) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M9 = cmpfun(model_M9, options=list(optimize=3))

loglike_M9_cs = function(par_M9){
  -model_M9(par_M9)
}
loglike_M9_cs = cmpfun(loglike_M9_cs, options=list(optimize=3)) 

## PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M9 = function( par_M9 ){
  n          <- par_M9[1]
  sig_n      <- par_M9[2]
  sig_mu     <- par_M9[3]
  b_av       <- par_M9[4]
  
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }

  ######################################
  ##  gamma prior on b_av 
  
  if( b_av>1 && b_av<1000 )
  {
    prior_b_av = log(0.001/100)
  }else{
    prior_b_av = -LARGE
  }
  
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_av 
  
  prior
}
prior_M9 = cmpfun(prior_M9, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M9           = matrix(NA, nrow=N_mcmc, ncol=6)
colnames(MCMC_par_M9) = c("n", "sig_n", "sig_mu", "b_av",  "loglike", "prior")

## Implement MCMC iterations
par_MC9 = c(150, 194, 67200, 10  )       ## (n, sig_n, sig_mu, b_av)
Sigma_MC9 = diag( (0.25*par_MC9)^2 )      ## Initial guess of covariance of MVN proposal dist
prior_MC9   = prior_M9( par_MC9 )                 
loglike_MC9 = loglike_M9_cs( par_MC9 ) + prior_MC9

for(mc in 1:N_mcmc){
  par_MCp9 = mvrnorm(n=1, mu=par_MC9, Sigma=step_scale*Sigma_MC9)
  
  prior_MCp9 = prior_M9( par_MCp9 )
  
  if( prior_MCp9 > -0.5*LARGE )
  {
    loglike_MCp9 = loglike_M9_cs( par_MCp9 ) + prior_MCp9 
    
    
    log_prob9 = min( loglike_MCp9-loglike_MC9, 0 )           
    
    if( log(runif(1)) < log_prob9 ) 
    {
      par_MC9 = par_MCp9
      
      loglike_MC9 = loglike_MCp9
      prior_MC9  = prior_MCp9
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob9)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M9[1:(mc-1),1:4] )
    }
  }
  
  MCMC_par_M9[mc,1:4] = par_MC9
  MCMC_par_M9[mc,5]   = loglike_MC9
  MCMC_par_M9[mc,6]   = prior_MC9
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and Posterior Distribution  
MCMC_burn_M9 <- MCMC_par_M9[floor(0.2*nrow(MCMC_par_M9)):(nrow(MCMC_par_M9)-1),]

par(mfrow=c(2,4))

hist(MCMC_burn_M9[,1], breaks=50, main="n")
hist(MCMC_burn_M9[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M9[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M9[,4], breaks=50, main="b_av")

plot(x=1:nrow(MCMC_burn_M9), y=MCMC_burn_M9[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M9), y=MCMC_burn_M9[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M9), y=MCMC_burn_M9[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M9), y=MCMC_burn_M9[,4], pch=19, cex=0.25, col="grey", main="b_av")

## Median parameter values and 95% Credible Intervals 
quantile( MCMC_burn_M9[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M9[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M9[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M9[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )


##---------------------------------------------------model_m10: ,	av-Hill--------------------------------------------------------
## LIKELIHOOD
model_M10 <- function(par_M10){
  n          <- par_M10[1]
  sig_n      <- par_M10[2]
  sig_mu     <- par_M10[3]
  b_av       <- par_M10[4]
  a_av       <- par_M10[5]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- ( 1/(1+(avidity/b_av)^a_av) )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M10 = cmpfun(model_M10, options=list(optimize=3))

## LIKELIHOOD
loglike_M10_cs = function(par_M10){
  -model_M10(par_M10)
}
loglike_M10_cs = cmpfun(loglike_M10_cs, options=list(optimize=3)) 

## 1.3 PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M10 = function( par_M10 ){
  n          <- par_M10[1]
  sig_n      <- par_M10[2]
  sig_mu     <- par_M10[3]
  b_av       <- par_M10[4]
  a_av       <- par_M10[5]
  
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  ######################################
  ##  guniform prior on b_av 
  
  if( b_av>1 && b_av<100 )
  {
    prior_b_av = log(0.001/100)
  }else{
    prior_b_av = -LARGE
  }
  ######################################
  ##  uniform prior on a_av ~ U(0,30)
  
  if( a_av>0 && a_av<30 )
  {
    prior_a_av = log(0.001/30)
  }else{
    prior_a_av = -LARGE
  }
  
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_b_av + prior_a_av
  
  prior
}
prior_M10 = cmpfun(prior_M10, options=list(optimize=3))

##  MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M10           = matrix(NA, nrow=N_mcmc, ncol=7)
colnames(MCMC_par_M10) = c("n", "sig_n", "sig_mu", "b_av", "a_av", "loglike", "prior")

## 2.3 Implement MCMC iterations
par_MC10 = c(150, 194, 67200, 10 , 1 )       ## (n, sig_n, sig_mu, )
Sigma_MC10 = diag( (0.8*par_MC10)^2 )       ## Initial guess of covariance of MVN proposal dist
prior_MC10   = prior_M10( par_MC10 )                 
loglike_MC10 = loglike_M10_cs( par_MC10 ) + prior_MC10

for(mc in 1:N_mcmc){
  par_MCp10 = mvrnorm(n=1, mu=par_MC10, Sigma=step_scale*Sigma_MC10)
  
  prior_MCp10 = prior_M10( par_MCp10 )
  
  if( prior_MCp10 > -0.5*LARGE )
  {
    loglike_MCp10 = loglike_M10_cs( par_MCp10 ) + prior_MCp10 
    
    
    log_prob10 = min( loglike_MCp10-loglike_MC10, 0 )           
    
    if( log(runif(1)) < log_prob10 ) 
    {
      par_MC10 = par_MCp10
      
      loglike_MC10 = loglike_MCp10
      prior_MC10   = prior_MCp10
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob10)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M10[1:(mc-1),1:5] )
    }
  }
  
  MCMC_par_M10[mc,1:5] = par_MC10
  MCMC_par_M10[mc,6]   = loglike_MC10
  MCMC_par_M10[mc,7]   = prior_MC10
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and Posterior Distributions
MCMC_burn_M10 <- MCMC_par_M10[floor(0.2*nrow(MCMC_par_M10)):(nrow(MCMC_par_M10)-1),]

par(mfrow=c(2,5))

hist(MCMC_burn_M10[,1], breaks=50, main="n")
hist(MCMC_burn_M10[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M10[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M10[,4], breaks=50, main="b_av")
hist(MCMC_burn_M10[,5], breaks=50, main="a_av")

plot(x=1:nrow(MCMC_burn_M10), y=MCMC_burn_M10[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M10), y=MCMC_burn_M10[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M10), y=MCMC_burn_M10[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M10), y=MCMC_burn_M10[,4], pch=19, cex=0.25, col="grey", main="b_av")
plot(x=1:nrow(MCMC_burn_M10), y=MCMC_burn_M10[,5], pch=19, cex=0.25, col="grey", main="a_av")

quantile( MCMC_burn_M10[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M10[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M10[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M10[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M10[,5], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

##------------------------------model_m11: base model - no antibody dose-response curves----------------------------------------
model_M11 <- function(par_M11){
  n        <- par_M11[1]
  sig_n    <- par_M11[2]
  sig_mu   <- par_M11[3]
  VE       <- par_M11[4]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  
  ############################
  ## Secondary Gamma parameters
  
  theta_g <- sig_mu*sig_mu/mu
  
  ##############################
  ## Dose response
  
  DR <- rep( 1-VE, NN )
  DR[which(DR<min_DR)] <- min_DR
  p_spz = n*DR/(n*DR + r)
  
  II <- 1:max_spz
  
  logL <- 0
  # Sporozoite infection model and likelihood equation  
  
  logL <- 0
  for(j in 1:NN){
    if( infect[j]==1 ){
      coeffs <- lgamma(II+r) - lgamma(II+1) - lgamma(r) + r*log(1-p_spz[j]) + II*log(p_spz[j])
      coeffs <- exp(coeffs)
      
      logL <- logL + log(sum( coeffs*dgamma(Q[j], shape=mu*II/theta_g, scale=theta_g) ))
    }
    if( infect[j]==0 ){
      logL <- logL + r*log(1-p_spz[j])
    }
  }
  
  -logL
}
model_M11 = cmpfun(model_M11, options=list(optimize=3))

## LIKELIHOOD
loglike_M11_cs = function(par_M11){
  -model_M11(par_M11)
}
loglike_M11_cs = cmpfun(loglike_M11_cs, options=list(optimize=3)) 

##  PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M11 = function( par_M11 ){
  n          <- par_M11[1]
  sig_n      <- par_M11[2]
  sig_mu     <- par_M11[3]
  VE         <- par_M11[4]
  
  ######################################
  ## Gamma prior on n 
  
  if( n>0 && n<500 )
  {
    prior_n = log(dgamma(n, shape=7.5, rate=0.05))
  }else{
    prior_n = -LARGE
  }
  
  ######################################
  ## Gamma prior on sig_n 
  
  if( sig_n>20 && sig_n<1000 )
  {
    prior_sig_n = log(20/1000)
  }else{
    prior_sig_n = -LARGE
  }
  
  
  ######################################
  ## Gamma prior on sig_mu 
  
  if( sig_mu>30000 && sig_mu<106800 )
  {
    prior_sig_mu = log(dgamma(sig_mu, shape=31.35, rate=0.0005))
  }else{
    prior_sig_mu = -LARGE
  }
  
  
  ######################################
  ##  Gamma prior on VE 
  
  if( VE >0 && VE<0.99 )
  {
    prior_VE = log(dgamma(VE, shape=5, rate=7))
  }else{
    prior_VE = -LARGE
  }
  
  prior = prior_n + prior_sig_n  + prior_sig_mu + prior_VE
  
  prior
}
prior_M11 = cmpfun(prior_M11, options=list(optimize=3))

## MCMC 
N_mcmc       = 200000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 5000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 6000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M11           = matrix(NA, nrow=N_mcmc, ncol=6)
colnames(MCMC_par_M11) = c("n", "sig_n", "sig_mu", "VE", "loglike", "prior")

##  Implement MCMC iterations
par_MC11 = c(150, 194, 67200, 0.5 )       ## (n, sig_n, sig_mu, VE)
Sigma_MC11 = diag( (0.95*par_MC11)^2 )      ## Initial guess of covariance of MVN proposal dist
prior_MC11   = prior_M11( par_MC11 )                 
loglike_MC11 = loglike_M11_cs( par_MC11 ) + prior_MC11

for(mc in 1:N_mcmc){
  par_MCp11 = mvrnorm(n=1, mu=par_MC11, Sigma=step_scale*Sigma_MC11)
  
  prior_MCp11 = prior_M11( par_MCp11 )
  
  if( prior_MCp11 > -0.5*LARGE )
  {
    loglike_MCp11 = loglike_M11_cs( par_MCp11 ) + prior_MCp11 
    
    
    log_prob11 = min( loglike_MCp11-loglike_MC11, 0 )           
    
    if( log(runif(1)) < log_prob11 ) 
    {
      par_MC11 = par_MCp11
      
      loglike_MC11 = loglike_MCp11
      prior_MC11   = prior_MCp11
      
      MCMC_accept = MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt )
    {
      step_scale = rm_scale( step_scale, mc, log_prob11)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC = cov( MCMC_par_M11[1:(mc-1),1:4] )
    }
  }
  
  
  MCMC_par_M11[mc,1:4] = par_MC11
  MCMC_par_M11[mc,5]   = loglike_MC11
  MCMC_par_M11[mc,6]   = prior_MC11
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and Posterior Distributions 
MCMC_burn_M11 <- MCMC_par_M11[floor(0.2*nrow(MCMC_par_M11)):(nrow(MCMC_par_M11)-1),]

par(mfrow=c(2,4))

hist(MCMC_burn_M11[,1], breaks=50, main="n")
hist(MCMC_burn_M11[,2], breaks=50, main="sig_n")
hist(MCMC_burn_M11[,3], breaks=50, main="sig_mu")
hist(MCMC_burn_M11[,4], breaks=50, main="VE")

plot(x=1:nrow(MCMC_burn_M11), y=MCMC_burn_M11[,1], pch=19, cex=0.25, col="grey", main="n")
plot(x=1:nrow(MCMC_burn_M11), y=MCMC_burn_M11[,2], pch=19, cex=0.25, col="grey", main="sig_n")
plot(x=1:nrow(MCMC_burn_M11), y=MCMC_burn_M11[,3], pch=19, cex=0.25, col="grey", main="sig_mu")
plot(x=1:nrow(MCMC_burn_M11), y=MCMC_burn_M11[,4], pch=19, cex=0.25, col="grey", main="VE")

## Median Parameter Values and 95% Credible Intervals 
quantile( MCMC_burn_M11[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M11[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M11[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M11[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

#-------------------------------------------SAVE DATA---------------------------------------------------------------------------
save.image("Thesis_march_update.RData")


