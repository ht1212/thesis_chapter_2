##
## BINARY INFECTION MODEL 
## 
## 
##
## 30-07-2020 
##---------------------------------------
## List of model names and their descriptions: 
## ab = titre
## av = avidity 
## Expon = Exponential dose-response curve
## Hill = Hill fucntion dose-response
##
## model_M2 : ab-Expon, av-Expon 
## model_M3 : ab-Expon, av-Hill 
## model_M4 : ab-Hill , av-Hill 
## model_M5 : ab-Hill , av-Expon
##
## Each model is run with a Robbin-Munro algorithm for adaptive tuning of the proposal distribution 

set.seed(123)

#--------------------Robbins-munro step scaler function------------------------------------------------------------------------
rm_scale <- function(step_scale, mc, log_prob){
  
  dd = exp(log_prob)
  if( dd < -30 ){ dd = 0 }
  dd = min( dd, 1 )
  
  rm_temp = ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
  
  out = step_scale*exp(rm_temp)
  
  out = max( out, 0.02 )
  out = min( out, 2)
  out
}

#--------------------model_M2 : ab-Expon, av:Expon-----------------------------------------------------------------------------
# MODEL LIKELIHOOD 
model_M2 <- function(par_M2){
  b_nanp     <- par_M2[1]
  b_av       <- par_M2[2]
  
  
  ##############################
  ## Dose response
  
  DR <- exp( -log(2)*(nanp/b_nanp) )*exp( -log(2)*(avidity/b_av) )
  DR[which(DR<min_DR)] <- min_DR

  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      
      logL <- logL + sum(log(DR))
    }
    
    if( infect[j]==0 ){
      logL <- logL + sum(log(1-DR))
    }
  }
  
  -logL
}
model_M2 = cmpfun(model_M2, options=list(optimize=3))

loglike_M2_cs = function(par_M2){
  -model_M2(par_M2)
}
loglike_M2_cs = cmpfun(loglike_M2_cs, options=list(optimize=3)) 

# PRIOR
LARGE = 1e10     ## Large value for rejecting parameters with prior
prior_M2 = function( par_M2 ){
  b_nanp     <- par_M2[1]
  b_av       <- par_M2[2]
  
  ######################################
  ## Gamma prior on b_nanp 
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
    # prior_b_nanp = log(500/150000)
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Uniform prior on b_av U(1,100)
  if( b_av>1 && b_av<100 )
  {
   # prior_b_av = log(dgamma(b_av, shape=6, rate=0.5))
    prior_b_av = log(1/100)
  }else{
    prior_b_av = -LARGE
  }
  
  prior =  prior_b_nanp + prior_b_av 
  
  prior
}
prior_M2 = cmpfun(prior_M2, options=list(optimize=3))

##  MCMC 
N_mcmc       = 800000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 10000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 11000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1             ## Scaler for step size
MCMC_accept = 0             ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par           = matrix(NA, nrow=N_mcmc, ncol=4)
colnames(MCMC_par) = c("b_nanp", "b_av", "loglike", "prior")

## Implement MCMC iterations
par_MC = c( 15000, 10 )                       ## (b_ab, bav)
Sigma_MC = diag( (0.65*par_MC)^2 )            ## Initial guess of covariance of MVN proposal dist
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
      cov_MC = cov( MCMC_par[1:(mc-1),1:2] )
    }
  }
  
  MCMC_par[mc,1:2] = par_MC
  MCMC_par[mc,3]   = loglike_MC
  MCMC_par[mc,4]   = prior_MC
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine posterior distributions  
## remove burn in 
MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]

par(mfrow=c(2,2))

hist(MCMC_burn[,1], breaks=50, main="b_nanp")
hist(MCMC_burn[,2], breaks=50, main="b_av")

plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,1], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn), y=MCMC_burn[,2], pch=19, cex=0.25, col="grey", main="b_av")

## Median and 95% Credible intervals 
quantile( MCMC_burn[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

#-------------------------------model_M3 : ab-Expon, av-Hill-----------------------------------------------------------------
# MODEL LIKELIHOOD
model_M3 <- function(par_M3){
  b_nanp     <- par_M3[1]
  b_av       <- par_M3[2]
  a_av       <- par_M3[3]
  
  ##############################
  ## Dose response
  
  DR <- exp( -log(2)*(nanp/b_nanp) )*( 1/(1+(avidity/b_av)^a_av) )
  DR[which(DR<min_DR)] <- min_DR
 
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      
      logL <- logL + sum(log(DR))
    }
    
    if( infect[j]==0 ){
      logL <- logL + sum(log(1-DR))
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
  b_nanp     <- par_M3[1]
  b_av       <- par_M3[2]
  a_av       <- par_M3[3]

  ######################################
  ## Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
   # prior_b_nanp = log(dgamma(b_nanp, shape=5, rate=0.0005))
    prior_b_nanp = log(500/150000)
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Uniform prior on b_av 
  
  if( b_av>1 && b_av<100 )
  {
   # prior_b_av = log(dgamma(b_av, shape=6, rate=0.5))
    prior_b_av = log(1/100)
  }else{
    prior_b_av = -LARGE
  }
  
  ######################################
  ##  Uniform prior on a_av ~ U(0,30)
  
  if( a_av>0 && a_av<30 )
  {
    prior_a_av = log(1/30)
  }else{
    prior_a_av = -LARGE
  }
  
  
  prior =  prior_b_nanp + prior_b_av + prior_a_av
  
  prior
}
prior_M3 = cmpfun(prior_M3, options=list(optimize=3))

##  MCMC 
N_mcmc       = 800000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 10000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 11000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1             ## Scaler for step size
MCMC_accept = 0             ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M3           = matrix(NA, nrow=N_mcmc, ncol=5)
colnames(MCMC_par_M3) = c("b_nanp", "b_av", "a_av", "loglike", "prior")

## Implement MCMC iterations
par_MC3 = c( 15000, 10 , 1 )            ## (b_nanp, b_av, a_av)
Sigma_MC3 = diag( (0.65*par_MC3)^2 )    ## Initial guess of covariance of MVN proposal dist
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
      cov_MC = cov( MCMC_par_M3[1:(mc-1),1:3] )
    }
  }
  
  
  MCMC_par_M3[mc,1:3] = par_MC3
  MCMC_par_M3[mc,4]   = loglike_MC3
  MCMC_par_M3[mc,5]   = prior_MC3
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and positerior distributions 
par(mfrow=c(2,3))
MCMC_burn_M3 <- MCMC_par_M3[floor(0.2*nrow(MCMC_par_M3)):(nrow(MCMC_par_M3)-1),]

hist(MCMC_burn_M3[,1], breaks=50, main="b_nanp")
hist(MCMC_burn_M3[,2], breaks=50, main="b_av")
hist(MCMC_burn_M3[,3], breaks=50, main="a_av")

plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,1], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,2], pch=19, cex=0.25, col="grey", main="b_av")
plot(x=1:nrow(MCMC_burn_M3), y=MCMC_burn_M3[,3], pch=19, cex=0.25, col="grey", main="a_av")

### median parameter values and 95% Credible Intervals #### 
quantile( MCMC_burn_M3[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M3[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )


#-------------------------------- model_M4 : ab-Hill, av-Hill--------------------------------------------------------------------
## LIKELIHOOD
model_M4 <- function(par_M4){
  b_nanp     <- par_M4[1]
  a_nanp     <- par_M4[2]
  b_av       <- par_M4[3]
  a_av       <- par_M4[4]
  
  
  ##############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )*( 1/(1+(avidity/b_av)^a_av) )
  DR[which(DR<min_DR)] <- min_DR
  
  logL <- 0
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      
      logL <- logL + sum(log(DR))
    }
    
    if( infect[j]==0 ){
      logL <- logL + sum(log(1-DR))
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
  b_nanp     <- par_M4[1]
  a_nanp     <- par_M4[2]
  b_av       <- par_M4[3]
  a_av       <- par_M4[4]
  
  ######################################
  ## Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
   # prior_b_nanp = log(500/150000)
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Uniform prior on a_nanp ~ U(0,30)
  
  if( a_nanp>0 && a_nanp<30 )
  {
    prior_a_nanp = log(1/30)
  }else{
    prior_a_nanp = -LARGE
  }
  
  ######################################
  ##  uniform prior on b_av 
  
  if( b_av>1 && b_av<100 )
  {
    # prior_b_av = log(dgamma(b_av, shape=6, rate=0.5))
    prior_b_av = log(1/100)
  }else{
    prior_b_av = -LARGE
  }
  
  ######################################
  ## Uniform prior on a_av ~ U(0,30)
  
  if( a_av>0 && a_av<30 )
  {
    prior_a_av = log(1/30)
  }else{
    prior_a_av = -LARGE
  }
  
  
  prior =  prior_b_nanp + prior_a_nanp + prior_b_av + prior_a_av
  
  prior
}
prior_M4 = cmpfun(prior_M4, options=list(optimize=3))

##  MCMC 
N_mcmc       = 800000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 10000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 11000       ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1             ## Scaler for step size
MCMC_accept = 0             ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M4           = matrix(NA, nrow=N_mcmc, ncol=6)
colnames(MCMC_par_M4) = c("b_nanp", "a_nanp", "b_av", "a_av", "loglike", "prior")

## Implement MCMC iterations
par_MC4 = c(4000, 1, 10 , 1 )       ## ( b_nanp, a_nanp, b_av, a_av)
Sigma_MC4 = diag( (0.85*par_MC4)^2 )  ## Initial guess of covariance of MVN proposal dist
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
      cov_MC = cov( MCMC_par_M4[1:(mc-1),1:4] )
    }
  }
  
  
  MCMC_par_M4[mc,1:4] = par_MC4
  MCMC_par_M4[mc,5]   = loglike_MC4
  MCMC_par_M4[mc,6]   = prior_MC4
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distributions 

MCMC_burn_M4 <- MCMC_par_M4[floor(0.2*nrow(MCMC_par_M4)):(nrow(MCMC_par_M4)-1),]

par(oma=c(1,1,1,1)) 
par(mar=c(4,1,1,1) + 0.1)
par(mfrow=c(2,4))
options(scipen=10)

hist(MCMC_burn_M4[,1], breaks=50, xlab=expression(beta["t"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M4[,2], breaks=50, xlab=expression(alpha["t"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M4[,3], breaks=50, xlab=expression(beta["ai"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)
hist(MCMC_burn_M4[,4], breaks=50, xlab=expression(alpha["ai"]), main = " ", ylab = " ", yaxt="n")
axis(2, labels = FALSE)

plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,1], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(beta["t"]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,2], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(alpha["t"]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,3], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(beta["ai"]), xlab="Iteration")
plot(x=1:nrow(MCMC_burn_M4), y=MCMC_burn_M4[,4], pch=19, cex=0.25, col="grey", main=" ", ylab = expression(alpha["ai"]), xlab="Iteration")

## median parameter values and 95% Credible Intervals
quantile( MCMC_burn_M4[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M4[,4], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

#------------------------------------------ model_M5 ab-Hill , av-Expon--------------------------------------------------
## LIKELIHOOD
model_M5 <- function(par_M5){
  b_nanp     <- par_M5[1]
  a_nanp     <- par_M5[2]
  b_av       <- par_M5[3]
  
  #############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )* exp( -log(2)*(avidity/b_av) )
  DR[which(DR<min_DR)] <- min_DR
  
  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      
      logL <- logL + sum(log(DR))
    }
    
    if( infect[j]==0 ){
      logL <- logL + sum(log(1-DR))
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
  b_nanp     <- par_M5[1]
  a_nanp     <- par_M5[2]
  b_av       <- par_M5[3]
  
  ######################################
  ## Gamma prior on b_nanp 
  
  if( b_nanp>500 && b_nanp<150000 )
  {
    prior_b_nanp = log(dgamma(b_nanp, shape=0.5, rate=0.00003))
    # prior_b_nanp = log(500/150000)
  }else{
    prior_b_nanp = -LARGE
  }
  
  ######################################
  ##  Uniform prior on a_nanp ~ U(0,20)
  
  if( a_nanp>0 && a_nanp<30 )
  {
    prior_a_nanp = log(1/30)
  }else{
    prior_a_nanp = -LARGE
  }
  
  ######################################
  ##  U prior on b_av 
  
  if( b_av>1 && b_av<100 )
  {
    #prior_b_av = log(dgamma(b_av, shape=6, rate=0.5))
    prior_b_av = log(1/100)
  }else{
    prior_b_av = -LARGE
  }
  
  
  prior =  prior_b_nanp + prior_a_nanp + prior_b_av 
  
  prior
}
prior_M5 = cmpfun(prior_M5, options=list(optimize=3))

##  MCMC 
N_mcmc       = 800000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 10000      ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 11000      ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M5           = matrix(NA, nrow=N_mcmc, ncol=5)
colnames(MCMC_par_M5) = c("b_nanp", "a_nanp", "b_av",  "loglike", "prior")

## Implement MCMC iterations
par_MC5     = c( 15000, 1, 10  )       ## (n, sig_n, sig_mu, b_nanp, a_nanp, b_av )
Sigma_MC5   = diag( (0.7*par_MC5)^2 )  ## Initial guess of covariance of MVN proposal dist
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
      cov_MC = cov( MCMC_par_M5[1:(mc-1),1:3] )
    }
  }
  
  MCMC_par_M5[mc,1:3] = par_MC5
  MCMC_par_M5[mc,4]   = loglike_MC5
  MCMC_par_M5[mc,5]   = prior_MC5
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distributions 
par(mfrow=c(2,3))

MCMC_burn_M5 <- MCMC_par_M5[floor(0.2*nrow(MCMC_par_M5)):(nrow(MCMC_par_M5)-1),]

hist(MCMC_burn_M5[,1], breaks=50, main="b_nanp")
hist(MCMC_burn_M5[,2], breaks=50, main="a_nanp")
hist(MCMC_burn_M5[,3], breaks=50, main="b_av")

plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,1], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,2], pch=19, cex=0.25, col="grey", main="a_nanp")
plot(x=1:nrow(MCMC_burn_M5), y=MCMC_burn_M5[,3], pch=19, cex=0.25, col="grey", main="b_av")

### Median parameter values and 95% Credible Interval 
quantile( MCMC_burn_M5[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M5[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

##------------------------------------------model_m6 : Interaction dose-response--------------------------------------------------
## LIKELIHOOD
model_M6 <- function(par_M6){
  b_nanp    <- par_M6[1]
  b_av      <- par_M6[2]
  gamma     <- par_M6[3]
  
  ##############################
  ## Dose response
  DR <- exp( -0.5*log(2)*( (nanp/b_nanp+avidity/b_av) +
                             sqrt( (nanp/b_nanp+avidity/b_av)^2 + 4*gamma*(nanp/b_nanp)*(avidity/b_av) ) ) )
  DR[which(DR<min_DR)] <- min_DR

  logL <- 0
  
  for(j in 1:NN){
    if( infect[j]==1 ){
      
      logL <- logL + sum(log(DR))
    }
    
    if( infect[j]==0 ){
      logL <- logL + sum(log(1-DR))
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
  b_nanp    <- par_M6[1]
  b_av      <- par_M6[2]
  gamma     <- par_M6[3]
  
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
  
  prior =  prior_b_nanp +  prior_b_av + prior_gamma
  
  prior
}
prior_M6 = cmpfun(prior_M6, options=list(optimize=3))

##  MCMC 
N_mcmc       = 800000     ## Number of MCMC iterations
N_tune_start = 500        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   = 10000      ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      = 11000      ## End of adaptive scaling of proposal size with rm_scale 

step_scale  = 1           ## Scaler for step size
MCMC_accept = 0           ## Track the MCMC acceptance rate

## Prepare object for MCMC fitting
MCMC_par_M6           = matrix(NA, nrow=N_mcmc, ncol=5)
colnames(MCMC_par_M6) = c("b_nanp", "b_av", "gamma", "loglike", "prior")

## Implement MCMC iterations
par_MC6 = c( 15000,  10 , 1  )          ## (n, sig_n, sig_mu, b_nanp, b_av, gamma)
Sigma_MC6 = diag( (0.8*par_MC6)^2 )     ## Initial guess of covariance of MVN proposal dist
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
      cov_MC = cov( MCMC_par_M6[1:(mc-1),1:3] )
    }
  }
  
  
  MCMC_par_M6[mc,1:3] = par_MC6
  MCMC_par_M6[mc,4]   = loglike_MC6
  MCMC_par_M6[mc,5]   = prior_MC6
}

MCMC_accept = MCMC_accept/N_mcmc

## Examine MCMC chains and posterior distrbutions  
MCMC_burn_M6 <- MCMC_par_M6[floor(0.2*nrow(MCMC_par_M6)):(nrow(MCMC_par_M6)-1),]

par(mfrow=c(2,3))

hist(MCMC_burn_M6[,1], breaks=50, main="b_nanp")
hist(MCMC_burn_M6[,2], breaks=50, main="b_av")
hist(MCMC_burn_M6[,3], breaks=50, main="gamma")

plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,1], pch=19, cex=0.25, col="grey", main="b_nanp")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,2], pch=19, cex=0.25, col="grey", main="b_av")
plot(x=1:nrow(MCMC_burn_M6), y=MCMC_burn_M6[,3], pch=19, cex=0.25, col="grey", main="gamma")

## Median paramter values and 95% Credible Intervals 
quantile( MCMC_burn_M6[,1], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,2], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )
quantile( MCMC_burn_M6[,3], probs=c(0.025, 0.5, 0.975), na.rm=TRUE )

#-----------------------------------------DIC----------------------------------------------------
## Mean of the posterior
theta_bar_M2 = apply( X=MCMC_burn[,1:2], FUN=mean, MARGIN=2)
theta_bar_M3 = apply( X=MCMC_burn_M3[,1:3], FUN=mean, MARGIN=2)
theta_bar_M4 = apply( X=MCMC_burn_M4[,1:4], FUN=mean, MARGIN=2)
theta_bar_M5 = apply( X=MCMC_burn_M5[,1:3], FUN=mean, MARGIN=2)
theta_bar_M6 = apply( X=MCMC_burn_M5[,1:3], FUN=mean, MARGIN=2)

## Deviance at mean of the posterior
D_theta_bar_M2 = -2*loglike_M2_cs( theta_bar_M2 )
D_theta_bar_M3 = -2*loglike_M3_cs( theta_bar_M3 )
D_theta_bar_M4 = -2*loglike_M4_cs( theta_bar_M4 )
D_theta_bar_M5 = -2*loglike_M5_cs( theta_bar_M5 )
D_theta_bar_M6 = -2*loglike_M5_cs( theta_bar_M6 )

## Mean deviance (averaged over posterior)
D_bar_M2 = -2*mean( MCMC_burn[,3] - MCMC_burn[,4] )
D_bar_M3 = -2*mean( MCMC_burn_M3[,4] - MCMC_burn_M3[,5] )
D_bar_M4 = -2*mean( MCMC_burn_M4[,5] - MCMC_burn_M4[,6] )
D_bar_M5 = -2*mean( MCMC_burn_M5[,4] - MCMC_burn_M5[,5] )
D_bar_M6 = -2*mean( MCMC_burn_M6[,4] - MCMC_burn_M6[,5] )

## Effective number of model parameters
pD_M2 = D_bar_M2 - D_theta_bar_M2
pD_M3 = D_bar_M3 - D_theta_bar_M3
pD_M4 = D_bar_M4 - D_theta_bar_M4
pD_M5 = D_bar_M5 - D_theta_bar_M5
pD_M6 = D_bar_M6 - D_theta_bar_M6

## Estimate of DIC
DIC_M2 = pD_M2 + D_bar_M2
DIC_M3 = pD_M3 + D_bar_M3
DIC_M4 = pD_M4 + D_bar_M4
DIC_M5 = pD_M5 + D_bar_M5
DIC_M6 = pD_M6 + D_bar_M6

min(DIC_M2, DIC_M3, DIC_M4, DIC_M5, DIC_M6)

### create object to store data ###
spz_fit_b <- matrix(NA, nrow=5, ncol=8)

colnames(spz_fit_b) <- c("b_nanp", "a_nanp",  "b_av", "a_av", "ML", "DIC", "dDIC", "gamma")  

rownames(spz_fit_b) <- c("exp exp", "exp hill", "hill hill", "hill exp" , "interaction")

spz_fit_b[1,1] <- paste0(round(median(MCMC_burn[,1]), digits=0)," (",round(quantile(MCMC_burn[,1], probs = c(0.025)), digits=0),"-", round(quantile(MCMC_burn[,1], probs = c(0.975)), digits=0),")")
spz_fit_b[1,3] <- paste0(round(median(MCMC_burn[,2]), digits=1)," (",round(quantile(MCMC_burn[,2], probs = c(0.025)), digits=1),"-", round(quantile(MCMC_burn[,2], probs = c(0.975)), digits=1),")")
spz_fit_b[1,5] <- -median(MCMC_burn[,3])
spz_fit_b[1,6] <- as.numeric(DIC_M2)

spz_fit_b[2,1] <- paste0(round(median(MCMC_burn_M3[,1]),digits=0)," (",  round(quantile(MCMC_burn_M3[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M3[,1],probs = c(0.975)),digits=0),")")
spz_fit_b[2,3] <- paste0(round(median(MCMC_burn_M3[,2]), digits=1)," (", round(quantile(MCMC_burn_M3[,2], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M3[,2], probs = c(0.975)),digits=1),")")
spz_fit_b[2,4] <- paste0(round(median(MCMC_burn_M3[,3]), digits=1)," (", round(quantile(MCMC_burn_M3[,3], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M3[,3], probs = c(0.975)),digits=1),")")
spz_fit_b[2,5] <- -median(MCMC_burn[,4])
spz_fit_b[2,6] <- as.numeric(DIC_M3)

spz_fit_b[3,1] <- paste0(round(median(MCMC_burn_M4[,1]),digits=0)," (",  round(quantile(MCMC_burn_M4[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M4[,1],probs = c(0.975)),digits=0),")")
spz_fit_b[3,2] <- paste0(round(median(MCMC_burn_M4[,2]),digits=1)," (",  round(quantile(MCMC_burn_M4[,2], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M4[,2],probs = c(0.975)),digits=1),")")
spz_fit_b[3,3] <- paste0(round(median(MCMC_burn_M4[,3]),digits=1)," (",  round(quantile(MCMC_burn_M4[,3], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M4[,3],probs = c(0.975)),digits=1),")")
spz_fit_b[3,4] <- paste0(round(median(MCMC_burn_M4[,4]),digits=1)," (",  round(quantile(MCMC_burn_M4[,4], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M4[,4],probs = c(0.975)),digits=1),")")
spz_fit_b[3,5] <- -median(MCMC_burn_M4[,5])
spz_fit_b[3,6] <- as.numeric(DIC_M4)

spz_fit_b[4,1] <- paste0(round(median(MCMC_burn_M5[,1]),digits=0)," (",  round(quantile(MCMC_burn_M5[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M5[,1],probs = c(0.975)),digits=0),")")
spz_fit_b[4,2] <- paste0(round(median(MCMC_burn_M5[,2]),digits=1)," (",  round(quantile(MCMC_burn_M5[,2], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M5[,2],probs = c(0.975)),digits=1),")")
spz_fit_b[4,3] <- paste0(round(median(MCMC_burn_M5[,3]),digits=1)," (",  round(quantile(MCMC_burn_M5[,3], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M5[,3],probs = c(0.975)),digits=1),")")
spz_fit_b[4,5] <- -median(MCMC_burn_M5[,4])
spz_fit_b[4,6] <- as.numeric(DIC_M5)

spz_fit_b[5,1] <- paste0(round(median(MCMC_burn_M6[,1]),digits=0)," (",  round(quantile(MCMC_burn_M6[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M6[,1],probs = c(0.975)),digits=0),")")
spz_fit_b[5,3] <- paste0(round(median(MCMC_burn_M6[,2]),digits=1)," (",  round(quantile(MCMC_burn_M6[,2], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M6[,2],probs = c(0.975)),digits=1),")")
spz_fit_b[5,8] <- paste0(round(median(MCMC_burn_M6[,3]),digits=1)," (",  round(quantile(MCMC_burn_M6[,3], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M6[,3],probs = c(0.975)),digits=1),")")
spz_fit_b[5,5] <- -median(MCMC_burn_M6[,4])
spz_fit_b[5,6] <- as.numeric(DIC_M6)

write.csv(spz_fit_b, "binary_fitting.csv")

#-------------------------------------estimated vaccine efficacy based on best fitting model-------------------------------------------------------------------
par_median = apply(X=MCMC_burn_M4[,1:4], MARGIN=2, FUN=median)
b_nanp <- as.numeric(par_median[1])
a_nanp <- as.numeric(par_median[2])
b_av   <- as.numeric(par_median[3])
a_av   <- as.numeric(par_median[4])

## VE calculation 
model_VE <- function(nanp, avidity, par_M4){
  b_nanp     <- par_M4[1]
  a_nanp     <- par_M4[2]
  b_av       <- par_M4[3]
  a_av       <- par_M4[4]
  
  
  ##############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )*( 1/(1+(avidity/b_av)^a_av) )
  DR[which(DR<min_DR)] <- min_DR
  
  VE_inf <-  1 - DR 
  VE_inf <- mean(VE_inf)
  VE_inf
  
}

N_sam = 20000

sam_seq = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M2_sam_test = matrix(NA, nrow=N_sam, ncol=length(avidity))

for(k in 1:N_sam){
  M2_sam_test[k,] = sapply(avidity, model_VE, nanp=nanp, par=MCMC_burn_M4[sam_seq[k],1:4])
}

M2_quant_test = matrix(NA, nrow=3, ncol=length(avidity))

for(j in 1:length(avidity)){
  M2_quant_test[,j] = quantile( M2_sam_test[,j], prob=c(0.025, 0.5, 0.975) )
}

## CIs and Point Estimate of VE across the study (all Volunteers)
# lower CI
mean(M2_quant_test[1,])
#value
mean(M2_quant_test[2,])
#upper CI
mean(M2_quant_test[3,])
