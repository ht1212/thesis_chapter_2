#########################################################################################################################
## Code to recreate the results from the section: Vaccine Efficacy Per Sporozoite  
## Per parasite efficacy is estimating the proportional reduction in the number of parasites emerging from the liver
#########################################################################################################################

### per parasite efficacy model 

model_VE_par <- function(nanp, avidity, par_M5){
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
  
  #############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )* exp( -log(2)*(avidity/b_av) )
  
  VE_par <- 1 - DR 
  
  VE_par <- mean(VE_par)
  
}

model_VE_par = cmpfun(model_VE_par, options=list(optimize=3))

N_sam = 20000

sam_seq = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M2_sam_par = matrix(NA, nrow=N_sam, ncol=length(avidity))

for(k in 1:N_sam){
  M2_sam_par[k,] = sapply(avidity, model_VE_par, nanp=nanp, par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_par = matrix(NA, nrow=3, ncol=length(avidity))

for(j in 1:length(avidity)){
  M2_quant_par[,j] = quantile( M2_sam_par[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_par[1,])
#value
mean(M2_quant_par[2,])
#upper CI
mean(M2_quant_par[3,])

