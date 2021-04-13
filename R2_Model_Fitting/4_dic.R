##################################################################################################################################
## Code contained in this file calculate the Deviance Inofrmation Criterion for all of the fitted sporozoite infection models.   #
## DIC is the Deviance Information Criterion used for model selection. Accounting for parameter number                           #
## Note that the deviance is -2*log(likelihood)                                                                                  # 
##################################################################################################################################

######################################
## Mean of the posterior
theta_bar_M2 = apply( X=MCMC_burn[,1:5], FUN=mean, MARGIN=2)
theta_bar_M3 = apply( X=MCMC_burn_M3[,1:6], FUN=mean, MARGIN=2)
theta_bar_M4 = apply( X=MCMC_burn_M4[,1:7], FUN=mean, MARGIN=2)
theta_bar_M5 = apply( X=MCMC_burn_M5[,1:6], FUN=mean, MARGIN=2)
theta_bar_M6 = apply( X=MCMC_burn_M6[,1:6], FUN=mean, MARGIN=2)
theta_bar_M7 = apply( X=MCMC_burn_M7[,1:4], FUN=mean, MARGIN=2)
theta_bar_M8 = apply( X=MCMC_burn_M8[,1:5], FUN=mean, MARGIN=2)
theta_bar_M9 = apply( X=MCMC_burn_M9[,1:4], FUN=mean, MARGIN=2)
theta_bar_M10 = apply( X=MCMC_burn_M10[,1:5], FUN=mean, MARGIN=2)
theta_bar_M11 = apply( X=MCMC_burn_M11[,1:4], FUN=mean, MARGIN=2)

######################################
## Deviance at mean of the posterior
D_theta_bar_M2 = -2*loglike_M2_cs( theta_bar_M2 )
D_theta_bar_M3 = -2*loglike_M3_cs( theta_bar_M3 )
D_theta_bar_M4 = -2*loglike_M4_cs( theta_bar_M4 )
D_theta_bar_M5 = -2*loglike_M5_cs( theta_bar_M5 )
D_theta_bar_M6 = -2*loglike_M6_cs( theta_bar_M6 )
D_theta_bar_M7 = -2*loglike_M7_cs( theta_bar_M7 )
D_theta_bar_M8 = -2*loglike_M8_cs( theta_bar_M8 )
D_theta_bar_M9 = -2*loglike_M9_cs( theta_bar_M9 )
D_theta_bar_M10 = -2*loglike_M10_cs( theta_bar_M10 )
D_theta_bar_M11 = -2*loglike_M11_cs( theta_bar_M11 )

####################################
## Mean deviance (averaged over posterior)
D_bar_M2 = -2*mean( MCMC_burn[,6] - MCMC_burn[,7] )
D_bar_M3 = -2*mean( MCMC_burn_M3[,7] - MCMC_burn_M3[,8] )
D_bar_M4 = -2*mean( MCMC_burn_M4[,8] - MCMC_burn_M4[,9] )
D_bar_M5 = -2*mean( MCMC_burn_M5[,7] - MCMC_burn_M5[,8] )
D_bar_M6 = -2*mean( MCMC_burn_M6[,7] - MCMC_burn_M6[,8] )
D_bar_M7 = -2*mean( MCMC_burn_M7[,5] - MCMC_burn_M7[,6] )
D_bar_M8 = -2*mean( MCMC_burn_M8[,6] - MCMC_burn_M8[,7] )
D_bar_M9 = -2*mean( MCMC_burn_M9[,5] - MCMC_burn_M9[,6] )
D_bar_M10 = -2*mean( MCMC_burn_M10[,6] - MCMC_burn_M10[,7] )
D_bar_M11 = -2*mean( MCMC_burn_M11[,5] - MCMC_burn_M11[,6] )

######################################
## Effective number of model parameters
pD_M2 = D_bar_M2 - D_theta_bar_M2
pD_M3 = D_bar_M3 - D_theta_bar_M3
pD_M4 = D_bar_M4 - D_theta_bar_M4
pD_M5 = D_bar_M5 - D_theta_bar_M5
pD_M6 = D_bar_M6 - D_theta_bar_M6
pD_M7 = D_bar_M7 - D_theta_bar_M7
pD_M8 = D_bar_M8 - D_theta_bar_M8
pD_M9 = D_bar_M9 - D_theta_bar_M9
pD_M10 = D_bar_M10 - D_theta_bar_M10
pD_M11 = D_bar_M11 - D_theta_bar_M11

######################################
## Estimate of DIC
DIC_M2 = pD_M2 + D_bar_M2
DIC_M3 = pD_M3 + D_bar_M3
DIC_M4 = pD_M4 + D_bar_M4
DIC_M5 = pD_M5 + D_bar_M5
DIC_M6 = pD_M6 + D_bar_M6
DIC_M7 = pD_M7 + D_bar_M7
DIC_M8 = pD_M8 + D_bar_M8
DIC_M9 = pD_M9 + D_bar_M9
DIC_M10 = pD_M10 + D_bar_M10
DIC_M11 = pD_M11 + D_bar_M11