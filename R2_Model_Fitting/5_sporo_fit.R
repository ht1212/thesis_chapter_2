### create object to store data ###
spz_fit <- matrix(NA, nrow=10, ncol=12)

colnames(spz_fit) <- c("n", "sig_n", "sig_mu", "b_nanp", "a_nanp", 
                       "b_av", "a_av", "gamma", "ML", "DIC", "dDIC" , "VE")  

rownames(spz_fit) <- c("exp exp", "exp hill", "hill hill", "hill exp", "interaction" , "ab exp" 
                       , "ab hil" , "av exp" , "av hill" , "base" )

spz_fit[1,1] <- paste0(round(median(MCMC_burn[,1]), digits=0)," (",round(quantile(MCMC_burn[,1], probs = c(0.025)), digits=0),"-", round(quantile(MCMC_burn[,1], probs = c(0.975)), digits=0),")")
spz_fit[1,2] <- paste0(round(median(MCMC_burn[,2]), digits=0)," (",round(quantile(MCMC_burn[,2], probs = c(0.025)), digits=0),"-", round(quantile(MCMC_burn[,2], probs = c(0.975)), digits=0),")")
spz_fit[1,3] <- paste0(round(median(MCMC_burn[,3]), digits=0)," (",round(quantile(MCMC_burn[,3], probs = c(0.025)), digits=0),"-", round(quantile(MCMC_burn[,3], probs = c(0.975)), digits=0),")")
spz_fit[1,4] <- paste0(round(median(MCMC_burn[,4]), digits=0)," (",round(quantile(MCMC_burn[,4], probs = c(0.025)), digits=0),"-", round(quantile(MCMC_burn[,4], probs = c(0.975)), digits=0),")")
spz_fit[1,6] <- paste0(round(median(MCMC_burn[,5]), digits=1)," (",round(quantile(MCMC_burn[,5], probs = c(0.025)), digits=1),"-", round(quantile(MCMC_burn[,5], probs = c(0.975)), digits=1),")")
spz_fit[1,9]           <- -median(MCMC_burn[,6])
spz_fit[1,10]           <- as.numeric(DIC_M2)


spz_fit[2,1] <- paste0(round(median(MCMC_burn_M3[,1]),digits=0)," (",  round(quantile(MCMC_burn_M3[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M3[,1],probs = c(0.975)),digits=0),")")
spz_fit[2,2] <- paste0(round(median(MCMC_burn_M3[,2]), digits=0)," (", round(quantile(MCMC_burn_M3[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M3[,2], probs = c(0.975)),digits=0),")")
spz_fit[2,3] <- paste0(round(median(MCMC_burn_M3[,3]), digits=0)," (", round(quantile(MCMC_burn_M3[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M3[,3], probs = c(0.975)),digits=0),")")
spz_fit[2,4] <- paste0(round(median(MCMC_burn_M3[,4]), digits=0)," (", round(quantile(MCMC_burn_M3[,4], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M3[,4], probs = c(0.975)),digits=0),")")
spz_fit[2,6] <- paste0(round(median(MCMC_burn_M3[,5]), digits=1)," (", round(quantile(MCMC_burn_M3[,5], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M3[,5], probs = c(0.975)),digits=1),")")
spz_fit[2,7] <- paste0(round(median(MCMC_burn_M3[,6]), digits=1)," (", round(quantile(MCMC_burn_M3[,6], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M3[,6], probs = c(0.975)),digits=1),")")
spz_fit[2,9]           <- -median(MCMC_burn[,7])
spz_fit[2,10]           <- as.numeric(DIC_M3)


spz_fit[3,1] <- paste0(round(median(MCMC_burn_M4[,1]),digits=0)," (",  round(quantile(MCMC_burn_M4[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M4[,1],probs = c(0.975)),digits=0),")")
spz_fit[3,2] <- paste0(round(median(MCMC_burn_M4[,2]),digits=0)," (",  round(quantile(MCMC_burn_M4[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M4[,2],probs = c(0.975)),digits=0),")")
spz_fit[3,3] <- paste0(round(median(MCMC_burn_M4[,3]),digits=0)," (",  round(quantile(MCMC_burn_M4[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M4[,3],probs = c(0.975)),digits=0),")")
spz_fit[3,4] <- paste0(round(median(MCMC_burn_M4[,4]),digits=0)," (",  round(quantile(MCMC_burn_M4[,4], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M4[,4],probs = c(0.975)),digits=0),")")
spz_fit[3,5] <- paste0(round(median(MCMC_burn_M4[,5]),digits=1)," (",  round(quantile(MCMC_burn_M4[,5], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M4[,5],probs = c(0.975)),digits=1),")")
spz_fit[3,6] <- paste0(round(median(MCMC_burn_M4[,6]),digits=1)," (",  round(quantile(MCMC_burn_M4[,6], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M4[,6],probs = c(0.975)),digits=1),")")
spz_fit[3,7] <- paste0(round(median(MCMC_burn_M4[,7]),digits=1)," (",  round(quantile(MCMC_burn_M4[,7], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M4[,7],probs = c(0.975)),digits=1),")")
spz_fit[3,9]             <-   -median(MCMC_burn_M4[,8])
spz_fit[3,10]           <-     as.numeric(DIC_M4)


spz_fit[4,1] <- paste0(round(median(MCMC_burn_M5[,1]),digits=0)," (",  round(quantile(MCMC_burn_M5[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M5[,1],probs = c(0.975)),digits=0),")")
spz_fit[4,2] <- paste0(round(median(MCMC_burn_M5[,2]),digits=0)," (",  round(quantile(MCMC_burn_M5[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M5[,2],probs = c(0.975)),digits=0),")")
spz_fit[4,3] <- paste0(round(median(MCMC_burn_M5[,3]),digits=0)," (",  round(quantile(MCMC_burn_M5[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M5[,3],probs = c(0.975)),digits=0),")")
spz_fit[4,4] <- paste0(round(median(MCMC_burn_M5[,4]),digits=0)," (",  round(quantile(MCMC_burn_M5[,4], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M5[,4],probs = c(0.975)),digits=0),")")
spz_fit[4,5] <- paste0(round(median(MCMC_burn_M5[,5]),digits=1)," (",  round(quantile(MCMC_burn_M5[,5], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M5[,5],probs = c(0.975)),digits=1),")")
spz_fit[4,6] <- paste0(round(median(MCMC_burn_M5[,6]),digits=1)," (",  round(quantile(MCMC_burn_M5[,6], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M5[,6],probs = c(0.975)),digits=1),")")
spz_fit[4,9]             <-   -median(MCMC_burn_M5[,7])
spz_fit[4,10]           <-     as.numeric(DIC_M5)

spz_fit[5,1] <- paste0(round(median(MCMC_burn_M6[,1]),digits=0)," (",  round(quantile(MCMC_burn_M6[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M6[,1],probs = c(0.975)),digits=0),")")
spz_fit[5,2] <- paste0(round(median(MCMC_burn_M6[,2]),digits=0)," (",  round(quantile(MCMC_burn_M6[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M6[,2],probs = c(0.975)),digits=0),")")
spz_fit[5,3] <- paste0(round(median(MCMC_burn_M6[,3]),digits=0)," (",  round(quantile(MCMC_burn_M6[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M6[,3],probs = c(0.975)),digits=0),")")
spz_fit[5,4] <- paste0(round(median(MCMC_burn_M6[,4]),digits=0)," (",  round(quantile(MCMC_burn_M6[,4], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M6[,4],probs = c(0.975)),digits=0),")")
spz_fit[5,6] <- paste0(round(median(MCMC_burn_M6[,5]),digits=1)," (",  round(quantile(MCMC_burn_M6[,5], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M6[,5],probs = c(0.975)),digits=1),")")
spz_fit[5,8] <- paste0(round(median(MCMC_burn_M6[,6]),digits=1)," (",  round(quantile(MCMC_burn_M6[,6], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M6[,6],probs = c(0.975)),digits=1),")")
spz_fit[5,9]             <-   -median(MCMC_burn_M6[,7])
spz_fit[5,10]           <-     as.numeric(DIC_M6)

spz_fit[6,1] <- paste0(round(median(MCMC_burn_M7[,1]),digits=0)," (",  round(quantile(MCMC_burn_M7[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M7[,1],probs = c(0.975)),digits=0),")")
spz_fit[6,2] <- paste0(round(median(MCMC_burn_M7[,2]),digits=0)," (",  round(quantile(MCMC_burn_M7[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M7[,2],probs = c(0.975)),digits=0),")")
spz_fit[6,3] <- paste0(round(median(MCMC_burn_M7[,3]),digits=0)," (",  round(quantile(MCMC_burn_M7[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M7[,3],probs = c(0.975)),digits=0),")")
spz_fit[6,4] <- paste0(round(median(MCMC_burn_M7[,4]),digits=0)," (",  round(quantile(MCMC_burn_M7[,4], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M7[,4],probs = c(0.975)),digits=0),")")
spz_fit[6,9]             <-   -median(MCMC_burn_M7[,5])
spz_fit[6,10]           <-     as.numeric(DIC_M7)


spz_fit[7,1] <- paste0(round(median(MCMC_burn_M8[,1]),digits=0)," (",  round(quantile(MCMC_burn_M8[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M8[,1],probs = c(0.975)),digits=0),")")
spz_fit[7,2] <- paste0(round(median(MCMC_burn_M8[,2]),digits=0)," (",  round(quantile(MCMC_burn_M8[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M8[,2],probs = c(0.975)),digits=0),")")
spz_fit[7,3] <- paste0(round(median(MCMC_burn_M8[,3]),digits=0)," (",  round(quantile(MCMC_burn_M8[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M8[,3],probs = c(0.975)),digits=0),")")
spz_fit[7,4] <- paste0(round(median(MCMC_burn_M8[,4]),digits=0)," (",  round(quantile(MCMC_burn_M8[,4], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M8[,4],probs = c(0.975)),digits=0),")")
spz_fit[7,5] <- paste0(round(median(MCMC_burn_M8[,5]),digits=1)," (",  round(quantile(MCMC_burn_M8[,5], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M8[,5],probs = c(0.975)),digits=1),")")
spz_fit[7,9]             <-   -median(MCMC_burn_M8[,6])
spz_fit[7,10]           <-     as.numeric(DIC_M8)

spz_fit[8,1] <- paste0(round(median(MCMC_burn_M9[,1]),digits=0)," (",  round(quantile(MCMC_burn_M9[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M9[,1],probs = c(0.975)),digits=0),")")
spz_fit[8,2] <- paste0(round(median(MCMC_burn_M9[,2]),digits=0)," (",  round(quantile(MCMC_burn_M9[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M9[,2],probs = c(0.975)),digits=0),")")
spz_fit[8,3] <- paste0(round(median(MCMC_burn_M9[,3]),digits=0)," (",  round(quantile(MCMC_burn_M9[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M9[,3],probs = c(0.975)),digits=0),")")
spz_fit[8,6] <- paste0(round(median(MCMC_burn_M9[,4]),digits=1)," (",  round(quantile(MCMC_burn_M9[,4], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M9[,4],probs = c(0.975)),digits=1),")")
spz_fit[8,9]             <-   -median(MCMC_burn_M9[,5])
spz_fit[8,10]           <-     as.numeric(DIC_M9)

spz_fit[9,1] <- paste0(round(median(MCMC_burn_M10[,1]),digits=0)," (",  round(quantile(MCMC_burn_M10[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M10[,1],probs = c(0.975)),digits=0),")")
spz_fit[9,2] <- paste0(round(median(MCMC_burn_M10[,2]),digits=0)," (",  round(quantile(MCMC_burn_M10[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M10[,2],probs = c(0.975)),digits=0),")")
spz_fit[9,3] <- paste0(round(median(MCMC_burn_M10[,3]),digits=0)," (",  round(quantile(MCMC_burn_M10[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M10[,3],probs = c(0.975)),digits=0),")")
spz_fit[9,6] <- paste0(round(median(MCMC_burn_M10[,4]),digits=1)," (",  round(quantile(MCMC_burn_M10[,4], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M10[,4],probs = c(0.975)),digits=1),")")
spz_fit[9,7] <- paste0(round(median(MCMC_burn_M10[,5]),digits=1)," (",  round(quantile(MCMC_burn_M10[,5], probs = c(0.025)),digits=1),"-",round(quantile(MCMC_burn_M10[,5],probs = c(0.975)),digits=1),")")
spz_fit[9,9]             <-   -median(MCMC_burn_M10[,6])
spz_fit[9,10]           <-     as.numeric(DIC_M10)


spz_fit[10,1]  <- paste0(round(median(MCMC_burn_M11[,1]),digits=0)," (",  round(quantile(MCMC_burn_M11[,1], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M11[,1],probs = c(0.975)),digits=0),")")
spz_fit[10,2]  <- paste0(round(median(MCMC_burn_M11[,2]),digits=0)," (",  round(quantile(MCMC_burn_M11[,2], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M11[,2],probs = c(0.975)),digits=0),")")
spz_fit[10,3]  <- paste0(round(median(MCMC_burn_M11[,3]),digits=0)," (",  round(quantile(MCMC_burn_M11[,3], probs = c(0.025)),digits=0),"-",round(quantile(MCMC_burn_M11[,3],probs = c(0.975)),digits=0),")")
spz_fit[10,12] <- paste0(round(median(MCMC_burn_M11[,4]),digits=2)," (",  round(quantile(MCMC_burn_M11[,4], probs = c(0.025)),digits=2),"-",round(quantile(MCMC_burn_M11[,4],probs = c(0.975)),digits=2),")")
spz_fit[10,9]             <-   -median(MCMC_burn_M11[,5])
spz_fit[10,10]           <-     as.numeric(DIC_M11)

write.csv(spz_fit, "model_fitting.csv")
### Difference in DIC ###
spz_fit[,11] <- spz_fit[,10] - min(spz_fit[,10])

### order table based on difference in DIC ### 
spz_fit <- spz_fit[order(spz_fit[,11]),]
