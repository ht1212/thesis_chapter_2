##########################################################################################################################################
## Modelling the Roles of Antibody Titre and Avidity in Protection from Malaria Infection Following RTS,S/AS01 Vaccination
##
## This paper represents further use of the White et al Sporozoite Infection Model first published in 2013 
## (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061395). In it we use the model to 
## characterise assocations between antibody number and antibody quality in RTS,S/AS01 vaccine induced protection
## against malaria infection via mosquito challenge. The model code used here was kindly developed from Michael White's
## previous work on this. 
##
## Any questions, errors I've missed or general comments please get in touch: h.thompson12@imperial.ac.uk
##
## MCMC Model Comparisons for results section Model Fitting. 
##
## Visualising Sporozoite Infection Model Comparisons  
## The M5 model parameters are used here to calculate vaccine efficacy as a funciton of individual's immune responses 
## and predicted time to onset of infection given immune measurements. Comparisons are made between models which took into account 
## both antibody titre and antibody avidity and that which only used titre to predict efficacy. 
###########################################################################################################################################

################
## Parameters ##
################

par_median = apply(X=MCMC_burn_M5[,1:7], MARGIN=2, FUN=median)

n <- as.numeric(par_median[1])
sig_n <- as.numeric(par_median[2])
sig_mu <- as.numeric(par_median[3])
b_nanp <- as.numeric(par_median[4])
a_nanp <- as.numeric(par_median[5])
b_av <- as.numeric(par_median[6])
LogL <- as.numeric(par_median[7])

p <- (sig_n^2-n)/(sig_n^2)
r <- (n^2)/(sig_n^2-n)

#################################
## Time To Onset of Infection  ##
#################################

#======================================================================================
## 1. Trial Data  

## creating survival objects  
status <- infect  
time <- time  
trt <- rep(NA,46)
trt[vaccine=="fx"] <- 1
trt[vaccine=="std"]<- 0 

df_observed <- data.frame(status, time, trt)

#=======================================================================================
## 2. Antibody Titre and Avidity Model Predictions of Time to Onset

## model predicted survival object  
mult <- 3.8
pfT <- 50000000
t_L <- 6.5

N_imm <- 2000
N_sam <- 500000

nanp_infect <- nanp[which(infect==1)]
avidity_infect <- avidity[which(infect==1)]
Time_inf <- time[which(infect==1)]

Time_infect <- rep(NA[1:10])

p_nb <- (sig_n^2-n)/(sig_n^2)
r_nb <- (n^2)/(sig_n^2-n)


for (i in 1:length(Time_infect)) {
  n_DR <- n * (( 1/(1+(nanp_infect[i]/b_nanp)^a_nanp)) * exp(-log(2)*avidity_infect[i]/b_av))
  n_DR <- max(n_DR, 1e-5)
  p_nb_DR = n_DR/( n_DR + r_nb )
  
  SPZ <- rnbinom(N_sam, size=r_nb, prob=1-p_nb_DR)
  if( 0 %in% SPZ ){
    SPZ <- SPZ[-which(SPZ==0)]
  }
  
  theta = (sig_mu^2)/mu
  k     = (mu/sig_mu)^2
  
  mero <- rgamma(length(SPZ), scale=theta, shape=SPZ*k )
  
  TT <- t_L - log( mero/pfT )/log(mult)
  
  Time_infect[i] <- quantile(TT, prob=c(0.5))
  
}

## creating new objects for survival plots 
Time_infect
Time_inf
time
time_estimate <- time 
time_estimate[which(infect==1)] <- Time_infect
time_estimate

df_estimate <- data.frame(status, time_estimate, trt)
colnames(df_estimate) <- c("status", "time", "trt")

#========================================================================================
## 3. Compute survival curves by combining preditcted and observed survival functions 
obs <- survfit(Surv(time, status) ~ trt, data=df_observed)
est <- survfit(Surv(time, status) ~ trt, data=df_estimate)

fit <- list(OBS = obs, EST = est)

p1 <- ggsurvplot(fit,
                 combine = TRUE,   
                 break.time.by = 7 ,
                 linetype = c(1,1, 3, 3), 
                 legend.labs=c("Observed", "  ", "Model predicted", " "),
                 censor = FALSE,                         # Remove censor points
                 palette = c( "#FBBB48","#C23A4B", "#fabd50","#bd3e4e" ),#lacroix_palette("PommeBaya"),
                 fun = "event" ,
                 xlab="Time from challenge (days)", 
                 ylab = "Proportion infected",
                 xlim = c(0,28), 
                 #  linetype = c(1,1),
                 legend.title=" ", 
                 legend="bottom", 
                 ggtheme = theme_bw(11)) 

p1

p1 <- p1$plot
#========================================================================================
## 4. Antibody Titre Only Model Predicted Time To Onset  

### model parameters  
ab_par_median = apply(X=MCMC_burn_M8[,1:5], MARGIN=2, FUN=median)

n_ab <- as.numeric(ab_par_median[1])
sig_n_ab <- as.numeric(ab_par_median[2])
sig_mu_ab <- as.numeric(ab_par_median[3])
b_nanp_ab <- as.numeric(ab_par_median[4])
a_nanp_ab <- as.numeric(ab_par_median[5])

p_ab <- (sig_n_ab^2-n_ab)/(sig_n_ab^2)
r_ab <- (n_ab^2)/(sig_n_ab^2-n_ab)

### model predicted survival object 
Time_infect_ab <- rep(NA[1:10])

for (i in 1:length(Time_infect_ab)) {
  n_DR_ab <- n_ab * (( 1/(1+(nanp_infect[i]/b_nanp_ab)^a_nanp_ab))) 
  n_DR_ab <- max(n_DR_ab, 1e-5)
  p_nb_DR_ab = n_DR_ab/( n_DR_ab + r_ab )
  
  SPZ_ab <- rnbinom(N_sam, size=r_ab, prob=1-p_nb_DR_ab)
  if( 0 %in% SPZ_ab ){
    SPZ_ab <- SPZ_ab[-which(SPZ_ab==0)]
  }
  
  theta_ab = (sig_mu_ab^2)/mu
  k_ab     = (mu/sig_mu_ab)^2
  
  mero_ab <- rgamma(length(SPZ_ab), scale=theta_ab, shape=SPZ_ab*k_ab )
  
  TT_ab <- t_L - log( mero_ab/pfT )/log(mult)
  
  Time_infect_ab[i] <- quantile(TT_ab, prob=c(0.5))
  
}

Time_infect_ab
Time_infect
Time_inf

time
time_estimate_ab <- time 
time_estimate_ab[which(infect==1)] <- Time_infect_ab
time_estimate_ab

df_estimate_ab <- data.frame(status, time_estimate_ab, trt)
colnames(df_estimate_ab) <- c("status", "time", "trt")

#### compute survival curves by combining preditcted and observed survival functions 

obs    <- survfit(Surv(time, status) ~ trt, data=df_observed)
est_ab <- survfit(Surv(time, status) ~ trt, data=df_estimate_ab)

fit_ab <- list(OBS = obs,  AB = est_ab)

p2 <- ggsurvplot(fit_ab,
                 combine = TRUE,   
                 break.time.by = 7 ,
                 linetype = c(1,1, 3, 3), 
                 legend.labs=c("Observed", "  ", "Model predicted", " "),
                 censor = FALSE,                         # Remove censor points
                 palette = c( "#FBBB48","#C23A4B", "#fabd50","#bd3e4e" ),#lacroix_palette("PommeBaya"),
                 fun = "event" ,
                 xlab="Time from challenge (days)", 
                 ylab = "Proportion infected",
                 xlim = c(0,28), 
               #  linetype = c(1,1),
                 legend.title=" ", 
                 legend="bottom", 
                 ggtheme = theme_bw(11)) 

p2 

p2 <- p2$plot

###########################
## Mulitplot Comparison  ##
###########################

# ggsurvlist <- list(
#   x = p1,
#   y = p2
# )
# 
# # Arrange multiple ggsurvplots and print the output
# arrange_ggsurvplots(ggsurvlist, print = TRUE, ncol = 1, nrow = 2)

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


#=====================================================================================================================
## Predicted Vaccine Efficacy Model Comparisons 

## Values needed 
observed_VE_fx <-  0.867 
observed_VE_std <- 0.625 

## swapping event codes for Hosmer-Lemeshow
protected <- infect
protected[which(infect=="1")]<- 0
protected[which(infect==0)] <-1

#=====================================================================================================================
## Anitbody Avidity and Titre Model 
## Individual Vounteer Efficacy Predictions 
DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )* exp( -log(2)*(avidity/b_av) )
VE_inf <- (r/(n*DR+r))^r
VE_inf <- 1 - (1-VE_inf)/(1-(1-p)^r)

## Hosmer-Lemeshow to bin into 5 groups based on predicted eff 
hl.results_M5 <- hoslem.test(protected, VE_inf, g=5 ) 
hl.results_M5

# model for predicted efficacy 
model_VE <- function(nanp, avidity, par_M5){
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
  
  VE_inf <- (r/(n*DR+r))^r
  VE_inf <-  1 - (1-VE_inf)/(1-(1-p)^r)
  VE_inf <- mean(VE_inf)
  
  VE_inf
  
  
}

model_VE = cmpfun(model_VE, options=list(optimize=3))


## Overall efficacy predicted for fractional and stanard arms based on immune responses 
## Fractional 
N_sam = 1000

sam_seq_ve = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M5_sam_test_fx = matrix(NA, nrow=N_sam, ncol=length(avidity_fx))

for(k in 1:N_sam){
  M5_sam_test_fx[k,] = sapply(avidity_fx, model_VE, nanp=nanp_fx, par=MCMC_burn_M5[sam_seq_ve[k],1:6])
}

M5_quant_fx = matrix(NA, nrow=3, ncol=length(avidity_fx))

for(j in 1:length(avidity_fx)){
  M5_quant_fx[,j] = quantile( M5_sam_test_fx[,j], prob=c(0.025, 0.5, 0.975) )
}

## Standard
M5_sam_test_std = matrix(NA, nrow=N_sam, ncol=length(avidity_std))

for(k in 1:N_sam){
  M5_sam_test_std[k,] = sapply(avidity_std, model_VE, nanp=nanp_std, par=MCMC_burn_M5[sam_seq_ve[k],1:6])
}

M5_quant_std = matrix(NA, nrow=3, ncol=length(avidity_std))

for(j in 1:length(avidity_std)){
  M5_quant_std[,j] = quantile( M5_sam_test_std[,j], prob=c(0.025, 0.5, 0.975) )
}

### Plotting  
plot.data_M5 <- data.frame(predicted_ve=prop.table(hl.results_M5$expected,1)[,2],
                           
                           observed_ve=prop.table(hl.results_M5$observed, 1)[,2] )
confint.df <- binom.confint(
  hl.results_M5$observed[,2],
  margin.table(hl.results_M5$observed, 1),
  methods="wilson")

plot.data_M5$observed.rate.upper <- confint.df$upper.Freq
plot.data_M5$observed.rate.lower <- confint.df $lower.Freq

data.1_M5 <- data.frame(observed_VE_std)
data.1_M5$estimate <- median(M5_quant_std[2,])
data.1_M5$lower <- median(M5_quant_std[1,])
data.1_M5$higher <- median(M5_quant_std[3,])

data.2_M5 <- data.frame(observed_VE_fx)
data.2_M5$estimate <- median(M5_quant_fx[2,])
data.2_M5$lower <- median(M5_quant_fx[1,])
data.2_M5$higher <- median(M5_quant_fx[3,])

ab_av_plot <- ggplot() + 
  geom_pointrange(plot.data_M5,
                  mapping=aes(x=predicted_ve, 
                      y=observed_ve, 
                      ymin=observed.rate.lower,
                      ymax=observed.rate.upper ),
                  colour = "black" ) + 
  geom_abline(intercept=0, 
              slope=1, 
              linetype='dashed') + 
  coord_cartesian(xlim=c(0.2, 1 ),
                  ylim=c( 0.2, 1))  + 
  labs(x='Predicted Vaccine Efficacy (%)',
       y='Observed Vaccine Efficacy (%)',
       col=" ") + 
  theme_bw(11)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = "bottom") +
  geom_pointrange( aes(y=observed_VE_std, x=estimate, ymin = 0.295, ymax = 0.801, colour="Standard"), data.1_M5) +
  geom_pointrange( aes(y=observed_VE_fx, x=estimate, ymin=0.668, ymax=0.946, colour="Delayed-Fractional"), data.2_M5) + 
  scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1), labels=c(20,40,60,80,100)) + 
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1), labels=c(20,40,60,80,100)) + 
  scale_colour_manual(values=lacroix_palette("PommeBaya"))

ab_av_plot

#=====================================================================================================================
### Anitbody Titre only Model Comparison 

### parameter values 
par_abonly = apply(X=MCMC_burn_M8[,1:5], MARGIN=2, FUN=median)

n_ab <- as.numeric(par_abonly[1])
sig_n_ab <- as.numeric(par_abonly[2])
sig_mu_ab <- as.numeric(par_abonly[3])
b_nanp_ab <- as.numeric(par_abonly[4])
a_nanp_ab <- as.numeric(par_abonly[5])

p_ab <- (sig_n_ab^2-n_ab)/(sig_n_ab^2)
r_ab <- (n_ab^2)/(sig_n_ab^2-n_ab)

### vaccine efficacy point estimates 
DR_abonly <- ( 1/(1+(nanp/b_nanp_ab)^a_nanp_ab) )

VE_inf_abonly <- (r_ab/(n_ab*DR_abonly+r_ab))^r_ab
VE_inf_abonly <- 1 - (1-VE_inf_abonly)/(1-(1-p_ab)^r_ab)

### Hosmer-Lemshow results  
hl.results_M8 <- hoslem.test(protected, VE_inf_abonly, g=5 ) 
hl.results_M8 

### predictd efficacy funcation for titre only model 
model_VE_ab <- function(nanp, par_M8){
  n          <- par_M8[1]
  sig_n      <- par_M8[2]
  sig_mu     <- par_M8[3]
  b_nanp     <- par_M8[4]
  a_nanp     <- par_M8[5]
  
  ############################
  ## Secondary NB parameters
  
  p = (sig_n^2-n)/(sig_n^2)
  r = (n^2)/(sig_n^2-n)
  #############################
  ## Dose response
  
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )
  
  VE_inf <- (r/(n*DR+r))^r
  VE_inf <-  1 - (1-VE_inf)/(1-(1-p)^r)
  VE_inf <- mean(VE_inf)
  
  VE_inf
  
}

model_VE_ab = cmpfun(model_VE_ab, options=list(optimize=3))

### run for preductions - fractional arm only 
N_sam = 1000

sam_seq_ab = round(seq(from=1, to=nrow(MCMC_burn_M8), length=N_sam))

M8_sam_test_fxab = matrix(NA, nrow=N_sam, ncol=length(avidity_fx))

for(k in 1:N_sam){
  M8_sam_test_fxab[k,] = sapply(nanp_fx, model_VE_ab,  par=MCMC_burn_M8[sam_seq_ab[k],1:5])
}

M8_quant_fxab = matrix(NA, nrow=3, ncol=length(avidity_fx))

for(j in 1:length(avidity_fx)){
  M8_quant_fxab[,j] = quantile( M8_sam_test_fxab[,j], prob=c(0.025, 0.5, 0.975) )
}

### run for preductions - standard arm only 
M8_sam_test_stdab = matrix(NA, nrow=N_sam, ncol=length(avidity_std))

for(k in 1:N_sam){
  M8_sam_test_stdab[k,] = sapply(nanp_std, model_VE_ab,  par=MCMC_burn_M8[sam_seq_ab[k],1:6])
}

M8_quant_stdab = matrix(NA, nrow=3, ncol=length(avidity_std))

for(j in 1:length(avidity_std)){
  M8_quant_stdab[,j] = quantile( M8_sam_test_stdab[,j], prob=c(0.025, 0.5, 0.975) )
}

### plotting 
plot.data_M8 <- data.frame(predicted_ve=prop.table(hl.results_M8$expected,1)[,2],                      
                           observed_ve=prop.table(hl.results_M8$observed, 1)[,2] )

confint.df <- binom.confint(
  hl.results_M8$observed[,2],
  margin.table(hl.results_M8$observed, 1),conf.level = 0.95,
  methods="wilson")

plot.data_M8$observed.rate.upper <- confint.df$upper.Freq
plot.data_M8$observed.rate.lower <- confint.df $lower.Freq

data.p <- data.frame(observed_VE_std)
data.p$estimate <- median(M8_quant_stdab[2,])
data.p$lower <- median(M8_quant_stdab[1,])
data.p$higher <- median(M8_quant_stdab[3,])

data.w <- data.frame(observed_VE_fx)
data.w$estimate <- median(M8_quant_fxab[2,])
data.w$lower <- median(M8_quant_fxab[1,])
data.w$higher <- median(M8_quant_fxab[3,])

### PLOTTING ELEMENT 
ab_plot <- ggplot() + 
  geom_pointrange( aes(x=predicted_ve, 
                       y=observed_ve, 
                       ymin=observed.rate.lower,
                       ymax=observed.rate.upper ),
                   plot.data_M8,  
                   colour = "black") + 
  geom_abline(intercept=0, 
              slope=1, 
              linetype='dashed') + 
  coord_cartesian(xlim=c(0.2, 1 ),
                  ylim=c( 0.2, 1))  + 
  labs(x='Predicted Vaccine Efficacy (%)',
       y='Observed Vaccine Efficacy (%)',
       col=" ") + 
  theme_bw(11)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom") +
  geom_pointrange( aes(y=observed_VE_std, x=estimate, ymin = 0.294, ymax = 0.801, col="Standard"), data.p) +
  geom_pointrange( aes(y=observed_VE_fx, x=estimate, ymin=0.668, ymax=0.946, col="Delayed-Fractional"), data.w) +
  scale_x_continuous(breaks=c(0.2,0.4,0.6,0.8,1), labels=c(20,40,60,80,100)) + 
  scale_y_continuous(breaks=c(0.2,0.4,0.6,0.8,1), labels=c(20,40,60,80,100)) +
  scale_colour_manual(values=lacroix_palette("PommeBaya"))

ab_plot

### MULTIPLOT 
# figure_4 <- ggarrange(ab_av_plot, ab_plot , 
#                       labels = c("A", "B"),
#                       ncol = 2, nrow = 1)
# 
# annotate_figure(figure_4,
#                 bottom = text_grob("Predicted Vaccine Efficacy (%)", color = "black"),
#                 left = text_grob("Observed Vaccine Efficacy (%)", color = "black", rot = 90))

(ab_av_plot + ab_plot) / (p1/p2) + 
  plot_layout(guides = "collect", heights = c(1,1.5)) & 
  theme(legend.position = "bottom") & 
  plot_annotation(tag_levels = "A")

ggsave("Fig_4.png", width=9, height=9, dpi=600)
