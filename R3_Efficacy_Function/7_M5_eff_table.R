######################################################################################################################################
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
## Code to generate the results from the section: Model Predicted Vaccine Efficacy as a Function of the Antibody Response  
## Section 1: Recreates Table.3 and the Efficacy Distributions 
## Section 2: Recreates Figure.5 - heatmap of predicted VE  
## Section 3: Recreates Figure.6 - VE dose-response 
#######################################################################################################################################

################################
##        SECTION 1           ##
################################


### make objects for each parameter ###

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

### 1 VE against Infection across the whole study ###

### VE calculation ###
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

N_sam = 20000

sam_seq = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M2_sam_test = matrix(NA, nrow=N_sam, ncol=length(avidity))

for(k in 1:N_sam){
  M2_sam_test[k,] = sapply(avidity, model_VE, nanp=nanp, par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_test = matrix(NA, nrow=3, ncol=length(avidity))

for(j in 1:length(avidity)){
  M2_quant_test[,j] = quantile( M2_sam_test[,j], prob=c(0.025, 0.5, 0.975) )
}

########### CIs and Point Estimate of VE across the study (all Volunteers)
# lower CI
mean(M2_quant_test[1,])
#value
mean(M2_quant_test[2,])
#upper CI
mean(M2_quant_test[3,])

##########################################################################################################
##########################################################################################################
### 2. VE stratified by terciles of immune measurements 

nanp_low  <- quantile(nanp, prob=0.333)
nanp_high <- quantile(nanp, prob=0.667)

avidity_low  <- quantile(avidity, prob=0.333)
avidity_high <- quantile(avidity, prob=0.667)

nanp_1 <- which(nanp <= nanp_low)
nanp_2 <- intersect( which(nanp>nanp_low), which(nanp<=nanp_high) )
nanp_3 <- which(nanp > nanp_high)

avidity_1 <- which(avidity <= avidity_low)
avidity_2 <- intersect( which(avidity>avidity_low), which(avidity<=avidity_high) )
avidity_3 <- which(avidity > avidity_high)

### storage matrix ###
VE_pred <- matrix(NA, nrow=4, ncol=4)
VE_obs <- matrix(NA, nrow=4, ncol=4)
count <- matrix(NA, nrow=3, ncol=3)

### naming matrix cols and rows ####
colnames(VE_pred) <- c("nanp_1", "nanp_2", "nanp_3", "mean")
rownames(VE_pred) <- c("avidity_1", "avidity_2", "avidity_3", "mean")

colnames(VE_obs) <- c("nanp_1", "nanp_2", "nanp_3", "mean")
rownames(VE_obs) <- c("avidity_1", "avidity_2", "avidity_3", "mean")

colnames(count) <- c("nanp_1", "nanp_2", "nanp_3")
rownames(count) <- c("avidity_1", "avidity_2", "avidity_3")

### identifiers ###
# All avidity low 
# titre low
C_11 <- intersect( avidity_1, nanp_1 )
# titre medium 
C_12 <- intersect( avidity_1, nanp_2 )
# titre high 
C_13 <- intersect( avidity_1, nanp_3 )

# All avidity medium 
# titre low
C_21 <- intersect( avidity_2, nanp_1 )
# titre medium 
C_22 <- intersect( avidity_2, nanp_2 )
# titre high 
C_23 <- intersect( avidity_2, nanp_3 )

# All avidity High 
# titre low 
C_31 <- intersect( avidity_3, nanp_1 )
# titre medium 
C_32 <- intersect( avidity_3, nanp_2 )
# titre high 
C_33 <- intersect( avidity_3, nanp_3 )

vaccine[C_33]
vaccine[C_32]
vaccine[C_31]
vaccine[C_33]

### counts ###
count[1,1] <- length(C_11)
count[1,2] <- length(C_12)
count[1,3] <- length(C_13)

count[2,1] <- length(C_21)
count[2,2] <- length(C_22)
count[2,3] <- length(C_23)

count[3,1] <- length(C_31)
count[3,2] <- length(C_32)
count[3,3] <- length(C_33)

### observed VE ###
VE_obs[1,1]  <- 1 - sum(infect[C_11])/length(infect[C_11]) 
VE_obs[1,2]  <- 1 - sum(infect[C_12])/length(infect[C_12])
VE_obs[1,3]  <- 1 - sum(infect[C_13])/length(infect[C_13])

VE_obs[2,1]  <- 1 - sum(infect[C_21])/length(infect[C_21])  
VE_obs[2,2]  <- 1 - sum(infect[C_22])/length(infect[C_22])  
VE_obs[2,3]  <- 1 - sum(infect[C_23])/length(infect[C_23])  

VE_obs[3,1]  <- 1 - sum(infect[C_31])/length(infect[C_31]) 
VE_obs[3,2]  <- 1 - sum(infect[C_32])/length(infect[C_32])  
VE_obs[3,3]  <- 1 - sum(infect[C_33])/length(infect[C_33])  

VE_obs[1,4]  <- 1 - sum(infect[avidity_1])/length(infect[avidity_1])  
VE_obs[2,4]  <- 1 - sum(infect[avidity_2])/length(infect[avidity_2])
VE_obs[3,4]  <- 1 - sum(infect[avidity_3])/length(infect[avidity_3])

VE_obs[4,1]  <- 1 - sum(infect[nanp_1])/length(infect[nanp_1])  
VE_obs[4,2]  <- 1 - sum(infect[nanp_2])/length(infect[nanp_2])
VE_obs[4,3]  <- 1 - sum(infect[nanp_3])/length(infect[nanp_3])

VE_obs[4,4]  <- 1 - sum(infect)/length(infect)

#######################
#### Predicted VE #####

######################################
#### C_11 avidity low titre low
N_sam = 20000

sam_seq = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M2_sam_ll = matrix(NA, nrow=N_sam, ncol=length(avidity[C_11]))

for(k in 1:N_sam){
  M2_sam_ll[k,] = sapply(avidity[C_11], model_VE, nanp=nanp[C_11], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_ll = matrix(NA, nrow=3, ncol=length(avidity[C_11]))


for(j in 1:length(avidity[C_11])){
  M2_quant_ll[,j] = quantile( M2_sam_ll[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_ll[1,])
mean(M2_quant_ll[2,])
mean(M2_quant_ll[3,])

#####################################
#### C_12 aidity low titre medium 
M2_sam_12 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_12]))

for(k in 1:N_sam){
  M2_sam_12[k,] = sapply(avidity[C_12], model_VE, nanp=nanp[C_12], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_12 = matrix(NA, nrow=3, ncol=length(avidity[C_12]))


for(j in 1:length(avidity[C_12])){
  M2_quant_12[,j] = quantile( M2_sam_12[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_12[1,])
mean(M2_quant_12[2,])
mean(M2_quant_12[3,])


#####################################
#### C_13 avidity low titre high

M2_sam_13 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_13]))

for(k in 1:N_sam){
  M2_sam_13[k,] = sapply(avidity[C_13], model_VE, nanp=nanp[C_13], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_13 = matrix(NA, nrow=3, ncol=length(avidity[C_13]))


for(j in 1:length(avidity[C_13])){
  M2_quant_13[,j] = quantile( M2_sam_13[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_13[1,])
mean(M2_quant_13[2,])
mean(M2_quant_13[3,])


####################################
## C_21 avidity medium titre low 

M2_sam_21 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_21]))

for(k in 1:N_sam){
  M2_sam_21[k,] = sapply(avidity[C_21], model_VE, nanp=nanp[C_21], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_21 = matrix(NA, nrow=3, ncol=length(avidity[C_21]))


for(j in 1:length(avidity[C_21])){
  M2_quant_21[,j] = quantile( M2_sam_21[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_21[1,])
mean(M2_quant_21[2,])
mean(M2_quant_21[3,])

#####################################
## C_22 avidity medium titre medium 

M2_sam_22 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_22]))

for(k in 1:N_sam){
  M2_sam_22[k,] = sapply(avidity[C_22], model_VE, nanp=nanp[C_22], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_22 = matrix(NA, nrow=3, ncol=length(avidity[C_22]))


for(j in 1:length(avidity[C_22])){
  M2_quant_22[,j] = quantile( M2_sam_22[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_22[1,])
mean(M2_quant_22[2,])
mean(M2_quant_22[3,])

#####################################
## C_23 avidity medium titre high  

M2_sam_23 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_23]))

for(k in 1:N_sam){
  M2_sam_23[k,] = sapply(avidity[C_23], model_VE, nanp=nanp[C_23], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_23 = matrix(NA, nrow=3, ncol=length(avidity[C_23]))


for(j in 1:length(avidity[C_23])){
  M2_quant_23[,j] = quantile( M2_sam_23[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_23[1,])
mean(M2_quant_23[2,])
mean(M2_quant_23[3,])

#######################################
## C_31 avidity high titre low 

M2_sam_31 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_31]))

for(k in 1:N_sam){
  M2_sam_31[k,] = sapply(avidity[C_31], model_VE, nanp=nanp[C_31], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_31 = matrix(NA, nrow=3, ncol=length(avidity[C_31]))


for(j in 1:length(avidity[C_31])){
  M2_quant_31[,j] = quantile( M2_sam_31[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_31[1,])
mean(M2_quant_31[2,])
mean(M2_quant_31[3,])


########################################
## C_32 avidity high titre medium 

M2_sam_32 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_32]))

for(k in 1:N_sam){
  M2_sam_32[k,] = sapply(avidity[C_32], model_VE, nanp=nanp[C_32], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_32 = matrix(NA, nrow=3, ncol=length(avidity[C_32]))


for(j in 1:length(avidity[C_32])){
  M2_quant_32[,j] = quantile( M2_sam_32[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_32[1,])
mean(M2_quant_32[2,])
mean(M2_quant_32[3,])

###########################################
## C_33 avidity high titre high 

M2_sam_33 = matrix(NA, nrow=N_sam, ncol=length(avidity[C_33]))

for(k in 1:N_sam){
  M2_sam_33[k,] = sapply(avidity[C_33], model_VE, nanp=nanp[C_33], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_33 = matrix(NA, nrow=3, ncol=length(avidity[C_33]))


for(j in 1:length(avidity[C_33])){
  M2_quant_33[,j] = quantile( M2_sam_33[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_33[1,])
mean(M2_quant_33[2,])
mean(M2_quant_33[3,])

########################################
## 14 low avidity across all titres 

M2_sam_14 = matrix(NA, nrow=N_sam, ncol=length(avidity_1))

for(k in 1:N_sam){
  M2_sam_14[k,] = sapply(avidity[avidity_1], model_VE, nanp=nanp[avidity_1], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_14 = matrix(NA, nrow=3, ncol=length(avidity_1))


for(j in 1:length(avidity_1)){
  M2_quant_14[,j] = quantile( M2_sam_14[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_14[1,])
mean(M2_quant_14[2,])
mean(M2_quant_14[3,])

########################################
## 24 avidity medium acorss all titres 

M2_sam_24 = matrix(NA, nrow=N_sam, ncol=length(avidity_2))

for(k in 1:N_sam){
  M2_sam_24[k,] = sapply(avidity[avidity_2], model_VE, nanp=nanp[avidity_2], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_24 = matrix(NA, nrow=3, ncol=length(avidity_2))


for(j in 1:length(avidity_2)){
  M2_quant_24[,j] = quantile( M2_sam_24[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_24[1,])
mean(M2_quant_24[2,])
mean(M2_quant_24[3,])

########################################
## 34 avidity high across all titres 

M2_sam_34 = matrix(NA, nrow=N_sam, ncol=length(avidity_3))

for(k in 1:N_sam){
  M2_sam_34[k,] = sapply(avidity[avidity_3], model_VE, nanp=nanp[avidity_3], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_34 = matrix(NA, nrow=3, ncol=length(avidity_3))


for(j in 1:length(avidity_3)){
  M2_quant_34[,j] = quantile( M2_sam_34[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_34[1,])
mean(M2_quant_34[2,])
mean(M2_quant_34[3,])

########################################
## 41 titre low across all avidities 

M2_sam_41 = matrix(NA, nrow=N_sam, ncol=length(nanp_1))

for(k in 1:N_sam){
  M2_sam_41[k,] = sapply(nanp[nanp_1], model_VE, avidity=avidity[nanp_1], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_41 = matrix(NA, nrow=3, ncol=length(nanp_1))


for(j in 1:length(nanp_1)){
  M2_quant_41[,j] = quantile( M2_sam_41[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_41[1,])
mean(M2_quant_41[2,])
mean(M2_quant_41[3,])

########################################
## 42 titre medium across all avidities 

M2_sam_42 = matrix(NA, nrow=N_sam, ncol=length(nanp_2))

for(k in 1:N_sam){
  M2_sam_42[k,] = sapply(nanp[nanp_2], model_VE, avidity=avidity[nanp_2], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_42 = matrix(NA, nrow=3, ncol=length(nanp_2))


for(j in 1:length(nanp_2)){
  M2_quant_42[,j] = quantile( M2_sam_42[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_42[1,])
mean(M2_quant_42[2,])
mean(M2_quant_42[3,])

########################################
## 43 titre high across all titres 

M2_sam_43 = matrix(NA, nrow=N_sam, ncol=length(nanp_3))

for(k in 1:N_sam){
  M2_sam_43[k,] = sapply(nanp[nanp_3], model_VE, avidity=avidity[nanp_3], par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_43 = matrix(NA, nrow=3, ncol=length(nanp_3))


for(j in 1:length(nanp_3)){
  M2_quant_43[,j] = quantile( M2_sam_43[,j], prob=c(0.025, 0.5, 0.975) )
}

mean(M2_quant_43[1,])
mean(M2_quant_43[2,])
mean(M2_quant_43[3,])

#################################################################################################################
################## put values into data frame = predicted values for Table.3 

### Predicted VE ###
VE_pred[1,1]  <- mean(M2_quant_ll[2,])
VE_pred[1,2]  <- mean(M2_quant_12[2,])
VE_pred[1,3]  <- mean(M2_quant_13[2,]) 

VE_pred[2,1]  <-mean(M2_quant_21[2,]) 
VE_pred[2,2]  <- mean(M2_quant_22[2,])  
VE_pred[2,3]  <- mean(M2_quant_23[2,]) 

VE_pred[3,1]  <- mean(M2_quant_31[2,])
VE_pred[3,2]  <- mean(M2_quant_32[2,])  
VE_pred[3,3]  <- mean(M2_quant_33[2,])

VE_pred[1,4]  <- mean(M2_quant_14[2,]) 
VE_pred[2,4]  <- mean(M2_quant_24[2,])
VE_pred[3,4]  <- mean(M2_quant_34[2,])

VE_pred[4,1]  <- mean(M2_quant_41[2,]) 
VE_pred[4,2]  <- mean(M2_quant_42[2,])
VE_pred[4,3]  <- mean(M2_quant_43[2,])

VE_pred[4,4]  <- mean(M2_quant_test[2,])

# view data frame for Table.3 predicted values 


########################################
## Histograms of individual level VE  ##
########################################
DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )* exp( -log(2)*(avidity/b_av) )

VE_inf <- (r/(n*DR+r))^r
VE_inf <- 1 - (1-VE_inf)/(1-(1-p)^r)

hist <- data.frame(VE_inf, vaccine)
infect_fx   <- infect[which(vaccine=="fx")]

hist$VE_inf <- hist$VE_inf*100

ggplot(hist, aes(x=VE_inf, fill=vaccine)) + 
  facet_wrap( ~ vaccine ) +
  geom_histogram(breaks=seq(20, 100, by=10), alpha=0.8, color="black" ,position="dodge" ) + 
  scale_x_continuous(name = "Vaccine Efficacy (%)") + 
  scale_fill_manual(values =lacroix_palette("PommeBaya"),name="Vaccine Schedule" , labels=c("Delayed-Fractional", "Standard") ) +
  theme_bw(11) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    panel.grid = element_blank(), 
    legend.position = "bottom"
  )

ggsave("Fig_6.png", width=8,height=4,dpi=600)
