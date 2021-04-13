##############################################################################################################################
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
###############################################################################################################################

### make objects for each parameter if not in enviroment already ###
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

## SECTION 2    
## Vaccine Efficacy Heatmap - Figure.5 

## sample iteration number  
NN <- 2000 

## sequences to compute VE  ###
nanp_seq <- exp(seq(from=log(500), to=log(2.0*max(nanp)), length=NN))
av_seq <- exp(seq(from=log(0.8), to=log(100), length=NN))

## Vaccine Efficacy Function ###
VE <- function(nanp, avidity){
  DR <- ( 1/(1+(nanp/b_nanp)^a_nanp) )*exp(-log(2)*avidity/b_av)
  VE_inf <- (r/(n*DR+r))^r
  VE_inf <- 1 - (1-VE_inf)/(1-(1-p)^r)
  
  VE_inf
}

### matrix to store VE prediction ###
VE_mat <- matrix(NA, nrow=NN, ncol=NN)

for(i in 1:nrow(VE_mat)){
  VE_mat[NN-i+1,] <- VE(nanp_seq, av_seq[NN-i+1])
}	

Ncol = 10000
pp <- 0.5
colours = colorRampPalette(brewer.pal(9,"YlOrRd"))(10000)

col.sample <- seq(from=1, to=Ncol^(1/pp), length=floor(sqrt(Ncol)) )
col.sample <- floor(col.sample^pp)

colours <- colours[col.sample]

par(oma=c(1,2,1,2)) 
par(mar=c(4,4.5,1,1) + 0.1)

#png("Figure_5.png", res=600)

image.plot(x=nanp_seq, y=av_seq, z=t(VE_mat), col=colours, log="x",
           xlab=expression(paste( " Anti-NANP Antibody Titre (ELISA Units) ")), ylab="Avidity Index",
           legend.lab="Vaccine Efficacy", xaxt='n', cex.axis=0.8, 
           axis.args=list(at=c(0.2,0.4,0.6,0.8,1),
                          labels=c("20%", "40%", "60%" , "80%", "100%"), cex.axis=0.7), las=1)
axis(1, at=c(500,1000,2000,5000,10000,20000,50000,100000), label=c("500","1000", "2000", "5000", 
                                                               "10000", "20000", "50000", "100000"), cex.axis=0.8)

### calculate interquartile ranges  ###
nanp_med <- quantile(nanp, prob=0.5)
nanp_low <- quantile(nanp, prob=0.1)
nanp_high <- quantile(nanp, prob=0.9)

av_med  <- quantile(avidity, prob=0.5)
av_low  <- quantile(avidity, prob=0.1)
av_high <- quantile(avidity, prob=0.9)

points(x=rep(nanp_med,2), y=c(0.001, 1000000), 
       type='l', lty="longdash",  col="grey39", cex=0.5)
points(x=rep(nanp_low,2), y=c(0.001, 1000000), 
       type='l', lty="longdash",  col="grey39", cex=0.5)
points(x=rep(nanp_high,2), y=c(0.001, 1000000), 
       type='l', lty="longdash",  col="grey39", cex=0.5)

points(x=c(0.001, 1000000), y=rep(av_med,2),
       type='l', lty="longdash",  col="grey39", cex=0.5)
points(x=c(0.001, 1000000), y=rep(av_low,2),
       type='l', lty="longdash",  col="grey39", cex=0.5)
points(x=c(0.001, 1000000), y=rep(av_high,2), 
       type='l', lty="longdash",  col="grey39", cex=0.5)

### data points ###
infect_col <- rep("white", length(infect_fx))

# infected fx 
infect_col[which(infect_fx==1)] <- "greenyellow"

infect_col2 <- rep("white", length(infect_std))

infect_col2[which(infect_std==1)] <- "greenyellow"

points(x=nanp_fx, y=avidity_fx, pch="*",  col=infect_col)
points(x=nanp_std, y=avidity_std, pch=19, cex=0.7, col=infect_col2)


### Calculate isoclines ###
av_long <- exp(seq(from=log(0.001), to=log(100), length=100000))

av_30 <- rep(NA, NN)
av_50 <- rep(NA, NN)
av_70 <- rep(NA, NN)
av_90 <- rep(NA, NN)

for(i in 1:NN){
  av_30[i] <- av_long[ which.min( (VE(nanp_seq[i], av_long) - 0.3)^2 ) ] 	
  av_50[i] <- av_long[ which.min( (VE(nanp_seq[i], av_long) - 0.5)^2 ) ]
  av_70[i] <- av_long[ which.min( (VE(nanp_seq[i], av_long) - 0.7)^2 ) ]
  av_90[i] <- av_long[ which.min( (VE(nanp_seq[i], av_long) - 0.9)^2 ) ]
}

points(x=nanp_seq, y=av_30, type='l', cex=0.5)
points(x=nanp_seq, y=av_50, type='l', cex=0.5)
points(x=nanp_seq, y=av_70, type='l', cex=0.5)
points(x=nanp_seq, y=av_90, type='l', cex=0.5)

### add labels for isoklines ###

text(y=3+av_30[1], x=750, labels="30%", cex=0.7)
text(y=3+av_50[1], x=750, labels="50%", cex=0.7) 
text(y=2+av_70[1], x=750, labels="70%", cex=0.7)
text(y=2+av_90[1], x=750, labels="90%", cex=0.7) 

dev.off()

##----------------------------------------------------------------------------------------------------------------------
## SECTION 3 
## VE dose response graphs for each immune measurement  

### extract posterior medians and calculate VE for each increasing immune measurement with 
### the other antibody characteristic held constant (have chosen the median value here, can also set at an arbitary value) ###

##=====================================================================================================================
## 1 VE for each unit increase in Avidity with titre held constant 

par(mfrow=c(2,2))

#############################################
## Extract posterior medians and calculate model prediction
avidity_seq = avidity_seq <- seq(from=0, to=100, length=200)

############################
## VE calculation 
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
  
  DR <-  exp( -log(2)*(avidity/b_av) ) * ( 1/(1+(nanp/b_nanp)^a_nanp) )
  
  #############################
  ## VE against infection 
  VE_inf <- (r/(n*DR+r))^r
  VE_inf <-  1 - (1-VE_inf)/(1-(1-p)^r)
  
  
  VE_inf
  
}

model_VE_s <- function(nanp, avidity, par_M5){
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
  
  DR <-  exp( -log(2)*(avidity/b_av) ) * ( 1/(1+(nanp/b_nanp)^a_nanp) )
  
  #############################
  ## VE against infection 
  # VE_inf <- (r/(n*DR+r))^r
  VE_inf <-  1 - DR #(1-VE_inf)/(1-(1-p)^r)
  
  
  VE_inf
  
}

M5_par_median_pred_av = sapply(avidity_seq, model_VE, nanp=median_nanp, par=par_median)
M5_par_median_pred_av_s = sapply(avidity_seq, model_VE_s, nanp=median_nanp, par=par_median)

#############################################
## Posterior prediction intervals 
median_nanp <-  3612#median(nanp) 
sporo_median <-  3612

N_sam = 10000
sam_seq = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M2_sam_av = matrix(NA, nrow=N_sam, ncol=length(avidity_seq))
M2_sam_av_s = matrix(NA, nrow=N_sam, ncol=length(avidity_seq))

for(k in 1:N_sam){
  M2_sam_av[k,] = sapply(avidity_seq, model_VE, nanp=median_nanp, par=MCMC_burn_M5[sam_seq[k],1:6])
}
for(k in 1:N_sam){
  M2_sam_av_s[k,] = sapply(avidity_seq, model_VE_s, nanp=sporo_median, par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_av = matrix(NA, nrow=3, ncol=length(avidity_seq))
M2_quant_av_s = matrix(NA, nrow=3, ncol=length(avidity_seq))

for(j in 1:length(avidity_seq)){
  M2_quant_av[,j] = quantile( M2_sam_av[,j], prob=c(0.025, 0.5, 0.975) )
}
for(j in 1:length(avidity_seq)){
  M2_quant_av_s[,j] = quantile( M2_sam_av_s[,j], prob=c(0.025, 0.5, 0.975) )
}

##########################################################
## Construct a histogram of observed avidity measurements
bars <- 10
breaks <- exp( seq( from=log(min(0.99*avidity)), to=log(max(avidity)), length=bars+1) )
counts <- rep(0, bars)
for(i in 1:length(avidity)){
  for(j in 1:bars){
    
    if( (avidity[i]>breaks[j]) && (avidity[i]<=breaks[j+1]) ){
      counts[j] <- counts[j] + 1 
    }
  }
}
counts <- 0.9*counts/max(counts)

#####################################################################
## plot infection blocking efficacy 1 and sporo blocking efficacy 2 
plot(x=avidity_seq, y=M2_quant_av[2,], type='l',  
     ylim=c(0,1.05), xlim=c(1, max(avidity_seq)), yaxt='n', xaxt='n',
     xlab=expression(paste( "Avidity Index")),
     ylab="Vaccine Efficacy (%)") 
#title("A" , adj  =)

axis(1, at=c(1,20,50,70,100), label=c(1,20,50,70,100)) #cex.axis=axis.size)
axis(2, at=c(0.0,0.2, 0.4, 0.6, 0.8, 1.0), label=c(0, 20, 40, 60, 80, 100), las=1)

for(j in 1:bars){
  polygon(x=c(breaks[j:(j+1)], rev(breaks[j:(j+1)]) ),
          y=c(0, 0, rep(counts[j], 2) ),
          col=rgb(190/256,190/256,190/256,0.4) )
}		

points(x=avidity_seq, y=M2_quant_av[2,], 
       type='l', lwd=3, col="darkorange") 

polygon(x=c(avidity_seq, rev(avidity_seq)), 
        y=c( M2_quant_av[1,], rev(M2_quant_av[3,]) ),
        col=rgb(255/256,165/256,0/256,0.4), border=NA)

#-----------------efficacy per sporzoite---------------------------------
plot(x=avidity_seq, y=M2_quant_av_s[2,], type='l',  
     ylim=c(0,1.05), xlim=c(1, max(avidity_seq)), yaxt='n', xaxt='n',
     xlab=expression(paste( "Avidity Index")),
     ylab="Vaccine Efficacy (%)") 
axis(1, at=c(1,20,50,70,100), label=c(1,20,50,70,100)) #cex.axis=axis.size)
axis(2, at=c(0.0,0.2, 0.4, 0.6, 0.8, 1.0), label=c(0, 20, 40, 60, 80, 100), las=1)

for(j in 1:bars){
  polygon(x=c(breaks[j:(j+1)], rev(breaks[j:(j+1)]) ),
          y=c(0, 0, rep(counts[j], 2) ),
          col=rgb(190/256,190/256,190/256,0.4) )
}	

points(x=avidity_seq, y=M2_quant_av_s[2,], 
       type='l', lwd=3, col="#1BB6AF") 

polygon(x=c(avidity_seq, rev(avidity_seq)), 
        y=c( M2_quant_av_s[1,], rev(M2_quant_av_s[3,]) ),
        col=rgb(27/256, 182/256, 175/256, 0.4), border=NA)



###===========================================================================================================================
## 2 VE for each unit increase in anti-nanp antibody titre
##

## best fitting median parameter values  
par_median = apply(X=MCMC_burn_M5[,1:6], MARGIN=2, FUN=median)

## nanp_sequence  
nanp_seq = exp( seq(from=log(1), to=log(1.5*max(nanp)), length=200) )

############################
## VE efficacy model to calculate  
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
  
  #####################################
  ## Vaccine Efficacy Against Infection 
  
  VE_inf <- (r/(n*DR+r))^r
  VE_inf <-  1 -(1-VE_inf)/(1-(1-p)^r)
  
  
  VE_inf
  
}

model_VE_s <- function(nanp, avidity, par_M5){
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
  
  #####################################
  ## Vaccine Efficacy Against Infection 
  
 # VE_inf <- (r/(n*DR+r))^r
  VE_inf <-  1 -DR #(1-VE_inf)/(1-(1-p)^r)
  
  
  VE_inf
  
}

## prediction ##
median_av <- 8.7#median(avidity) 
av_sporo  <- 8.7

M5_par_median_pred = sapply(nanp_seq, model_VE, avidity=median_av, par=par_median)
M5_par_median_pred_s = sapply(nanp_seq, model_VE_s, avidity=av_sporo, par=par_median)


### posterior prediction intervals ###
N_sam = 10000
sam_seq = round(seq(from=1, to=nrow(MCMC_burn_M5), length=N_sam))

M2_sam_cs = matrix(NA, nrow=N_sam, ncol=length(nanp_seq))
M2_sam_cs_s = matrix(NA, nrow=N_sam, ncol=length(nanp_seq))

for(k in 1:N_sam){
  M2_sam_cs[k,] = sapply(nanp_seq, model_VE, avidity=median_av, par=MCMC_burn_M5[sam_seq[k],1:6])
}
for(k in 1:N_sam){
  M2_sam_cs_s[k,] = sapply(nanp_seq, model_VE_s, avidity=av_sporo, par=MCMC_burn_M5[sam_seq[k],1:6])
}

M2_quant_cs = matrix(NA, nrow=3, ncol=length(nanp_seq))
M2_quant_cs_s = matrix(NA, nrow=3, ncol=length(nanp_seq))

for(j in 1:length(nanp_seq)){
  M2_quant_cs[,j] = quantile( M2_sam_cs[,j], prob=c(0.025, 0.5, 0.975) )
}
for(j in 1:length(nanp_seq)){
  M2_quant_cs_s[,j] = quantile( M2_sam_cs_s[,j], prob=c(0.025, 0.5, 0.975) )
}

###############################
## Plotting 

## construct a histogram out of polygons ##
bars <- 10
breaks <- exp( seq( from=log(min(0.99*nanp)), to=log(max(nanp)), length=bars+1) )
counts <- rep(0, bars)
for(i in 1:length(nanp)){
  for(j in 1:bars){
    
    if( (nanp[i]>breaks[j]) && (nanp[i]<=breaks[j+1]) ){
      counts[j] <- counts[j] + 1 
    }
  }
}
counts <- 0.9*counts/max(counts)

### plot ###
plot(x=nanp_seq, y=M2_quant_cs[2,], type='l', log="x",
     ylim=c(0,1.05), xlim=c(500, max(nanp_seq)), yaxt = "n", xaxt='n',
     xlab=expression(paste( "Total IgG Titre (ELISA Units)")),
     ylab="Vaccine Efficacy (%)", las=1)
axis(1, at=c(1000,  10000,  100000), label=c("1,000", "10,000", "100,000")) 
axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), label=c(0, 20, 40, 60,80,100), las=1)

for(j in 1:bars){
  polygon(x=c(breaks[j:(j+1)], rev(breaks[j:(j+1)]) ),
          y=c(0, 0, rep(counts[j], 2) ),
          col=rgb(190/256,190/256,190/256,0.4) )
}
points(x=nanp_seq, y=M2_quant_cs[2,], 
       type='l', lwd=3, col="darkorange")
polygon(x=c(nanp_seq, rev(nanp_seq)), 
        y=c( M2_quant_cs[1,], rev(M2_quant_cs[3,]) ),
        col=rgb(255/256,165/256,0/256,0.4), border=NA)



plot(x=nanp_seq, y=M2_quant_cs_s[2,], type='l', log="x",
     ylim=c(0,1.05), xlim=c(500, max(nanp_seq)), yaxt = "n", xaxt='n',
     xlab=expression(paste( "Total IgG Titre (ELISA Units)")),
     ylab="Vaccine Efficacy (%)", las=1)
axis(1, at=c(1000,  10000,  100000), label=c("1,000", "10,000", "100,000")) 
axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), label=c(0, 20, 40, 60,80,100), las=1)

for(j in 1:bars){
  polygon(x=c(breaks[j:(j+1)], rev(breaks[j:(j+1)]) ),
          y=c(0, 0, rep(counts[j], 2) ),
          col=rgb(190/256,190/256,190/256,0.4) )
}
polygon(x=c(nanp_seq, rev(nanp_seq)), 
        y=c( M2_quant_cs_s[1,], rev(M2_quant_cs_s[3,]) ),
        col=rgb(27/256, 182/256, 175/256, 0.4), border=NA)

points(x=nanp_seq, y=M2_quant_cs_s[2,], 
       type='l', lwd=3, col="#1BB6AF")

# legend(x="topleft", 
#        legend = c("efficacy against infection" ), 
#        col = c("darkorange"),
#        lwd = 3,  bty='n', lty=c(1,1) )




