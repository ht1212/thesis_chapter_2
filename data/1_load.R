#######################################################################
# Code to load data into R and create necessary objects for analysis  #
#######################################################################

MAL71 <- read.csv("MAL71_data.csv") 

############################################
### create objects in global environment ###
############################################

nanp    <-  MAL71$nanp_titre
avidity <-  MAL71$nanp_avidity * 100
vaccine <-  as.character(MAL71$vaccine)
infect  <-  MAL71$infect
time       <-  MAL71$time_to_onset

nanp_std      <- nanp[which(vaccine=="std")]
avidity_std   <- avidity[which(vaccine=="std")]
infect_std    <- infect[which(vaccine=="std")]
time_std      <- time[which(vaccine=="std")]

nanp_fx     <- nanp[which(vaccine=="fx")]
avidity_fx  <- avidity[which(vaccine=="fx")]
infect_fx   <- infect[which(vaccine=="fx")]
time_fx     <- time[which(vaccine=="fx")]

############################################ 
## Fixed Values to load into environment  ##
############################################

mu <- 30000
mult <- 3.8
pfT <- 50000000
t_L <- 6.5
Q = pfT*mult^(-(time-t_L))

max_spz <- 1000
min_DR <- 1e-6
log2 <- log(2)

NN_fx <-length(nanp_fx)
NN_std <-length(nanp_std)

NN <- NN_fx + NN_std

#############################################
## Packages to Have installed for analysis ##
#############################################

#install.packages(c('ggplot2', 'wesanderson', 'ggpubr', 'compiler', 'MASS', 'binom', 'survival', 'survminer', 'ResourceSelection', 'RColorBrewer', 'fields' , 'NumDeriv' ))

library(fields)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(compiler)
library(wesanderson) 
library(survminer)
library(survMisc)
library(numDeriv)
library(MASS)
library(binom) 
library(compiler)
library(ResourceSelection)
library(ghibli)
library(patchwork)
library(LaCroixColoR)
