
## Statistical Model checks  

## linear regression models of the time to onset of parasitemia agaist antibody titre, avidity and vaccine arm  
id <- which(infect==1)

time <- T
tyme <- vector(length=46)
tyme <- rep(NA, 46)
tyme[id] <- time[id]

plot(avidity, tyme,     
     xlab = "avidity", ylab = "time",
     pch = 19, frame = FALSE)
abline(lm(tyme ~ avidity), col = "blue")

cor.test(avidity, tyme,  method = "spearman")

#linear regression models 
summary(lm( tyme ~ log(nanp) + avidity ) )

summary(lm( tyme ~ log(nanp) + avidity + vaccine) )

summary(lm( time ~ log(nanp) + avidity + log(nanp)*avidity) ) 


# logistic regression models of protection status vs titre, avidity and vaccine arm 
protect <- vector(length = 46)
protect <- rep(1, 46)  
protect[id] <- 0 

summary(glm( protect ~ log(nanp) + avidity, family=binomial("logit") ))

summary(glm( protect ~ log(nanp) + avidity + vaccine, family=binomial("logit") ))

summary(glm( protect ~ log(nanp) + avidity + log(nanp)*avidity, family=binomial("logit") ))





