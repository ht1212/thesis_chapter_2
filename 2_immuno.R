#############################################################################################################
# Code to Reproduce the Results from the section: Challenge Study Data 
# Within this section we summarise the immunological data we have been provided, produce display plots
# and the associated statistical tests
# Need to have loaded the data and installed the packages in 1.Data 
#############################################################################################################

### Immune measurements boxplots  
bp <- data.frame(vaccine, avidity, infect) 
bp$infection[which(bp$infect=="0")] <- "Protected"
bp$infection[which(bp$infect=="1")] <- "Infected"

bp$vaccine2[which(bp$vaccine=="std")] <- "Standard"
bp$vaccine2[which(bp$vaccine=="fx")] <- "Delayed-Fractional"

bp$titre <- nanp
bp$time <- time 

#### avidity plot  ####
p <- 
  ggplot(bp,aes(x=vaccine2, y=avidity, fill=infection)) + 
  geom_boxplot(outlier.colour="black", 
               outlier.shape=20,
               position=position_dodge(0.8)) +
  labs(x=" ", 
       y = "IgG Avidity Index" ,
       fill = "Infection Status") +
  theme_bw(11) +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values=lacroix_palette("Lime"))
 # scale_fill_ghibli_d(name="MarnieMedium2", direction=-1)

p

#### titre plot ####
p2 <- 
  ggplot(bp, aes(x=vaccine2, y=titre, fill=infection)) + 
  geom_boxplot(outlier.colour="black", 
               outlier.shape=20,
               position=position_dodge(0.8)) +
  labs(x=" ", 
       y = "IgG Titre (ELISA Units)" ,
       fill = "Infection Status") +
  theme_bw(11) +
  theme(legend.position = "bottom") + 
  scale_fill_manual(values=lacroix_palette("Lime"))
 # scale_fill_ghibli_d(name="MarnieMedium2", direction=-1)

p2

#### together ##### 
p2 + p + plot_annotation(tag_levels = c("A")) + plot_layout(guides = "collect")& theme(legend.position = 'bottom') 
                
ggsave("Figure_2-2.png", width=8, height=4, dpi=600)

##======================================================================
## Statistical tests

## Testing for normality of the data  
## 1. Density Plots  
ggdensity((nanp), 
          main = "Density plot of antibody titre",
          xlab = "Anti-NANP Antibody Titre")

ggdensity(avidity, 
          main = "Density plot of avidity index",
          xlab = "Anti-NANP Avidity Index")

### sharpiro-Wilk normality test for nanp 
shapiro.test(bp$titre)
# p.value 0.00903 = reject the null that data are normally distributed 

### sharpiro-wilk normality test for avidity 
shapiro.test(bp$avidity)
# p.value 0.4198 = data ~ normally distributed 

#### visual inspection of normality  
ggqqplot(log(nanp), ylab="anti-NANP titre")
ggqqplot(log(avidity), ylab="avidity index")

#### correlations between individuals measurements  
ggscatter(bp, x = "titre", y = "avidity" ,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "anti-NANP antibody titre", ylab = "avidity index")

### spearman rank correlation co-efficient  
cor <-cor.test(nanp, avidity,  method = "spearman")
cor

##=============================================================================================================
# Non-parametric statistical tests along with Figure.1 
# Due to slight non-normality of ab data 

wilcox.test(nanp_fx, nanp_std)
wilcox.test(avidity_fx, avidity_std)

wilcox.test(bp$avidity[which(bp$vaccine2=="Delayed-Fractional" & bp$infection == "Protected")], bp$avidity[which(bp$vaccine2=="Delayed-Fractional" & bp$infection == "Infected")])
wilcox.test(bp$titre[which(bp$vaccine2=="Delayed-Fractional" & bp$infection == "Protected")], bp$tire[which(bp$vaccine2=="Delayed-Fractional" & bp$infection == "Infected")])

wilcox.test(bp$avidity[which(bp$vaccine2=="Standard" & bp$infection == "Protected")], bp$avidity[which(bp$vaccine2=="Standard" & bp$infection == "Infected")])
wilcox.test(bp$titre[which(bp$vaccine2=="Standard" & bp$infection == "Protected")], bp$titre[which(bp$vaccine2=="Standard" & bp$infection == "Infected")])

#==============================================================================================================
## adding plot for combined figure 1 
obs <- survfit(Surv(time, infect) ~ vaccine2, data=bp)

surv_raw <- 
  ggsurvplot(obs,
             #risk.table = TRUE,  
             break.time.by = 7 ,
            # conf.int = T,
             legend.labs=c("Delayed-Fractional", "Standard"),
            # tables.theme = theme_cleantable(),       # Clean risk table
            censor = FALSE,                         # Remove censor points
            palette =  lacroix_palette("PommeBaya"),
            fun = "event" ,
            xlab="Time from challenge (days)", 
            ylab = "Proportion infected",
            xlim = c(0,28), 
            linetype = c(1,1),
            legend.title="Vaccination Schedule", 
            legend="bottom", 
           ggtheme = theme_bw(11)) 

surv_raw <- surv_raw$plot 

(p2 + p + plot_layout(guides = "collect")& theme(legend.position = 'bottom') ) / surv_raw + plot_annotation(tag_levels = c("A")) 

ggsave("Figure_2-2.png", width=8, height=7, dpi=600)


