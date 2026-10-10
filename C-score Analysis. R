#########Co-occurence tests##############
library(tidyverse)
library(ggeffects)
library(ggpubr)
library(viridis)
library(bipartite)
library(ecospat)
library(HuraultMisc)
library(cowplot)
library(lavaan)
library(semPlot)
library(semptools)
library(paletteer)
library(OpenMx)

#########################
##Checkerboard analysis##
#########################

theme_set(theme_pubr())

coi2017.dat <- read.csv(file="Coinfection_data_2017.csv")
env.dat <- read.csv(file="SummarizedData_2017.csv")
coi2017.dat <- left_join(coi2017.dat, env.dat, by="Lake")
coi2017.dat$Year <- "2017"


coi2019.dat <- read.csv(file="Parasite_Data_2019.csv")
coi2019.dat$Year <- "2019"
env.dat <- read.csv(file="SummarizedData_2019.csv")
#coi2019.dat <- left_join(coi2019.dat, env.dat, by="Lake")
par.data.2019 <- merge(coi2019.dat, env.dat, by="Lake")


ISPO.nb <- MASS::glm.nb(Gregarines ~ 1, data=filter(coi2019.dat, Species=="ISPO"))
#theta=0.08127322236
#b0=0.4488 

cbbPalette1 <- c("#000000", "#56B4E9","#009E73", "#D55E00")
cbbPalette2 <- c("#56B4E9","#009E73", "#D55E00")

#add variable for mite scars
par.data.2019$holes.plus.mites<-par.data.2019$final_mite+par.data.2019$total_holes

#add variable for ln(mites) and ln(greg)
par.data.2019$ln.mites<-log(par.data.2019$holes.plus.mites+1)
par.data.2019$ln.greg<-log(par.data.2019$Gregarines+1)

#add variable for coinfection
par.data.2019<- par.data.2019%>%
  mutate(
    coinfected=case_when(
      Mite_Infected == 1 & Greg_Infected == 1 ~ 1,
      Mite_Infected == 0 & Greg_Infected == 1 ~ 0,
      Mite_Infected == 1 & Greg_Infected == 0 ~ 0,
      Mite_Infected == 0 & Greg_Infected == 0 ~ 0,
    )
  )



####simulations
MASS::glm.nb(Gregarines ~ total_holes + I(total_holes^2), data=filter(coi2019.dat, Species=="ISPO"))

x <- seq(0, 30)
b0 = 0.24152 #9.210e-01
b1 = 0.18833
b2 = -0.00684

None <- data.frame(x=x, model=exp(b0), Type="No association")
Linear <- data.frame(x=x, model=exp(b0 + x*0.25), Type="Linear")
Weak.Quad <- data.frame(x=x, model=exp(b0 + x*0.25 + x^2*(-0.001)), Type="Weak quadratic")
Strong.Quad <- data.frame(x=x, model=exp(b0 + x*0.25 + x^2*(-0.01)), Type="Strong quadratic")

model.df <- rbind(None, Linear, Weak.Quad, Strong.Quad)
model.df$Type <- factor(model.df$Type, levels=c("No association", "Linear", "Weak quadratic", "Strong quadratic"))

model.plot <- ggplot(data=model.df) +
  geom_line(aes(x=x, y=model, color=Type)) +
  scale_y_log10() +
  scale_color_manual(values=cbbPalette1) +
  theme(legend.position="none") + 
  labs(x="Mite load", y="Gregarine load", color="")
  

plot(model.plot)

sim.func <- function(greg_lambda=exp(0.4488), n=1e3, b0= 9.210e-01, b1=1.068e-01, b2=-2.088e-03, theta=0.08127322236) {
  
  set.seed(1)
  
  Greg_load <- rnbinom(n=n, mu=exp(0.4488), size=theta)
  log_lambda <- b0 + b1*Greg_load + b2*Greg_load^2
  
  Mite_load <- rnbinom(n=n, mu=exp(log_lambda), size=theta)
  
  return( data.frame(Greg_load=Greg_load, Mite_load=Mite_load, Greg_Infected=as.numeric(Greg_load>0), Mite_Infected=as.numeric(Mite_load>0)) )
  
}

#sim.dat <- data.frame(
#  Greg_Infected = rbinom(n=1000, size=1, prob=0.4)
#)

sim.dat.NoCoi <- sim.func(b1=0, b2=0)

sim.means.NoCoi <- replicate(1000, {
  sim.dat.NoCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
                Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="No association")

sim.means.NoCoi$Observed <- sim.dat.NoCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()
  
print(sim.means.NoCoi)

#######

sim.dat.LinCoi <- sim.func(b1=0.25, b2=0)

sim.means.LinCoi <- replicate(1000, {
  sim.dat.LinCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="Linear association")

sim.means.LinCoi$Observed <- sim.dat.LinCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()

print(sim.means.LinCoi)


#######weak quadratic

sim.dat.WeakQuadCoi <- sim.func(b1=0.25, b2=-0.0005)

sim.means.WeakQuadCoi <- replicate(1000, {
  sim.dat.WeakQuadCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="Weak quadratic association")

sim.means.WeakQuadCoi$Observed <- sim.dat.WeakQuadCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()

print(sim.means.WeakQuadCoi)

#######strong quadratic

sim.dat.QuadCoi <- sim.func(b1=0.25, b2=-0.01)

sim.means.QuadCoi <- replicate(1000, {
  sim.dat.QuadCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="Strong quadratic association")

sim.means.QuadCoi$Observed <- sim.dat.QuadCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()

print(sim.means.QuadCoi)

sim.dat <- rbind(sim.means.NoCoi, sim.means.LinCoi, sim.means.WeakQuadCoi, sim.means.QuadCoi)

sim.dat$Species <- factor(sim.dat$Species, levels=c("No association", "Linear association", "Weak quadratic association", "Strong quadratic association"))

cscore.plot <- ggplot(data=sim.dat) +
  geom_errorbar(aes(x=Species, ymax=upper, ymin=lower, col=Species), width=0) +
  geom_point(aes(x=Species, y=Observed, col=Species), shape=17, size=3) +
  geom_point(aes(x=Species, y=mean, col=Species), size=3) +
  scale_color_manual(values=cbbPalette1) +
  labs(x="", y="Checkerboard score") + 
  theme(legend.position="right") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plot(cscore.plot)

#combined.plot <- plot_grid(model.plot, cscore.plot)
combined.plot <- ggarrange(model.plot, cscore.plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

save_plot(combined.plot, file="./Figures/CScoreSimulations.pdf", base_asp=2.2)




#####################################
###SEM: Structural equation models###
#####################################

coinf.data.2017<-read.csv("./Data/Coinfection_data_2017.csv",stringsAsFactors = T) #coinfection data (subset of mite data)
lake.data.2017<-read.csv("./Data/SummarizedData_2017.csv",
                         stringsAsFactors = T) #2017 env data

#scale environmental variables
lake.data.2017$Fish.Densitym2<-as.numeric(scale(lake.data.2017$Fish.Densitym2, center=TRUE, scale=TRUE))
lake.data.2017$pH<-as.numeric(scale(lake.data.2017$pH, center=TRUE, scale=TRUE))
lake.data.2017$Prey.DensityL<-as.numeric(scale(lake.data.2017$Prey.DensityL, center=TRUE, scale=TRUE))
lake.data.2017$Shoot.Countm2<-as.numeric(scale(lake.data.2017$Shoot.Countm2, center=TRUE, scale=TRUE))

#merge env data with parasitism data
par.data.2017<-merge(coinf.data.2017, lake.data.2017, by="Lake")

#add variable for mite scars
par.data.2017$holes.plus.mites<-par.data.2017$Total_Mites+par.data.2017$Total_Holes

#add variable for ln(mites) and ln(gregarines)
par.data.2017$ln.mites<-log(par.data.2017$holes.plus.mites+1)
par.data.2017$ln.greg<-log(par.data.2017$Total_gregarines+1)
#add variable for coinfection
par.data.2017<- par.data.2017%>%
  mutate(
    coinfected=case_when(
      Mite_Infected == 1 & Greg_Infected == 1 ~ 1,
      Mite_Infected == 0 & Greg_Infected == 1 ~ 0,
      Mite_Infected == 1 & Greg_Infected == 0 ~ 0,
      Mite_Infected == 0 & Greg_Infected == 0 ~ 0,
    )
  )

#make dataframe for SEM analyses of prevalence
par.data.2017.for.sem.prev<- par.data.2017

#create dataframe for the lake prevalence summaries for both parasite types
pro.inf.all.2017<-as.data.frame(par.data.2017 %>%
                                  dplyr::group_by(Lake,Species) %>%
                                  dplyr::summarize(n.inf.mite = sum(Mite_Infected),
                                                   n.inf.greg = sum(Greg_Infected),
                                                   n.coinf = sum(coinfected),
                                                   total = n())%>%
                                  mutate(n.not.inf.mite=total-n.inf.mite,
                                         pro.inf.mite=n.inf.mite/total,
                                         n.not.inf.greg=total-n.inf.greg,
                                         pro.inf.greg=n.inf.greg/total,
                                         n.not.coinf=total-n.coinf,
                                         pro.coinf=n.coinf/total
                                  ))
pro.inf.all.2017$Year<-as.factor("2017")

pro.inf.2017<-subset(pro.inf.all.2017,pro.inf.all.2017$total>9) #remove lake x species with n < 10

#loop to produce exact ci's for parasite prevalence
#mites
pro.inf.2017[c('lower.ci.mite', 'upper.ci.mite')] <- t(mapply(function(x, y) 
  PropCIs::exactci(x,y, conf.level=0.95)$conf.int, pro.inf.2017$n.inf.mite, pro.inf.2017$total))
#gregarines
pro.inf.2017[c('lower.ci.greg', 'upper.ci.greg')] <- t(mapply(function(x, y) 
  PropCIs::exactci(x,y, conf.level=0.95)$conf.int, pro.inf.2017$n.inf.greg, pro.inf.2017$total))

pro.inf.2017<-merge(pro.inf.2017,lake.data.2017,by="Lake")


lake.data.2019<-read.csv("./Data/SummarizedData_2019.csv",stringsAsFactors = T)
#scale environmental variables
lake.data.2019$Fish.Densitym2<-as.numeric(scale(lake.data.2019$Fish.Densitym2, center=TRUE, scale=TRUE))
lake.data.2019$pH<-as.numeric(scale(lake.data.2019$pH, center=TRUE, scale=TRUE))
lake.data.2019$Prey.DensityL<-as.numeric(scale(lake.data.2019$Prey.DensityL, center=TRUE, scale=TRUE))
lake.data.2019$Shoot.Countm2<-as.numeric(scale(lake.data.2019$Shoot.Countm2, center=TRUE, scale=TRUE))

#merge env data with parasitism data
par.data.2019<-read.csv("./Data/Parasite_Data_2019.csv",stringsAsFactors = T)
par.data.2019<-merge(par.data.2019, lake.data.2019,by="Lake")

#add variable for mite scars
par.data.2019$holes.plus.mites<-par.data.2019$final_mite+par.data.2019$total_holes

par.data.2019$holes.plus.mites
#add variable for ln(mites) and ln(greg)
par.data.2019$ln.mites<-log(par.data.2019$holes.plus.mites+1)
par.data.2019$ln.greg<-log(par.data.2019$Gregarines+1)

#add variable for coinfection
par.data.2019<- par.data.2019%>%
  mutate(
    coinfected=case_when(
      Mite_Infected == 1 & Greg_Infected == 1 ~ 1,
      Mite_Infected == 0 & Greg_Infected == 1 ~ 0,
      Mite_Infected == 1 & Greg_Infected == 0 ~ 0,
      Mite_Infected == 0 & Greg_Infected == 0 ~ 0,
    )
  )


#make dataframe for SEM analyses of prevalence
par.data.2019.for.sem.prev<- par.data.2019

#create dataframe for the lake prevalence summaries for both parasite types
pro.inf.all.2019<-as.data.frame(par.data.2019 %>%
                                  dplyr::group_by(Lake,Species) %>%
                                  dplyr::summarize(n.inf.mite = sum(Mite_Infected),
                                                   n.inf.greg = sum(Greg_Infected),
                                                   n.coinf = sum(coinfected),
                                                   total = n())%>%
                                  mutate(n.not.inf.mite=total-n.inf.mite,
                                         pro.inf.mite=n.inf.mite/total,
                                         n.not.inf.greg=total-n.inf.greg,
                                         pro.inf.greg=n.inf.greg/total,
                                         n.not.coinf=total-n.coinf,
                                         pro.coinf=n.coinf/total
                                  ))
pro.inf.all.2019$Year<-as.factor("2019")

pro.inf.2019<-subset(pro.inf.all.2019,pro.inf.all.2019$total>9) #remove lake x species with n < 10

#loop to produce exact ci's for parasite prevalence
#mites
pro.inf.2019[c('lower.ci.mite', 'upper.ci.mite')] <- t(mapply(function(x, y) 
  PropCIs::exactci(x,y, conf.level=0.95)$conf.int, pro.inf.2019$n.inf.mite, pro.inf.2019$total))
#gregarines
pro.inf.2019[c('lower.ci.greg', 'upper.ci.greg')] <- t(mapply(function(x, y) 
  PropCIs::exactci(x,y, conf.level=0.95)$conf.int, pro.inf.2019$n.inf.greg, pro.inf.2019$total))

pro.inf.2019<-merge(pro.inf.2019,lake.data.2019,by="Lake")

#create dataframe for the lake intensity summaries mites
mite.int.all.2019<-subset(par.data.2019,Mite_Infected=="1") %>%
  dplyr::group_by(Lake,Species) %>%
  dplyr::summarize(mean.int.mite = median(ln.mites),
                   sd.int.mite = sd(ln.mites),
                   n.int.mite = n()
  )%>%
  mutate(se.int.mite=sd.int.mite/sqrt(n.int.mite)
  )

#create dataframe for the lake intensity summaries gregarines
greg.int.all.2019<-subset(par.data.2019,Greg_Infected=="1") %>%
  dplyr::group_by(Lake,Species) %>%
  dplyr::summarize(mean.int.greg = median(ln.greg),
                   sd.int.greg = sd(ln.greg),
                   n.int.greg = n()
  )%>%
  mutate(se.int.greg=sd.int.greg/sqrt(n.int.greg)
  )




#############################################################
## SPLIT BY SPECIES
#############################################################

##########################################################
# ENBA
##########################################################

#first, fix colnames so they match up
colnames(par.data.2019)[5]<-"Total_Mites"
colnames(par.data.2019)[6]<-"Total_Resist"
colnames(par.data.2019)[7]<-"Total_Holes"
colnames(par.data.2019)[4]<-"Total_gregarines"

colnames(par.data.2017)[5]<-"Total_Mites"
colnames(par.data.2017)[6]<-"Total_Resist"
colnames(par.data.2017)[7]<-"Total_Holes"
colnames(par.data.2017)[4]<-"Total_gregarines"

#add variable for year
par.data.2017$Year<-"2017"
par.data.2019$Year<-"2019"


par.data<-rbind(par.data.2017[,c(1:3, 7:17)], par.data.2019[,c(1:3, 7:17)])

par.data$ln.mites2 <- par.data$ln.mites^2
par.data$Mites <- exp(par.data$ln.mites) - 1
par.data$Mites2 <- par.data$Mites^2

par.data$Greg_Infected


ggplot(data=par.data) + 
  geom_point(aes(x=exp(ln.mites)-1, y=exp(ln.greg)-1, col=Species)) +
  geom_smooth(aes(x=exp(ln.mites)-1, y=exp(ln.greg)-1, col=Species), se=F) +
  lims(x=c(0,10))

par.data.enba<-subset(par.data,Species=="ENBA")

par.data.enba$Year<-as.factor(par.data.enba$Year)

par.data.enba%>%group_by(Lake,Year)%>%tally() #all lakes have at least 10 individuals

cor(par.data.enba[,c("Fish.Densitym2","Prey.DensityL","pH","Shoot.Countm2",
                     "ln.mites","ln.greg")]) #looks good

#mites affecting greg
#ln.greg~ln.mites+pH+Prey.DensityL+Year

#env affecting mites
#ln.mites~Fish.Densitym2+Prey.DensityL+pH+Year
# add prey
model<-'
#mites affecting greg
ln.greg ~ ln.mites + ln.mites2 + Prey.DensityL + Year

#env affecting mites
ln.mites ~ pH + Prey.DensityL + Fish.Densitym2 + Shoot.Countm2 + Year
'

fit<-sem(model, data=par.data.enba, estimator="MLR")
summary(fit, fit.measures=TRUE,standardized=TRUE,
        modindices=T)


#fit<-cfa(model, data=par.data.enba, estimator="MLR")
summary(fit, fit.measures=TRUE,standardized=TRUE,
        modindices=T)

enba.data <- par.data.enba
enba.data$Fish.Densitym2 <- mean(enba.data$Fish.Densitym2)
enba.data$Shoot.Countm2 <- mean(enba.data$Shoot.Countm2)
enba.data$Prey.DensityL <- mean(enba.data$Prey.DensityL)
enba.data$pH <- mean(enba.data$pH)
enba.data$ln.mites <- seq(min(enba.data$ln.mites), max(enba.data$ln.mites), length.out=length(enba.data$ln.mites))
enba.data$ln.mites2 <- enba.data$ln.mites^2
enba.data$Year <- "2017"


enba.data$predict <- lavPredictY(fit, ynames="ln.greg", xnames=c("ln.mites2", "ln.mites"), newdata=enba.data)

# CFI = .969, RMSEA = .089, AIC= 679


my_label_list <- list(list(node = "P.D", to = "Prey\n density"),
                      list(node = "F.D", to = "Fish\n density"),
                      list(node = "S.C", to = "Shoot\n count"),
                      list(node = "Yer", to = "Year"), 
                      list(node = "ln.g", to = "Gregarines"), 
                      list(node = "ln.m", to = "Mites"),
                      list(node = "l.2", to = expression(paste("Mites"^2)) ))

p_sem <- semPaths(fit, what="std", residuals=FALSE, sizeMan=8, label.cex=1.5, label.font=7,
         edge.label.cex=1.5, exoCov = FALSE, fade=T, thresholds=F, posCol = c("blue", "red"))
p_sem <- mark_sig(p_sem, fit, alphas = c(`*` = 0.05))
p_sem_ENBA <- change_node_label(
  semPaths_plot = p_sem,
  label_list = my_label_list
)

######## ##################################################
# ENEX
##########################################################
par.data.enex<-subset(par.data,Species=="ENEX")
par.data.enex$Year<-as.factor(par.data.enex$Year)

fit<-sem(model,data=par.data.enex,estimator="MLR")
summary(fit,fit.measures=TRUE,standardized=TRUE,modindices=T)
# CFI = .991, RMSEA = .029, AIC= 2475

p_sem <- semPaths(fit,what="std",residuals=FALSE,sizeMan=8,label.cex=1.5,label.font=7,
         edge.label.cex=1.5, exoCov = FALSE, fade=T, posCol = c("blue", "red"))
p_sem <- mark_sig(p_sem, fit, alphas = c(`*` = 0.05))
p_sem_ENEX <- change_node_label(
  semPaths_plot = p_sem,
  label_list = my_label_list
)

enex.data <- par.data.enex
enex.data$Fish.Densitym2 <- mean(enex.data$Fish.Densitym2)
enex.data$Shoot.Countm2 <- mean(enex.data$Shoot.Countm2)
enex.data$Prey.DensityL <- mean(enex.data$Prey.DensityL)
enex.data$pH <- mean(enex.data$pH)
enex.data$ln.mites <- seq(min(enex.data$ln.mites), max(enex.data$ln.mites), length.out=length(enex.data$ln.mites))
enex.data$ln.mites2 <- enex.data$ln.mites^2
enex.data$Year <- "2017"


enex.data$predict <- lavPredictY(fit, ynames="ln.greg", xnames=c("ln.mites2", "ln.mites"), newdata=enex.data)


##########################################################
# ENSI
##########################################################
par.data.ensi<-subset(par.data,Species=="ENSI")
par.data.ensi$Year<-as.factor(par.data.ensi$Year)

fit<-sem(model,data=par.data.ensi)
summary(fit, fit.measures=TRUE, standardized=TRUE, modindices=T)
# CFI = 1, RMSEA = 0, AIC= 2206

p_sem <- semPaths(fit,what="std",residuals=FALSE,sizeMan=8,label.cex=1.5,label.font=7,
         edge.label.cex=1.5, exoCov = FALSE, fade=T, posCol = c("blue", "red"))
p_sem <- mark_sig(p_sem, fit, alphas = c(`*` = 0.05))
p_sem_ENSI<- change_node_label(
  semPaths_plot = p_sem,
  label_list = my_label_list
)

ensi.data <- par.data.ensi
ensi.data$Fish.Densitym2 <- mean(ensi.data$Fish.Densitym2)
ensi.data$Shoot.Countm2 <- mean(ensi.data$Shoot.Countm2)
ensi.data$Prey.DensityL <- mean(ensi.data$Prey.DensityL)
ensi.data$pH <- mean(ensi.data$pH)
ensi.data$ln.mites <- seq(min(ensi.data$ln.mites), max(ensi.data$ln.mites), length.out=length(ensi.data$ln.mites))
ensi.data$ln.mites2 <- ensi.data$ln.mites^2
ensi.data$Year <- "2017"


ensi.data$predict <- lavPredictY(fit, ynames="ln.greg", xnames=c("ln.mites2", "ln.mites"), newdata=ensi.data)


##########################################################
# ENTR
##########################################################
par.data.entr<-subset(par.data,Species=="ENTR")
par.data.entr$Year<-as.factor(par.data.entr$Year)

fit<-sem(model,data=par.data.entr)
summary(fit,fit.measures=TRUE,standardized=TRUE, modindices=T)
# CFI = 1, RMSEA = 0, AIC= 1751.842

p_sem <- semPaths(fit,what="std",residuals=FALSE,sizeMan=8,label.cex=1.5,label.font=7,
         edge.label.cex=1.5, exoCov = FALSE, fade=T, posCol = c("blue", "red"))
p_sem <- mark_sig(p_sem, fit, alphas = c(`*` = 0.05))
p_sem_ENTR <- change_node_label(
  semPaths_plot = p_sem,
  label_list = my_label_list
)


entr.data <- par.data.entr
entr.data$Fish.Densitym2 <- mean(entr.data$Fish.Densitym2)
entr.data$Shoot.Countm2 <- mean(entr.data$Shoot.Countm2)
entr.data$Prey.DensityL <- mean(entr.data$Prey.DensityL)
entr.data$pH <- mean(entr.data$pH)
entr.data$ln.mites <- seq(min(entr.data$ln.mites), max(entr.data$ln.mites), length.out=length(entr.data$ln.mites))
entr.data$ln.mites2 <- entr.data$ln.mites^2
entr.data$Year <- "2017"


entr.data$predict <- lavPredictY(fit, ynames="ln.greg", xnames=c("ln.mites2", "ln.mites"), newdata=entr.data)

##########################################################
# ISPO
##########################################################
par.data.ispo<-subset(par.data,Species=="ISPO")
par.data.ispo$Year<-as.factor(par.data.ispo$Year)

fit <- sem(model,data=par.data.ispo)
summary(fit,fit.measures=TRUE, standardized=F, modindices=T)
# CFI = .996, RMSEA = .021, AIC= 3455
 x <- seq(0, 4, length.out=100)
 
# plot(x, x*0.164 - 0.046*x^2, type='l')

ispo.data <- par.data.ispo
ispo.data$Fish.Densitym2 <- mean(ispo.data$Fish.Densitym2)
ispo.data$Shoot.Countm2 <- mean(ispo.data$Shoot.Countm2)
ispo.data$Prey.DensityL <- mean(ispo.data$Prey.DensityL)
ispo.data$pH <- mean(ispo.data$pH)
ispo.data$ln.mites <- seq(min(ispo.data$ln.mites), max(ispo.data$ln.mites), length.out=length(ispo.data$ln.mites))
ispo.data$ln.mites2 <- ispo.data$ln.mites^2
ispo.data$Year <- "2017"


ispo.data$predict <- lavPredictY(fit, ynames=lavNames(fit, "ov.y"), xnames=c("ln.mites2", "ln.mites"), newdata=ispo.data)

p_sem <- semPaths(fit, what="std", residuals=FALSE, sizeMan=8, label.cex=1.5, label.font=7, edge.label.cex=1.5, exoCov = FALSE, fade=F, posCol = c("blue", "red"))
p_sem <- mark_sig(p_sem, fit, alphas = c(`*` = 0.05))

p_sem_ISPO <- change_node_label(
  semPaths_plot = p_sem,
  label_list = my_label_list
)

pdf(file="./Figures/SEM_plot.pdf", width=21, height=14)
par(mfrow=c(2, 3))
plot(p_sem_ENBA)
plot(p_sem_ENEX)
plot(p_sem_ENSI)
plot(p_sem_ENTR)
plot(p_sem_ISPO)

dev.off()
predict.dat <- rbind(enba.data, enex.data, ensi.data, entr.data, ispo.data)

p.predict <- ggplot(data=predict.dat) +
  geom_line(aes(x=exp(ln.mites)-1, y=exp(predict)-1, col=Species)) +
  #scale_x_log10() +
  #scale_y_log10() +
  labs(x="Mites", y="Gregarines") +
  lims(x=c(0,25), y=c(0,1.5)) +
  scale_color_paletteer_d("lisa::FridaKahlo") 

plot(p.predict)
save_plot(filename="./Figures/PredictGregs.pdf", p.predict)
stop()
##log-effects
#ENBA: ln.m -0.07 , ln.m2 +0.14
#ENEX: ln.m -0.04 , ln.m2 +0.01
#ENSI: ln.m  0.00 , ln.m2 -0.03 
#ENTR: ln.m -0.06 , ln.m2  0.00
#ISPO: ln.m +0.21*, ln.m2 -0.16*


##linear effects
#ENBA: Mites 0.18 , Mites2 -0.13 
#ENEX: Mites -0.03 , Mites2 +0.02 
#ENSI: Mites  0.00 , Mites2 -0.03 
#ENTR: Mites -0.08*, Mites2  0.04
#ISPO: Mites +0.13*, Mites2 -0.12*


#TMBmodel
#ENBA: Mites  0.0357350 , Mites2 -0.0001525
#ENEX: Mites  7.822e-03 , Mites2 -5.871e-05
#ENSI: Mites  0.103019*, Mites2 -0.002438 
#ENTR: Mites  0.015189, Mites2  -0.002185
#ISPO: Mites  1.068e-01*, -2.088e-03
