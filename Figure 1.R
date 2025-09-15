##using the Hasik et al data to explore the variation in infection by 
#environmental factors and number/presence of mites
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(ggpubr)
library(MuMIn)
library(viridis)
library(mgcv)
library(MuMIn)
library(ggpubfigs)
library(performance)
theme_set(theme_pubr())


coi2017.dat <- read.csv(file="./Data/Coinfection_data_2017.csv")
env.dat <- read.csv(file="./Data/SummarizedData_2017.csv")
coi2017.dat <- left_join(coi2017.dat, env.dat, by="Lake")
coi2017.dat$Year <- "2017"

coi2019.dat <- read.csv(file="./Data/Parasite_Data_2019.csv")
coi2019.dat$Year <- "2019"
env.dat <- read.csv(file="./Data/SummarizedData_2019.csv")
coi2019.dat <- left_join(coi2019.dat, env.dat, by="Lake")


#rename stuff for consistency between the two datasets
coi2017.dat <- coi2017.dat %>% rename(total_holes  = "Total_Holes", Gregarines = "Total_gregarines")
coi2019.dat <- coi2019.dat %>% rename(Total_Mites="final_mite")


coi.dat <- rbind(coi2017.dat, coi2019.dat) 

coi.dat <- coi.dat %>% mutate(Mites=total_holes + Total_Mites)

coi.dat$Species <- as.factor(coi.dat$Species)
coi.dat$Lake    <- as.factor(coi.dat$Lake)
coi.dat$Year    <- as.factor(coi.dat$Year)

coi.dat$COI_Infected <- coi.dat$Greg_Infected*coi.dat$Mite_Infected

nSpecies <- nlevels(coi.dat$Species)

ggplot(data=coi.dat, aes(x=total_holes, y=Gregarines, col=Species)) +
  geom_point() +
  geom_smooth()

quadenv.bin.mod <- env.bin.mod <- linenv.bin.mod <- gamm.mod <- list()
quadenv.nbin.mod <- env.nbin.mod <- linenv.nbin.mod <- list()
quadenv.zip.mod <- env.zip.mod <- linenv.zip.mod <- list()

quadenv.nbin.pred <- NULL

for(i in 1:nlevels(coi.dat$Species)) {
  #this is our analysis
  env.nbin.mod[[i]] <- glmmTMB(Gregarines ~ Shoot.Countm2 + Prey.DensityL + pH + Year + (1|Lake), data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), dispformula = ~ 1,  family=nbinom2, REML=F)
  
  linenv.nbin.mod[[i]] <- glmmTMB(Gregarines ~ Mites + Shoot.Countm2 + Prey.DensityL + pH + Year + (1|Lake), data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), dispformula = ~  1, family=nbinom2, REML=F)
  
  quadenv.nbin.mod[[i]] <- glmmTMB(Gregarines ~ Mites + I(Mites^2) + Shoot.Countm2 + Prey.DensityL + pH + Year + (1|Lake), dispformula = ~  1,  data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) #AICc selects nonlinear effect in ENSI and ISPO
  
  #quadenv.nbin.mod[[i]] <- glmmTMB(Gregarines ~ Mites + I(Mites^2) + Year + (1|Lake), dispformula = ~  1,  data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) #AICc selects nonlinear effect in ENSI and ISPO
  r2(quadenv.nbin.mod[[i]])
  AICc.vec <- c(AICc(env.nbin.mod[[i]]), AICc(linenv.nbin.mod[[i]]), AICc(quadenv.nbin.mod[[i]]))
  
  print(AICc.vec)
  print(AICc.vec - min(AICc.vec))

  
  #compare to zip model
  env.zip.mod[[i]] <- glmmTMB(Gregarines ~ Year  + (1|Lake), ziformula=~1, dispformula = ~  1, data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) 
  
  linenv.zip.mod[[i]] <- glmmTMB(Gregarines ~ Mites + Year + (1|Lake), ziformula=~1, dispformula = ~  1, data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) 
  
  quadenv.zip.mod[[i]] <- glmmTMB(Gregarines ~ Mites + I(Mites^2) + Year + (1|Lake), ziformula=~ 1, dispformula = ~  1, data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) 
  
  AICc.zip <- c(AICc(env.zip.mod[[i]]), AICc(linenv.zip.mod[[i]]), AICc(quadenv.zip.mod[[i]]))
  
  print('-----')
}

#get predictions of best models
x <- 0:600
lin.dat <- data.frame(Mites=x, `I(Mites^2)`= x^2, Year="2017", Shoot.Countm2=mean(coi.dat$Shoot.Countm2), Prey.DensityL=mean(coi.dat$Prey.DensityL), pH=mean(coi.dat$pH), Lake="New" )

predict.dat <- rbind(data.frame(Mites=x, Gregarines=exp(predict(linenv.nbin.mod[[1]], newdata=lin.dat)), species="ENBA"),
                     data.frame(Mites=x, Gregarines=exp(predict(linenv.nbin.mod[[2]], newdata=lin.dat)), species="ENEX"), 
                     data.frame(Mites=x, Gregarines=exp(predict(env.nbin.mod[[3]], newdata=lin.dat)), species="ENSI"), 
                     data.frame(Mites=x, Gregarines=exp(predict(linenv.nbin.mod[[4]], newdata=lin.dat)), species="ENTR"),
                     data.frame(Mites=x, Gregarines=exp(predict(quadenv.nbin.mod[[5]], newdata=lin.dat)), species="ISPO"))


p.predict <- ggplot(data=coi.dat) +
  geom_point(data=coi.dat, aes(x=Mites, y=Gregarines, col=Species), alpha=0.25) +
  geom_line(data=predict.dat, aes(x=Mites, y=Gregarines, col=species)) +
  scale_y_sqrt() +
  scale_x_sqrt() +
  #lims(x=c(0, 60), y=c(c(0, 10))) +
  labs(x="Mite load", y="Gregarine load", col="Species") +
  scale_color_manual(values = friendly_pal("ibm_five"))
plot(p.predict)
save_plot(p.predict, filename="./Figures/GLM_predictions.pdf")


species1 <- data.frame(as.data.frame(ggpredict(linenv.nbin.mod[[1]], terms=c("Mites", "Year"))), Species="ENBA")
temp <- as.data.frame(ggpredict(env.nbin.mod[[1]], terms="Year"))
species1$predicted[species1$group == "2017"] <- temp$predicted[1]
species1$predicted[species1$group == "2019"] <- temp$predicted[2]

species2 <- data.frame(as.data.frame(ggpredict(linenv.zip.mod[[2]], terms=c("Mites", "Year"))), Species="ENEX")
species3 <- data.frame(as.data.frame(ggpredict(linenv.nbin.mod[[3]], terms=c("Mites", "Year"))), Species="ENSI")
temp <- as.data.frame(ggpredict(env.nbin.mod[[3]], terms="Year"))
species3$predicted[species3$group == "2017"] <- temp$predicted[1]
species3$predicted[species3$group == "2019"] <- temp$predicted[2]


species4 <- data.frame(as.data.frame(ggpredict(linenv.zip.mod[[4]], terms=c("Mites", "Year"))), Species="ENTR")
species5 <- data.frame(as.data.frame(ggpredict(quadenv.zip.mod[[5]], terms=c("Mites", "Year"))), Species="ISPO")

df <- rbind(species1, species2, species3, species4, species5)
df <- df %>% rename(Mites=x, Gregarines=predicted, Year=group)

p.predict <- ggplot(data=coi.dat) +
  geom_point(data=coi.dat, aes(x=Mites, y=Gregarines, col=Species, shape=Year), alpha=0.25) +
  geom_line(data=df, aes(x=Mites, y=Gregarines, col=Species, linetype=Year), linewidth=0.75) +
  scale_y_sqrt() +
  scale_x_sqrt() +
  #lims(x=c(0, 60), y=c(c(0, 10))) +
  labs(x="Mite load", y="Gregarine load", col="Species") +
  scale_color_manual(values = friendly_pal("ibm_five")) +
  facet_wrap(~ Species) 
plot(p.predict)
  save_plot(p.predict, filename="./Figures/GLM_predictions_species.pdf", base_aspect_ratio=1.6, base_height = 6)
