##using the Hasik et al data to explore the variation in infection by 
#environmental factors and number/presence of mites
library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(ggpubr)
library(MuMIn)
library(viridis)
library(mgcv)

# Hasik: look at within population analyses where data are suitbale...
# 

#model by species. Incorporate all environmental data.
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

coi.dat$holesplusmites <- coi.dat$Total_Mites + coi.dat$total_holes
#coi.dat <- coi.dat %>% filter(Species != "ENEX") #%>% filter(Species != "ENEX", Species != "ENTR") #very littlle coinfection for this species

#coi.dat$total_holes <- log(coi.dat$total_holes + 1)
coi.dat$Species <- as.factor(coi.dat$Species)
coi.dat$Lake <- as.factor(coi.dat$Lake)
coi.dat$Year <- as.factor(coi.dat$Year)

coi.dat$COI_Infected <- coi.dat$Greg_Infected*coi.dat$Mite_Infected

nSpecies <- nlevels(coi.dat$Species)

ggplot(data=coi.dat, aes(x=total_holes, y=Gregarines, col=Species)) +
  geom_point() +
  geom_smooth()

quadenv.bin.mod <- list()
quadenv.bin.pred <- NULL
quadenv.nbin.mod <- quadenv.nbin.mod2 <- list()
    quadenv.zip.mod <-quadenv.zip.mod2 <- quadenv.bin.mod2 <- list()
    
    coi.bin.mod <- list()
    coi.bin.mod2 <- list()
    coi.bin.mod.hasik <- list()
    gam.bin.mod <-gam.nbin.mod <- quadenv.nbin.mod0 <- list()
    
    greg.bin.mod <- list()
    greg.bin.mod.hasik <- list()
    
    mite.bin.mod <- list()
    mite.bin.mod.hasik <- list()
    
    quadenv.nbin.pred <- NULL
    
    for(i in 1:nlevels(coi.dat$Species)) {
      #this is our analysis
      
      #count models
      #environmental model
      quadenv.nbin.mod0[[i]] <- glmmTMB(Gregarines ~ Fish.Densitym2 + pH + Prey.DensityL + Shoot.Countm2 + Year + (1|Lake), data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), dispformula = ~ 1, family=nbinom2, REML=F)
      
      #linear mite model
      quadenv.nbin.mod[[i]] <- glmmTMB(Gregarines ~ Fish.Densitym2  + pH + Prey.DensityL + Shoot.Countm2 + Year + holesplusmites + (1|Lake), dispformula = ~  1,  data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) #AICc selects nonlinear effect in ENSI and ISPO
      
      
      quadenv.nbin.mod2[[i]] <- glmmTMB(Gregarines ~ Fish.Densitym2 + pH + Prey.DensityL + Shoot.Countm2 + Year + holesplusmites + I(holesplusmites^2) + (1|Lake), data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), dispformula = ~ 1, family=nbinom2, REML=F)
      

      
      print(c(AICc(quadenv.nbin.mod0[[i]]), AICc(quadenv.nbin.mod[[i]]), AICc(quadenv.nbin.mod2[[i]])))
      
      print(c(AICc(quadenv.nbin.mod0[[i]]), AICc(quadenv.nbin.mod[[i]]), AICc(quadenv.nbin.mod2[[i]])) - min(c(AICc(quadenv.nbin.mod0[[i]]), AICc(quadenv.nbin.mod[[i]]), AICc(quadenv.nbin.mod2[[i]]))))
      
      #add zip model
      quadenv.zip.mod[[i]] <- glmmTMB(Gregarines ~ pH + holesplusmites + I(holesplusmites^2) + Year + (1|Lake), ziformula=~  pH + Year + (1|Lake), data=filter(coi.dat, Species==levels(coi.dat$Species)[i]), family=nbinom2, REML=F) #nonlinear effect in 2 (positive)
      
  temp <- ggpredict(quadenv.nbin.mod[[i]], c("holesplusmites [all]"))
  quadenv.nbin.pred <- rbind(quadenv.nbin.pred, data.frame(Probability=temp$predicted, holesplusmites=temp$x, Species=levels(coi.dat$Species)[i], conf.low=temp$conf.low, conf.hi=temp$conf.high))

}

ggplot(quadenv.bin.pred, aes(x=holesplusmites, y=Probability, col=Species)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.hi, fill=Species, col=NULL), alpha=0.2) +
  geom_line() +
  geom_point(data=coi.dat, aes(x=holesplusmites, y=Greg_Infected, col=Species)) +
  labs(x="Number of mites", y="Probability of gregarine", col="Species", fill="Species") + 
  scale_color_viridis(discrete=TRUE, option="H") +
  scale_fill_viridis(discrete=TRUE, option="H") +
  scale_x_sqrt()


ggplot(quadenv.nbin.pred, aes(x=holesplusmites, y=Probability, col=Species)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.hi, fill=Species, col=NULL), alpha=0.2) +
  geom_line() +
  #geom_jitter(data=coi.dat, aes(x=holesplusmites, y=Gregarines, col=Species), height=0, width=2) +
  labs(x="Number of mites", y="Number of gregarines", col="Species", fill="Species") + 
  scale_color_viridis(discrete=TRUE, option="H") +
  scale_fill_viridis(discrete=TRUE, option="H") +
  coord_cartesian(ylim=c(0, 10), xlim=c(0, 100)) #+
  #scale_x_sqrt() #+
  #scale_y_sqrt()
# +
  #coord_cartesian(ylim=c(0, 100))


#reproduce models from Hasik
stop()

names(quadenv.bin.mod) <- levels(coi.dat$Species)
names(quadenv.nbin.mod) <- levels(coi.dat$Species)


null.mod <- glmmTMB(Greg_Infected ~ Species + Year + (1|Lake),
                               data=coi.dat,
                               family=binomial, REML=F)

env.mod <- glmmTMB(Greg_Infected ~ pH*Species + Year + (1|Lake),
                    data=coi.dat,
                    family=binomial, REML=F)

lin.mod <- glmmTMB(Greg_Infected ~ holesplusmites*Species + Year + (1|Lake),
       data=coi.dat,
       family=binomial, REML=F)

quad.mod <- glmmTMB(Greg_Infected ~ holesplusmites*Species + I(holesplusmites^2)*Species + Year + (1|Lake),
                   data=coi.dat, family=binomial, REML=F)


linenv.mod <- glmmTMB(Greg_Infected ~ holesplusmites*Species + pH*Species + Year + (1|Lake),
                   data=coi.dat,
                   family=binomial, REML=F)

quadenv.mod <- glmmTMB(Greg_Infected ~ holesplusmites*Species + pH*Species + I(holesplusmites^2)*Species + Year + (1|Lake),
                    data=coi.dat, family=binomial, REML=F)

#A more flexible non-parametric model
library(mgcv)
#quadenv.gam <- gamm(Greg_Infected ~ s(total_holes, by=Species) + s(pH, by=Species), data=coi.dat, family=binomial, random=list(Lake=~1))


AICc.vec <- c(AICc(null.mod), AICc(env.mod), AICc(lin.mod), AICc(quad.mod), AICc(linenv.mod), AICc(quadenv.mod))
#AICc.vec: 3131.216 3092.900 3126.869 3119.040 3093.605 3084.402 total_holes #quadratic model is best
#AICc.vec: 3131.216 3092.900 3130.949 3134.535 3093.668 3097.450 Total_Mites #fit is a bit worse. Linear model is best
(AICc.vec - min(AICc.vec))

pr1 <- ggpredict(quadenv.mod, c("holesplusmites", "Species"))


